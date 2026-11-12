"""Check sector composition, physical transforms, time alignment and resumption."""

from pathlib import Path

import h5py
import meshio
import numpy as np
import pytest
from numpy.typing import ArrayLike

from porous_media.data.xdmf_tools import (
    XDMFTimeSeriesReader,
    XDMFTimeSeriesWriter,
)
from research.geometry_conversion import sector_transforms
from research.spt.cache import ANALYSIS_FIELDS, load_processed_run
from research.spt.catalog import discover_runs
from research.spt.conditions import (
    MIXED_EXCLUDED_PATTERN,
    pattern_name,
    simulation_conditions,
)
from research.spt.data_processing import main, process_studies
from research.spt.fields import GEOMETRY_CELL_FIELDS
from research.spt.geometry import (
    CONVERGENCE_MESHES,
    SECTORS,
    geometry_jobs,
    geometry_path,
    lobulus_sources,
    process_geometries,
    process_geometry,
    write_composition,
)
from research.spt.storage import COMPLEVEL


def _source(root: Path, name: str, count: int, times: tuple[float, ...]) -> Path:
    """Write a physically valid extruded sixty-degree fan with directional fields."""
    study = (
        "convergence_sixth" if name.startswith("convergence") else "simulations_sixth"
    )
    path = root / study / "xdmf" / f"{name}.xdmf"
    path.parent.mkdir(parents=True, exist_ok=True)
    edges = np.linspace([0.0, 1.0, 0.0], [-np.sqrt(3) / 2, 0.5, 0.0], count + 1)
    bottom = np.vstack(([0.0, 0.0, 0.0], edges))
    points = np.vstack((bottom, bottom + np.array([0, 0, 0.2])))
    lower = np.array([[0, index + 1, index + 2] for index in range(count)])
    connectivity = np.column_stack((lower, lower + len(bottom)))
    centers = points[connectivity].mean(axis=1)
    with XDMFTimeSeriesWriter(path) as writer:
        writer.write_points_cells(points, [("wedge", connectivity)])
        for time in times:
            data: dict[str, list[ArrayLike]] = {
                field: [np.full((count, 1), count + 2 * time, dtype=np.float32)]
                for field in ANALYSIS_FIELDS
            }
            data["element_volume_point_TPM"] = [
                np.full((count, 1), np.sqrt(3) * 0.05 / count)
            ]
            data["rr_necrosis"] = [np.full((count, 1), time >= 2, dtype=np.float32)]
            data["rr_zonation_pattern"] = [np.full((count, 1), 2, dtype=np.float32)]
            data["cell_type"] = [np.full((count, 1), int(time), dtype=np.int32)]
            data["flux"] = [np.tile([1.0, 2.0, 3.0], (count, 1))]
            data["stress"] = [
                np.tile(
                    [[2.0, 1.0, 0.0], [1.0, 4.0, 0.0], [0.0, 0.0, 8.0]], (count, 1, 1)
                )
            ]
            for index, coordinate in enumerate("xyz"):
                data[f"position_{coordinate}"] = [centers[:, index : index + 1]]
            writer.write_data(
                time,
                point_data={"displacement": np.tile([1.0, 0.0, 0.0], (len(points), 1))},
                cell_data=data,
            )
    material = (
        root
        / study
        / "material"
        / (
            f"{name}.feb"
            if study == "convergence_sixth"
            else "lobule_zonation_pattern.feb.template"
        )
    )
    material.parent.mkdir(parents=True, exist_ok=True)
    material.write_text(
        "<febio_spec><Control><time_steps>4</time_steps><step_size>1</step_size></Control></febio_spec>"
    )
    return path


@pytest.fixture
def convergence(tmp_path: Path) -> tuple[Path, Path]:
    """Provide unequal meshes with different output grids and truncated histories."""
    root = tmp_path / "raw"
    for index, key in enumerate((*CONVERGENCE_MESHES, "50000")):
        _source(
            root,
            f"convergence_{key}",
            index + 1,
            (0.0, 1.0, 2.0) if index == 0 else (0.0, 1.5, 3.0),
        )
    return root, tmp_path / "derived"


def _determinants(points: np.ndarray, cells: list[meshio.CellBlock]) -> np.ndarray:
    values = np.concatenate([points[block.data] for block in cells])
    return np.linalg.det(
        np.stack(
            [
                values[:, 1] - values[:, 0],
                values[:, 2] - values[:, 0],
                values[:, 3] - values[:, 0],
            ],
            axis=-1,
        )
    )


def test_mixed_mesh_preserves_physics_and_discrete_states(
    convergence: tuple[Path, Path],
) -> None:
    """Rotate vectors/tensors, preserve orientation, and interpolate only continuous data."""
    root, target = convergence
    job = geometry_jobs(root, study="convergence", lobulus=False)[0]
    path = target / "mixed.xdmf"
    metadata = write_composition(job.sources, path, all_fields=True)
    assert metadata["tend_s"] == 2
    assert metadata["num_cells"] == 21
    with XDMFTimeSeriesReader(path) as reader:
        points, cells = reader.read_points_cells()
        assert np.all(_determinants(points, cells) > 0)
        assert reader.num_steps == 3
        time, point_data, cell_data = reader.read_data(1)
    assert time == 1
    sectors = np.concatenate(cell_data["sector_id"])
    matrices = sector_transforms(2 * np.pi / 3)
    tensor = np.array([[2, 1, 0], [1, 4, 0], [0, 0, 8]])
    flux = np.concatenate(cell_data["flux"])
    stress = np.concatenate(cell_data["stress"])
    centers = np.concatenate([points[block.data].mean(axis=1) for block in cells])
    for index, matrix in enumerate(matrices):
        selected = sectors == index
        np.testing.assert_allclose(
            flux[selected], np.tile(matrix @ [1, 2, 3], (index + 1, 1)), atol=1e-6
        )
        np.testing.assert_allclose(
            stress[selected],
            np.tile(matrix @ tensor @ matrix.T, (index + 1, 1, 1)),
            atol=1e-6,
        )
        np.testing.assert_allclose(np.trace(stress[selected], axis1=1, axis2=2), 14)
        np.testing.assert_allclose(
            np.concatenate(cell_data["rr_(T)"])[selected], index + 3
        )
        np.testing.assert_array_equal(
            np.concatenate(cell_data["cell_type"])[selected], 1 if index == 0 else 0
        )
        np.testing.assert_array_equal(
            np.concatenate(cell_data["rr_necrosis"])[selected], 0
        )
    for index, coordinate in enumerate("xyz"):
        np.testing.assert_allclose(
            np.concatenate(cell_data[f"position_{coordinate}"]).reshape(-1),
            centers[:, index],
            atol=1e-6,
        )
    assert point_data["displacement"].shape[1] == 3
    np.testing.assert_allclose(
        np.concatenate(cell_data["element_volume_point_TPM"]).sum(),
        6 * np.sqrt(3) * 0.05,
    )
    # The pair remains readable after relocation without access to the source tree.
    moved = target / "moved"
    moved.mkdir()
    path.rename(moved / path.name)
    path.with_suffix(".h5").rename(moved / path.with_suffix(".h5").name)
    with XDMFTimeSeriesReader(moved / path.name) as reader:
        reader.read_points_cells()
        reader.read_data(2)


def test_lobulus_pipeline_and_numeric_analysis(
    convergence: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Analyze full lobules with six times the volume and unchanged weighted means."""
    root, target = convergence
    results = process_geometries(
        root, target, study="convergence", combinations=False, workers=1
    )
    assert len(results) == 7
    assert all(result.status == "processed" for result in results)
    numerical = process_studies(
        root, target, study="convergence", geometry="lobulus", workers=1
    )
    assert all(result.status == "processed" for result in numerical)
    source = lobulus_sources(root, target, study="convergence")[0]
    assert source.study == "convergence_lobulus"
    cells, points, summary = load_processed_run(source, target)
    try:
        assert cells.sizes["cell"] == 6
        np.testing.assert_allclose(
            cells["element_volume_point_TPM"].sum("cell"), 6 * np.sqrt(3) * 0.05
        )
        np.testing.assert_allclose(summary["rr_(T)_mean"], [1, 3, 5])
    finally:
        cells.close()
        points.close()
        summary.close()
    monkeypatch.setattr(
        "research.spt.geometry.inspect_run",
        lambda *args: pytest.fail(
            "Current geometry must not reopen descriptors or numerical data"
        ),
    )
    cached = process_geometries(
        root, target, study="convergence", combinations=False, workers=1
    )
    assert all(result.status == "cached" for result in cached)


def test_failed_rebuild_preserves_geometry(
    convergence: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep the prior generation published if reconstruction or figure creation fails."""
    root, target = convergence
    job = geometry_jobs(
        root, study="convergence", combinations=False, run_ids=("convergence_06250",)
    )[0]
    assert process_geometry(job, target).status == "processed"
    path = geometry_path(job.directory(target))
    before = (job.directory(target) / "manifest.json").read_bytes()

    def fail(*args: object) -> None:
        raise RuntimeError("Synthetic interruption")

    monkeypatch.setattr("research.spt.geometry.write_composition", fail)
    assert process_geometry(job, target, force=True).status == "failed"
    assert (job.directory(target) / "manifest.json").read_bytes() == before
    assert geometry_path(job.directory(target)) == path
    assert not list(job.directory(target).glob(".partial-*"))


def test_comparison_figures_and_hdf_storage(
    convergence: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Render comparisons headlessly and store repeated fields only once per sidecar."""
    root, target = convergence
    result = process_geometries(
        root, target, study="convergence", lobulus=False, workers=1
    )[0]
    assert result.status == "processed", result.message
    path = Path(str(result.row["xdmf"]))
    assert path.with_suffix(".png").stat().st_size > 1000
    # The comparison preview is a raster figure; its vector copy was 10 times the
    # size of the PNG and nothing reads it.
    assert not path.with_suffix(".pdf").exists()
    monkeypatch.chdir(target.parent)
    assert (
        geometry_path(Path("derived/meshes/convergence_combinations/convergence_mixed"))
        == path
    )
    with h5py.File(path.with_suffix(".h5"), "r") as handle:
        # Each attribute is one dataset spanning the series rather than one
        # array per step, so a chunk covers consecutive steps.
        stacked = handle["Cell"]
        assert isinstance(stacked, h5py.Group)
        assert set(ANALYSIS_FIELDS) <= set(stacked)
        steps = int(str(result.row["num_steps"]))
        for name in stacked:
            dataset = stacked[name]
            assert isinstance(dataset, h5py.Dataset)
            assert dataset.shape[0] == steps
            assert dataset.compression == "gzip"
            assert dataset.chunks is not None and dataset.chunks[0] > 1


def test_pattern_selection_and_cli_dry_run(tmp_path: Path) -> None:
    """Group patterns by archived dose identity rather than adjacent file order."""
    root, target = tmp_path / "raw", tmp_path / "derived"
    conditions = simulation_conditions()
    for condition in conditions:
        _source(root, condition.sim_id, 1, (0.0, 4.0))
    doses = len({condition.boundary_flow_key for condition in conditions})
    jobs = geometry_jobs(root, study="simulations", lobulus=False)
    assert len(jobs) == doses
    # One sector per pattern: the comparison shows every pattern beside the
    # measured CYP2E1 profile, which keeps its own lobule and figure row.
    assert [source.run_id for source in jobs[0].sources] == [
        f"sim{1 + pattern * doses:03}" for pattern in range(SECTORS)
    ]
    compared = [
        source.condition.pattern_key for source in jobs[-1].sources if source.condition
    ]
    assert compared == list(range(SECTORS))
    excluded = {
        condition.pattern_key
        for condition in conditions
        if pattern_name(condition.pattern_key) == MIXED_EXCLUDED_PATTERN
    }
    assert excluded and not excluded.intersection(compared)
    assert (
        main(
            [
                "--input-root",
                str(root),
                "--output-root",
                str(target),
                "--study",
                "simulations",
                "--step",
                "geometry",
                "--dry-run",
            ]
        )
        == 0
    )
    assert not target.exists()
    (root / f"simulations_sixth/xdmf/{conditions[-1].sim_id}.xdmf").unlink()
    with pytest.raises(ValueError, match="Missing expected"):
        geometry_jobs(root, study="simulations", lobulus=False)


def test_geometry_rejects_wrong_sector_and_output_links(
    convergence: tuple[Path, Path],
) -> None:
    """Reject unsafe output directories and non-sector geometry before publication."""
    root, target = convergence
    runs = discover_runs(root, ("convergence_sixth",))
    job = geometry_jobs(root, study="convergence", combinations=False)[0]
    target.mkdir()
    (target / "meshes").symlink_to(root / "convergence_sixth", target_is_directory=True)
    assert process_geometry(job, target).status == "failed"
    with pytest.raises(ValueError, match="requires six"):
        write_composition(
            tuple(runs[:5]), tmp_path_placeholder := target / "invalid.xdmf"
        )
    assert not tmp_path_placeholder.exists()
    unsafe = root / "convergence_sixth" / "new.xdmf"
    with pytest.raises(ValueError, match="outside the input study"):
        write_composition((runs[0],) * 6, unsafe)
    assert not unsafe.exists()
    # A malformed source angle must not be accepted as a sixth of a lobule.
    with h5py.File(runs[0].xdmf.with_suffix(".h5"), "r+") as handle:
        coordinates = handle["data0"]
        values = coordinates[()]
        values[:, 1] *= 0.5
        coordinates[:] = values
    with pytest.raises(ValueError, match="60-degree"):
        write_composition((runs[0],) * 6, target / "wrong-sector.xdmf")


def test_geometry_spawn_and_source_invalidation(convergence: tuple[Path, Path]) -> None:
    """Run bounded spawned workers and reject derived geometries after a raw update."""
    root, target = convergence
    run_ids = ("convergence_06250", "convergence_09375")
    results = process_geometries(
        root,
        target,
        study="convergence",
        run_ids=run_ids,
        combinations=False,
        workers=2,
    )
    assert all(result.status == "processed" for result in results)
    raw = root / "convergence_sixth/xdmf/convergence_06250.h5"
    raw.touch()
    with pytest.raises(ValueError, match="Missing or stale"):
        lobulus_sources(root, target, study="convergence", run_ids=run_ids)
    status = process_geometries(
        root,
        target,
        study="convergence",
        run_ids=run_ids,
        combinations=False,
        dry_run=True,
    )
    assert [result.status for result in status] == ["pending", "cached"]


def test_complete_cli_analyzes_lobules_only(convergence: tuple[Path, Path]) -> None:
    """Run the user workflow and exclude visualization mixtures from numerical caches."""
    root, target = convergence
    assert (
        main(
            [
                "--input-root",
                str(root),
                "--output-root",
                str(target),
                "--study",
                "convergence",
                "--step",
                "all",
                "--workers",
                "1",
            ]
        )
        == 0
    )
    assert (
        len(list((target / "cache/cells/convergence_lobulus").glob("*/manifest.json")))
        == 7
    )
    assert not (target / "cache/cells/convergence_sixth").exists()
    assert not (target / "cache/cells/convergence_combinations").exists()
    assert "figure_png" in (target / "logs/geometry-runs.tsv").read_text()


def test_interrupted_staging_preserves_published_pair(
    convergence: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Keep the previous geometry readable and remove temporary copies on I/O failure."""
    root, target = convergence
    job = geometry_jobs(
        root, study="convergence", combinations=False, run_ids=("convergence_06250",)
    )[0]
    assert process_geometry(job, target).status == "processed"
    descriptor = geometry_path(job.directory(target))

    def fail_copy(source: Path, destination: Path) -> None:
        raise OSError("Synthetic network interruption")

    monkeypatch.setattr("research.spt.geometry.shutil.copyfile", fail_copy)
    failure = process_geometry(job, target, force=True)
    assert failure.status == "failed"
    assert "network interruption" in failure.message
    assert geometry_path(job.directory(target)) == descriptor
    assert not list(target.rglob(".sources-*"))
    assert not list(target.rglob(".partial-*"))
    with XDMFTimeSeriesReader(descriptor) as reader:
        reader.read_points_cells()
        reader.read_data(0)


def test_preview_style_refresh_reuses_simulation_arrays(
    convergence: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Refresh mixed previews atomically without reconstructing large HDF files."""
    root, target = convergence
    job = geometry_jobs(root, study="convergence", lobulus=False)[0]
    assert process_geometry(job, target).status == "processed"
    previous = geometry_path(job.directory(target))
    original_hdf = previous.with_suffix(".h5").stat()

    def unexpected(*args: object, **kwargs: object) -> None:
        pytest.fail("A preview style change must not reconstruct simulation arrays")

    monkeypatch.setattr("research.spt.geometry.write_composition", unexpected)
    monkeypatch.setattr("research.spt.geometry_figures.COMPARISON_PLOT_VERSION", 999)
    results = process_geometries(root, target, study="convergence", lobulus=False)
    assert [result.status for result in results] == ["cached"]
    current = geometry_path(job.directory(target))
    assert current != previous
    current_hdf = current.with_suffix(".h5").stat()
    assert (current_hdf.st_ino, current_hdf.st_size, current_hdf.st_mtime_ns) == (
        original_hdf.st_ino,
        original_hdf.st_size,
        original_hdf.st_mtime_ns,
    )
    manifest = job.directory(target) / "manifest.json"
    before = manifest.read_bytes()

    def fail(*args: object, **kwargs: object) -> None:
        raise RuntimeError("Synthetic preview failure")

    monkeypatch.setattr("research.spt.geometry_figures.COMPARISON_PLOT_VERSION", 1000)
    monkeypatch.setattr("research.spt.geometry_figures.render_comparison", fail)
    results = process_geometries(root, target, study="convergence", lobulus=False)
    assert [result.status for result in results] == ["failed"]
    assert manifest.read_bytes() == before
    assert geometry_path(job.directory(target)) == current
    assert not list(target.rglob(".partial-*"))


def test_composition_writes_only_the_fields_the_pipeline_reads(
    convergence: tuple[Path, Path],
) -> None:
    """Drop unread vectors, tensors and coordinates, and keep every analysis field."""
    root, target = convergence
    job = geometry_jobs(root, study="convergence", lobulus=False)[0]
    selected = target / "selected.xdmf"
    complete = target / "complete.xdmf"
    write_composition(job.sources, selected)
    write_composition(job.sources, complete, all_fields=True)
    with XDMFTimeSeriesReader(selected) as reader:
        reader.read_points_cells()
        _, point_data, cell_data = reader.read_data(1)
    # Every analysis field of this fixture is selected; the composition adds its
    # own sector metadata. GEOMETRY_CELL_FIELDS also names display-only fields
    # that these synthetic sources do not provide.
    assert set(ANALYSIS_FIELDS) <= set(GEOMETRY_CELL_FIELDS)
    assert set(cell_data) == set(ANALYSIS_FIELDS) | {"sector_id", "source_cell_id"}
    assert set(point_data) == {"sector_id"}
    for dropped in ("stress", "flux", "position_x", "position_y", "position_z"):
        assert dropped not in cell_data
    assert "displacement" not in point_data
    # Every analysis field survives, so the downstream stages still have inputs.
    for required in ANALYSIS_FIELDS:
        assert required in cell_data
    assert selected.with_suffix(".h5").stat().st_size < (
        complete.with_suffix(".h5").stat().st_size
    )


def test_composition_compresses_every_written_array(
    convergence: tuple[Path, Path],
) -> None:
    """Keep gzip and the shuffle filter on each sidecar dataset, as the manifest claims."""
    root, target = convergence
    job = geometry_jobs(root, study="convergence", lobulus=False)[0]
    path = target / "compressed.xdmf"
    write_composition(job.sources, path)
    datasets: list[h5py.Dataset] = []
    with h5py.File(path.with_suffix(".h5"), "r") as handle:
        # Stacked attributes live in a group per center, so the walk recurses.
        handle.visititems(
            lambda _, node: (
                datasets.append(node) if isinstance(node, h5py.Dataset) else None
            )
        )
        assert datasets
        for dataset in datasets:
            assert dataset.compression == "gzip"
            assert dataset.compression_opts == COMPLEVEL
            assert dataset.shuffle
    assert f"gzip level {COMPLEVEL}" in str(job.configuration()["compression"])
