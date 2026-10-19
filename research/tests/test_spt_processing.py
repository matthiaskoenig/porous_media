"""Exercise SPT discovery, immutable inputs, cache reuse, and interruption recovery."""

import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr

from porous_media.data.xdmf_calculations import mesh_datasets_from_xdmf
from porous_media.data.xdmf_tools import XDMFTimeSeriesWriter
from research.spt import cache
from research.spt.cache import (
    ANALYSIS_FIELDS,
    cache_status,
    load_processed_run,
    process_run,
    run_directory,
)
from research.spt.catalog import RunSource, discover_runs, inspect_run
from research.spt.conditions import conditions_from_source
from research.spt.data_processing import main, process_studies
from research.spt.information import main as information_main
from research.spt.information import simulation_conditions_df


def write_run(root: Path, study: str, run_id: str, *, truncated: bool = False) -> None:
    """Create small mixed-cell HDF fixtures with varying volumes and five frames."""
    path = root / study / "xdmf" / f"{run_id}.xdmf"
    path.parent.mkdir(parents=True, exist_ok=True)
    material = root / study / "material"
    material.mkdir(exist_ok=True)
    model = material / (
        "lobule_zonation_pattern.feb.template"
        if study == "simulations_sixth"
        else f"{run_id}.feb"
    )
    model.write_text(
        "<febio_spec><Control><time_steps>4</time_steps><step_size>1</step_size></Control></febio_spec>"
    )
    with XDMFTimeSeriesWriter(path, data_format="HDF") as writer:
        writer.write_points_cells(
            np.array(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]
            ),
            [("triangle", np.array([[0, 1, 2]])), ("quad", np.array([[0, 1, 3, 2]]))],
        )
        for step in range(3 if truncated else 5):
            values = np.array([1.0, 3.0], dtype="float32") + step
            data = {name: [values[:1], values[1:]] for name in ANALYSIS_FIELDS}
            volumes = np.array([1.0 + step, 3.0 + step], dtype="float32")
            data["element_volume_point_TPM"] = [volumes[:1], volumes[1:]]
            count = (
                np.array([1, 2], dtype="int32") if step == 0 else np.array([1.5, 2.5])
            )
            data["counter"] = [count[:1], count[1:]]
            writer.write_data(
                float(step),
                cell_data=data,
                point_data={
                    "pressure": np.arange(step, step + 4, dtype="float32"),
                    "displacement": np.zeros((4, 3)),
                },
            )


@pytest.fixture
def study_data(tmp_path: Path) -> tuple[Path, Path]:
    """Provide complete simulation settings and an independently named convergence run."""
    root = tmp_path / "source"
    write_run(root, "simulations_sixth", "sim001")
    write_run(root, "simulations_sixth", "sim002")
    write_run(root, "convergence_sixth", "convergence_06250", truncated=True)
    (root / "simulations_sixth/material/feb_manipulation_v2.py").write_text(
        'zonation_pattern = [0, 1]\nboundary_flow = [-2.5]\nraise RuntimeError("Never execute archived code")\nfor pattern in zonation_pattern:\n    for flow in boundary_flow:\n        pass\n'
    )
    return root, tmp_path / "output"


def sources(root: Path) -> list[RunSource]:
    """Discover the synthetic study through the public workflow."""
    return discover_runs(root)


def modify_hdf(source: RunSource, field: str, step: int, value: float) -> None:
    """Change one source scalar without changing its shape or time coordinates."""
    tree = ET.parse(source.xdmf)
    grids = tree.findall(".//Grid[@CollectionType='Temporal']/Grid")
    attribute = next(
        attr for attr in grids[step].findall("Attribute") if attr.get("Name") == field
    )
    item = attribute.find("DataItem")
    assert item is not None and item.text is not None
    filename, dataset = item.text.strip().split(":", 1)
    with h5py.File(source.xdmf.parent / filename, "r+") as handle:
        handle[dataset][0] = value


def test_catalog_uses_generator_and_actual_convergence_ids(
    study_data: tuple[Path, Path],
) -> None:
    """Read settings without execution and keep missing/short histories explicit."""
    root, _ = study_data
    conditions, origin = conditions_from_source(root)
    assert [item.boundary_flow for item in conditions] == [-2.5, -2.5]
    assert origin.endswith("feb_manipulation_v2.py")
    frame = simulation_conditions_df(conditions)
    assert frame.index.tolist() == ["sim001", "sim002"]
    assert frame.boundary_flow.dtype == np.dtype("float64")
    rows = {source.run_id: inspect_run(source).row for source in sources(root)}
    assert rows["sim001"]["pattern_key"] == 0
    assert rows["sim002"]["pattern_key"] == 1
    assert rows["sim001"]["num_cells"] == 2
    assert rows["sim001"]["complete"] is True
    assert rows["convergence_06250"]["complete"] is False
    assert rows["convergence_06250"]["resolution_key"] == "06250"


@pytest.mark.parametrize("problem", ["time", "schema", "sidecar"])
def test_invalid_descriptors_fail_explicitly(
    study_data: tuple[Path, Path], problem: str
) -> None:
    """Reject duplicate times, changing fields, and missing HDF data."""
    root, _ = study_data
    source = sources(root)[0]
    tree = ET.parse(source.xdmf)
    grids = tree.findall(".//Grid[@CollectionType='Temporal']/Grid")
    if problem == "time":
        item = grids[1].find("Time")
        assert item is not None
        item.set("Value", "0")
    elif problem == "schema":
        grids[-1].remove(grids[-1].findall("Attribute")[0])
    else:
        source.xdmf.with_suffix(".h5").unlink()
    tree.write(source.xdmf)
    with pytest.raises(ValueError):
        inspect_run(source)


def test_full_hdf_processing_preserves_types_blocks_and_weighting(
    study_data: tuple[Path, Path],
) -> None:
    """Read mixed cells in order and promote changing scalar types without truncation."""
    root, output = study_data
    source = sources(root)[0]
    before = {
        path: path.stat().st_mtime_ns for path in root.rglob("*") if path.is_file()
    }
    result = process_run(source, output, profile="full", legacy="ignore")
    assert result.status == "processed", result.message
    cells, points, summary = load_processed_run(source, output)
    np.testing.assert_array_equal(cells.rr_protein.values[0], [1, 3])
    np.testing.assert_array_equal(cells.counter.values[1], [1.5, 2.5])
    assert cells.rr_protein.dtype == np.dtype("float32")
    assert points.pressure.shape == (5, 4)
    assert "displacement" not in points
    assert summary.rr_protein_mean.values[0] == pytest.approx(2.5)
    assert summary.rr_protein_std.values[0] == pytest.approx(np.sqrt(0.75))
    assert summary.rr_protein_mean.dtype == np.dtype("float64")
    assert before == {path: path.stat().st_mtime_ns for path in before}


def test_legacy_reuse_selects_fields_and_verifies_values(
    study_data: tuple[Path, Path],
) -> None:
    """Adopt matching legacy caches while excluding unused numerical fields."""
    root, output = study_data
    source = sources(root)[0]
    mesh_datasets_from_xdmf(source.xdmf)
    result = process_run(source, output)
    assert result.status == "processed", result.message
    assert "legacy cache (sample validation)" in result.message
    cells, points, _ = load_processed_run(source, output)
    assert set(cells.data_vars) == set(ANALYSIS_FIELDS)
    assert not points.data_vars
    assert cells.rr_protein.values[-1, 1] == 7


def test_full_legacy_validation_detects_unsampled_corruption(
    study_data: tuple[Path, Path],
) -> None:
    """Verify every frame when a stronger audit than first/middle/last is requested."""
    root, output = study_data
    source = sources(root)[0]
    mesh_datasets_from_xdmf(source.xdmf)
    path = source.xdmf.with_name(f"{source.run_id}_cells.nc")
    data = xr.load_dataset(path)
    data.rr_protein.values[1, 0] = 999
    data.to_netcdf(path, engine="scipy")
    result = process_run(source, output, legacy="full")
    assert result.status == "processed", result.message
    assert "legacy rejected" in result.message
    cells, _, _ = load_processed_run(source, output)
    assert cells.rr_protein.values[1, 0] == 2


def test_unchanged_rerun_reads_no_arrays_or_descriptor(
    study_data: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Resume from file fingerprints without reopening XDMF or HDF data."""
    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"

    def unexpected(*args: object, **kwargs: object) -> None:
        raise AssertionError("Unchanged runs must not read source arrays or XML")

    monkeypatch.setattr(cache, "inspect_run", unexpected)
    monkeypatch.setattr(h5py, "File", unexpected)
    result = process_run(source, output)
    assert result.status == "cached", result.message


def test_source_changes_invalidate_cached_data(study_data: tuple[Path, Path]) -> None:
    """Detect changed HDF values even when the descriptor and dimensions stay fixed."""
    root, output = study_data
    source = sources(root)[0]
    mesh_datasets_from_xdmf(source.xdmf)
    assert process_run(source, output).status == "processed"
    modify_hdf(source, "rr_protein", 1, 42)
    assert cache_status(source, output, "analysis", "sample")[0] == "source changed"
    with pytest.raises(ValueError, match="stale"):
        load_processed_run(source, output)
    result = process_run(source, output)
    assert result.status == "processed", result.message
    assert result.message == "HDF source"
    cells, _, _ = load_processed_run(source, output)
    assert cells.rr_protein.values[1, 0] == 42


def test_failed_rebuild_preserves_completed_generation(
    study_data: tuple[Path, Path],
) -> None:
    """Never publish invalid volumes or replace a completed manifest on failure."""
    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"
    manifest = run_directory(source, output) / "manifest.json"
    previous = manifest.read_bytes()
    modify_hdf(source, "element_volume_point_TPM", 0, -1)
    result = process_run(source, output, force=True)
    assert result.status == "failed"
    assert "volumes" in result.message
    assert manifest.read_bytes() == previous
    assert (manifest.parent / json.loads(previous)["generation"] / "cells.nc").is_file()


@pytest.mark.parametrize("output_kind", ["same", "child", "parent", "symlink"])
def test_output_cannot_modify_raw_tree(
    study_data: tuple[Path, Path], output_kind: str
) -> None:
    """Reject overlapping roots and nested output symlinks into raw data."""
    root, output = study_data
    source = sources(root)[0]
    if output_kind == "same":
        output = root
    elif output_kind == "child":
        output = root / "derived"
    elif output_kind == "parent":
        output = root.parent
    else:
        output.mkdir()
        (output / source.study).symlink_to(
            root / source.study, target_is_directory=True
        )
    before = set(root.rglob("*"))
    result = process_run(source, output)
    assert result.status == "failed"
    assert set(root.rglob("*")) == before


def test_cli_dry_run_and_information_are_read_only(
    study_data: tuple[Path, Path],
) -> None:
    """Discover both studies without creating outputs, then export reports locally."""
    root, output = study_data
    args = ["--input-root", str(root), "--output-root", str(output)]
    assert main([*args, "--dry-run"]) == 0
    assert information_main([*args, "--dry-run"]) == 0
    assert not output.exists()
    assert information_main(args) == 0
    table = (output / "catalog.tsv").read_text()
    assert "convergence_06250" in table
    assert "86400" not in table


def test_cli_failure_does_not_discard_other_runs(study_data: tuple[Path, Path]) -> None:
    """Return failure while leaving independent successful runs ready for resumption."""
    root, output = study_data
    source = sources(root)[0]
    modify_hdf(source, "element_volume_point_TPM", 0, -1)
    args = ["--input-root", str(root), "--output-root", str(output), "--workers", "1"]
    assert main(args) == 1
    report = json.loads((output / "processing-report.json").read_text())
    assert sorted(row["status"] for row in report) == [
        "failed",
        "processed",
        "processed",
    ]
    results = process_studies(
        root, output, study="simulations", run_ids=("sim002",), workers=1
    )
    assert results[0].status == "cached"


def test_spawn_workers_process_both_studies(study_data: tuple[Path, Path]) -> None:
    """Run separate HDF handles in spawned workers and reuse completed results."""
    root, output = study_data
    results = process_studies(root, output, workers=2)
    assert [result.status for result in results] == ["processed"] * 3
    assert all(
        result.status == "cached" for result in process_studies(root, output, workers=2)
    )


@pytest.mark.parametrize("extra", [["--runs", "unknown"], ["--workers", "0"]])
def test_cli_rejects_bad_selection(
    study_data: tuple[Path, Path], extra: list[str]
) -> None:
    """Fail clearly rather than silently processing an empty glob."""
    root, output = study_data
    assert main(["--input-root", str(root), "--output-root", str(output), *extra]) == 1
    assert not output.exists()


def test_numerical_imports_do_not_load_rendering() -> None:
    """Keep metadata and cache commands independent of VTK and pyplot imports."""
    code = 'import sys; import research.spt.information; import research.spt.data_processing; assert "pyvista" not in sys.modules; assert "matplotlib.pyplot" not in sys.modules'
    subprocess.run([sys.executable, "-c", code], check=True, env=os.environ.copy())


def test_damaged_output_is_rebuilt(study_data: tuple[Path, Path]) -> None:
    """Detect incomplete local output rather than trusting a manifest alone."""
    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"
    directory = run_directory(source, output)
    manifest = json.loads((directory / "manifest.json").read_text())
    (directory / manifest["generation"] / "cells.nc").unlink()
    assert cache_status(source, output, "analysis", "sample")[0] == "output changed"
    assert process_run(source, output).status == "processed"
    assert load_processed_run(source, output)[0].sizes["cell"] == 2


@pytest.mark.skipif(os.name != "posix", reason="Uses an external POSIX lock holder")
def test_busy_run_fails_then_resumes(study_data: tuple[Path, Path]) -> None:
    """Reject concurrent writers and allow processing after the OS lock is released."""
    import fcntl

    root, output = study_data
    source = sources(root)[0]
    directory = run_directory(source, output)
    directory.mkdir(parents=True)
    with (directory / ".lock").open("a+b") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        assert process_run(source, output).status == "failed"
        assert not (directory / "manifest.json").exists()
    assert process_run(source, output).status == "processed"


def test_missing_simulation_is_reported(study_data: tuple[Path, Path]) -> None:
    """Show expected missing runs and require explicit subsets for partial studies."""
    root, output = study_data
    (root / "simulations_sixth/xdmf/sim002.xdmf").unlink()
    args = ["--input-root", str(root), "--output-root", str(output)]
    assert information_main(args) == 1
    assert "missing" in (output / "catalog.tsv").read_text()
    assert main(args) == 1
    assert (
        main([*args, "--study", "simulations", "--runs", "sim001", "--workers", "1"])
        == 0
    )


def test_summary_does_not_load_spatial_fields(
    study_data: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Read compact summaries without loading larger local cell caches."""
    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"
    original = xr.load_dataset

    def summary_only(path: Path, **kwargs: object) -> xr.Dataset:
        assert path.name == "summary.nc"
        return original(path, engine="scipy")

    monkeypatch.setattr(xr, "load_dataset", summary_only)
    summary = cache.load_summary(source, output)
    assert summary.rr_protein_mean.values[0] == pytest.approx(2.5)


def test_condition_plots_are_headless_and_closed(study_data: tuple[Path, Path]) -> None:
    """Save metadata graphics locally and leave no figures open."""
    from matplotlib import pyplot as plt

    root, output = study_data
    assert (
        information_main(
            [
                "--input-root",
                str(root),
                "--output-root",
                str(output),
                "--plots",
                "--excel",
            ]
        )
        == 0
    )
    assert (output / "figures/boundary_flux.png").is_file()
    assert (output / "figures/boundary_color.png").is_file()
    assert (output / "catalog.xlsx").is_file()
    assert not plt.get_fignums()


def test_failed_publication_keeps_previous_cache(
    study_data: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Recover from a failed pointer replacement without losing completed outputs."""
    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"
    pointer = run_directory(source, output) / "manifest.json"
    previous = pointer.read_bytes()
    original = Path.replace

    def fail_publication(path: Path, target: Path) -> Path:
        if target.name == "manifest.json":
            raise OSError("Simulated publication failure")
        return original(path, target)

    monkeypatch.setattr(Path, "replace", fail_publication)
    assert process_run(source, output, force=True).status == "failed"
    assert pointer.read_bytes() == previous
    assert load_processed_run(source, output)[0].rr_protein.values[0, 0] == 1


def test_information_reuses_current_metadata(
    study_data: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Avoid descriptor reads when an unchanged processed run supplies the catalogue."""
    from research.spt import catalog

    root, output = study_data
    source = sources(root)[0]
    assert process_run(source, output).status == "processed"

    def unexpected(source: RunSource) -> None:
        raise AssertionError("Current catalogue must not reread descriptors")

    monkeypatch.setattr(catalog, "inspect_run", unexpected)
    assert (
        information_main(
            [
                "--input-root",
                str(root),
                "--output-root",
                str(output),
                "--runs",
                source.run_id,
                "--dry-run",
            ]
        )
        == 0
    )


def test_malformed_xml_reports_failure(study_data: tuple[Path, Path]) -> None:
    """Report an interrupted descriptor upload without an uncaught parser traceback."""
    root, output = study_data
    source = sources(root)[0]
    source.xdmf.write_text("<Xdmf>")
    assert (
        information_main(["--input-root", str(root), "--output-root", str(output)]) == 1
    )
    assert process_run(source, output).status == "failed"


def test_relative_output_paths_load_successfully(
    study_data: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Use the documented relative output directory for both processing and loading."""
    root, output = study_data
    source = sources(root)[0]
    monkeypatch.chdir(output.parent)
    relative = Path(output.name)
    assert process_run(source, relative).status == "processed"
    assert load_processed_run(source, relative)[0].sizes["cell"] == 2
    assert cache.load_summary(source, relative).sizes["time"] == 5


def test_changed_generator_order_is_not_silently_relabelled(
    study_data: tuple[Path, Path],
) -> None:
    """Reject a scan with different loop ordering instead of assigning wrong conditions."""
    root, _ = study_data
    script = root / "simulations_sixth/material/feb_manipulation_v2.py"
    script.write_text(
        "zonation_pattern = [0, 1]\nboundary_flow = [-2.5]\nfor flow in boundary_flow:\n    for pattern in zonation_pattern:\n        pass\n"
    )
    with pytest.raises(ValueError, match="ordering"):
        discover_runs(root)
