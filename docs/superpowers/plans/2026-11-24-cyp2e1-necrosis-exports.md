# CYP2E1 Necrosis Exports Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish two CYP2E1 tables and their figures from the analysis stage: necrosis over lobular position at 24 h, and necrosis paired with CYP2E1 protein in the same position bin.

**Architecture:** A new focused module `research/spt/analysis_cyp2e1.py` turns the already concatenated profile dataset into two tidy frames and renders their figures. `render_analysis` calls it in the simulations branch, next to the existing convergence-only branch, so both files publish inside the same atomic generation directory.

**Tech Stack:** Python 3.13+, xarray, pandas, matplotlib (headless `Figure` + `FigureCanvasAgg` via `research.spt.plot_render`), pytest.

**Spec:** `docs/superpowers/specs/2026-11-24-cyp2e1-necrosis-exports-design.md`

## Global Constraints

- CYP2E1 is selected by `pattern_name(condition.pattern_key) == "cyp2e1"`, using `pattern_name` from `research.spt.conditions`. Never hard-code run ids, and do not reuse `MIXED_EXCLUDED_PATTERN`, which names the same pattern for an unrelated reason.
- Both tables take necrosis and protein from the same, final spatial time: `float(profiles.time.values[-1])`.
- Position bins: 20 in production, centers 0.025 to 0.975; tests use fewer via `AnalysisOptions(bins=...)`.
- Field names come from `research.spt.fields`: necrosis is `rr_necrosis`, protein is `rr_protein`.
- Dose columns follow the existing convention in `analysis_plots.py`: `dose` is one-based (`source.dose + 1`).
- Curve color is `simulation_color(pattern_key, dose_key)` from `research.spt.plot_style`, with the zero-based dose key.
- Figures are written with `figure(...)` and `save_figure(...)` from `research.spt.plot_render`, under the `RenderOptions(formats=options.formats, dpi=options.dpi)` already built in `render_analysis`.
- Files are written into the staging directory `render_analysis` already holds, never directly into the published generation.
- No em dash anywhere in code, comments, docs or commit messages.
- No attribution lines in commit messages.

## Review Focus

- A selection with no CYP2E1 run (for example `--runs sim001`): both exports are skipped and no empty TSV or figure appears. Covered in Task 1 Step 1 and Task 3 Step 5.
- A run whose catalogue row carries no condition (`source.run.condition is None`): it is ignored rather than raising. Covered in Task 1 Step 6.
- A position bin holding no cells (`cell_count == 0`, `mean` is NaN): the row stays in the table with NaN, so bin indexing is not silently shortened. Covered in Task 1 Step 10.
- Exactly one CYP2E1 run in the selection: both figures still render, with a one-entry legend. Covered in Task 2 Step 6.
- Several output formats at once (`formats=("svg", "png")`): each figure appears in both. Covered in Task 3 Step 5.

---

### Task 1: Build the two frames

**Files:**
- Create: `research/spt/analysis_cyp2e1.py`
- Test: `research/tests/test_spt_cyp2e1.py`

**Interfaces:**
- Consumes: `PlotSource` from `research.spt.plot_data` (`run.run_id`, `run.condition`, `dose`), `pattern_name` from `research.spt.conditions`.
- Produces: `cyp2e1_frames(profiles: xr.Dataset, sources: list[PlotSource]) -> tuple[pd.DataFrame, pd.DataFrame]`, returning `(necrosis_over_position, necrosis_with_protein)`. Both are empty `pd.DataFrame()` when the selection holds no CYP2E1 run.

- [ ] **Step 1: Write the failing test for an empty selection**

```python
def test_no_cyp2e1_runs_produce_no_frames() -> None:
    profiles, sources = _profiles(patterns=[0], doses=[0], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert necrosis.empty and paired.empty
```

Write `_profiles(patterns, doses, bins)` in the test module as a helper that builds an `xr.Dataset` with dims `run`, `time`, `field`, `position`, data variables `mean`, `std`, `minimum`, `maximum`, `cell_count`, `volume`, `bin_left`, `bin_right`, times `(0.0, 86400.0)`, fields `("rr_necrosis", "rr_protein")`, and a matching list of `PlotSource` objects whose `run.condition` carries the pattern and dose. Build the `PlotSource` objects with `dataclasses.replace` on a source from `research.tests.test_spt_geometry._source`-backed catalogue rows, or construct `SimulationCondition` directly and attach it to a `RunSource`.

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: FAIL with `ImportError` or `NameError` for `cyp2e1_frames`.

- [ ] **Step 3: Implement `cyp2e1_frames(profiles: xr.Dataset, sources: list[PlotSource]) -> tuple[pd.DataFrame, pd.DataFrame]` in `research/spt/analysis_cyp2e1.py`**

Select the CYP2E1 sources, return two empty frames when none remain, otherwise select `profiles.sel(run=<those ids>, time=<final time>)` and build the two frames from `rr_necrosis` and `rr_protein`.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS

- [ ] **Step 5: Write the failing test for the necrosis frame's columns and shape**

```python
def test_necrosis_frame_columns_and_rows() -> None:
    profiles, sources = _profiles(patterns=[6], doses=[0, 1], bins=4)
    necrosis, _ = cyp2e1_frames(profiles, sources)
    assert list(necrosis.columns) == [
        "run",
        "pattern",
        "pattern_name",
        "dose",
        "boundary_flow",
        "time",
        "position",
        "bin_left",
        "bin_right",
        "necrosis_mean",
        "necrosis_std",
        "necrosis_minimum",
        "necrosis_maximum",
        "cell_count",
        "volume",
    ]
    assert len(necrosis) == 8
    assert set(necrosis.pattern_name) == {"cyp2e1"}
    assert list(necrosis.dose.unique()) == [1, 2]
    assert necrosis.time.unique().tolist() == [86400.0]
```

- [ ] **Step 6: Write the failing test for a source without a condition**

```python
def test_source_without_condition_is_ignored() -> None:
    profiles, sources = _profiles(patterns=[6, 6], doses=[0, 1], bins=4)
    sources[1] = replace(sources[1], run=replace(sources[1].run, condition=None))
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert list(necrosis.run.unique()) == [sources[0].run.run_id]
    assert list(paired.run.unique()) == [sources[0].run.run_id]
```

- [ ] **Step 7: Run both tests to verify they fail**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: FAIL on the column list and on the ignored source.

- [ ] **Step 8: Extend `cyp2e1_frames` so both tests pass**

- [ ] **Step 9: Run the tests to verify they pass**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS

- [ ] **Step 10: Write the failing tests for the paired frame and for empty bins**

```python
def test_paired_frame_pairs_protein_and_necrosis_in_one_bin() -> None:
    profiles, sources = _profiles(patterns=[6], doses=[0], bins=4)
    _, paired = cyp2e1_frames(profiles, sources)
    assert list(paired.columns) == [
        "run",
        "pattern",
        "pattern_name",
        "dose",
        "boundary_flow",
        "time",
        "position",
        "bin_left",
        "bin_right",
        "protein_mean",
        "protein_std",
        "necrosis_mean",
        "necrosis_std",
    ]
    assert len(paired) == 4


def test_empty_bins_keep_their_rows() -> None:
    profiles, sources = _profiles(patterns=[6], doses=[0], bins=4, empty_bin=1)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    assert len(necrosis) == 4 and len(paired) == 4
    assert int(necrosis.cell_count.iloc[1]) == 0
    assert np.isnan(necrosis.necrosis_mean.iloc[1])
    assert np.isnan(paired.protein_mean.iloc[1])
```

Give `_profiles` an `empty_bin: int | None = None` argument that sets `cell_count` to 0 and every statistic to NaN in that bin.

- [ ] **Step 11: Run the tests to verify they fail**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: FAIL on the paired frame columns.

- [ ] **Step 12: Extend `cyp2e1_frames` so every test passes**

- [ ] **Step 13: Run the full test module**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS, 5 tests

- [ ] **Step 14: Commit**

```bash
git add research/spt/analysis_cyp2e1.py research/tests/test_spt_cyp2e1.py
git commit -m "feat: build CYP2E1 necrosis and protein export frames"
```

---

### Task 2: Render the two figures

**Files:**
- Modify: `research/spt/analysis_cyp2e1.py`
- Test: `research/tests/test_spt_cyp2e1.py`

**Interfaces:**
- Consumes: `cyp2e1_frames` from Task 1, `figure` and `save_figure` and `RenderOptions` from `research.spt.plot_render`, `simulation_color` from `research.spt.plot_style`.
- Produces: `render_cyp2e1_figures(necrosis: pd.DataFrame, paired: pd.DataFrame, directory: Path, options: RenderOptions) -> list[Path]`, returning the files it wrote, empty when either frame is empty.

- [ ] **Step 1: Write the failing test for the rendered files**

```python
def test_figures_are_written_in_every_format(tmp_path: Path) -> None:
    profiles, sources = _profiles(patterns=[6], doses=[0, 1], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    written = render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("svg", "png"), dpi=45)
    )
    assert sorted(path.name for path in written) == [
        "cyp2e1-necrosis-position.png",
        "cyp2e1-necrosis-position.svg",
        "cyp2e1-necrosis-protein.png",
        "cyp2e1-necrosis-protein.svg",
    ]
    assert all(path.is_file() for path in written)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: FAIL with `NameError` for `render_cyp2e1_figures`.

- [ ] **Step 3: Implement `render_cyp2e1_figures(necrosis: pd.DataFrame, paired: pd.DataFrame, directory: Path, options: RenderOptions) -> list[Path]`**

One axis per figure. `cyp2e1-necrosis-position` plots `necrosis_mean` against `position`, one line per `dose`. `cyp2e1-necrosis-protein` plots `necrosis_mean` against `protein_mean`, one line per `dose`, rows ordered by `position` so each series runs periportal to perivenous. Color each line with `simulation_color(6, dose - 1)`. Label the axes with the lobular position, the necrosis fraction and the protein level, and add a `D1` to `D8` legend built from the `dose` column.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS

- [ ] **Step 5: Write the failing test for axis labels and the legend**

```python
def test_figure_labels_and_legend(tmp_path: Path) -> None:
    profiles, sources = _profiles(patterns=[6], doses=[0, 7], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("svg",), dpi=45)
    )
    document = (tmp_path / "cyp2e1-necrosis-position.svg").read_text()
    assert "D1" in document and "D8" in document
    assert "Necrosis" in document and "position" in document.lower()
    assert "CYP2E1" in (tmp_path / "cyp2e1-necrosis-protein.svg").read_text()
```

- [ ] **Step 6: Write the failing test for a single run and for empty frames**

```python
def test_single_run_and_empty_frames(tmp_path: Path) -> None:
    profiles, sources = _profiles(patterns=[6], doses=[3], bins=4)
    necrosis, paired = cyp2e1_frames(profiles, sources)
    written = render_cyp2e1_figures(
        necrosis, paired, tmp_path, RenderOptions(formats=("png",), dpi=45)
    )
    assert len(written) == 2
    assert (
        render_cyp2e1_figures(
            pd.DataFrame(),
            pd.DataFrame(),
            tmp_path,
            RenderOptions(formats=("png",), dpi=45),
        )
        == []
    )
```

- [ ] **Step 7: Run the tests to verify they fail**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: FAIL on the labels and on the empty-frame guard.

- [ ] **Step 8: Extend `render_cyp2e1_figures` so every test passes**

- [ ] **Step 9: Run the full test module**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS, 8 tests

- [ ] **Step 10: Commit**

```bash
git add research/spt/analysis_cyp2e1.py research/tests/test_spt_cyp2e1.py
git commit -m "feat: render CYP2E1 necrosis position and protein figures"
```

---

### Task 3: Publish both files from the analysis stage

**Files:**
- Modify: `research/spt/analysis_plots.py` (the staging block inside `render_analysis`, around lines 683 to 737)
- Modify: `research/spt/analysis.py:197-199` (the `row` dict that records published table paths)
- Test: `research/tests/test_spt_cyp2e1.py`

**Interfaces:**
- Consumes: `cyp2e1_frames` and `render_cyp2e1_figures` from Tasks 1 and 2.
- Produces: `cyp2e1-necrosis-position.tsv` and `cyp2e1-necrosis-protein.tsv` in the published generation directory, plus the two figures; `analysis-<study>.tsv` gains the columns `cyp2e1_position` and `cyp2e1_protein` holding those paths when they exist.

- [ ] **Step 1: Write the failing end-to-end test**

```python
def test_analysis_publishes_cyp2e1_exports(tmp_path: Path) -> None:
    root, output = _cyp2e1_study(tmp_path, ("sim049", "sim050"))
    assert [
        item.status
        for item in process_analysis(
            root,
            output,
            output,
            study="simulations",
            run_ids=("sim049", "sim050"),
            options=AnalysisOptions(time_points=3, bins=4),
            plots=PlotOptions(formats=("svg", "png"), dpi=45),
            workers=1,
        )
    ] == ["processed", "processed"]
    directory = output / "analysis/simulations_lobulus"
    generation = (
        directory / json.loads((directory / "manifest.json").read_text())["generation"]
    )
    table = pd.read_csv(generation / "cyp2e1-necrosis-position.tsv", sep="\t")
    assert len(table) == 8
    assert sorted(table.dose.unique()) == [1, 2]
    assert set(table.pattern_name) == {"cyp2e1"}
    paired = pd.read_csv(generation / "cyp2e1-necrosis-protein.tsv", sep="\t")
    assert {"protein_mean", "necrosis_mean"} <= set(paired.columns)
    for name in ("cyp2e1-necrosis-position", "cyp2e1-necrosis-protein"):
        assert (generation / f"{name}.svg").is_file()
        assert (generation / f"{name}.png").is_file()
```

Write `_cyp2e1_study(tmp_path, run_ids)` in the test module, modelled on the `study` fixture in `research/tests/test_spt_analysis.py:73`: build each run with `_input` from that module, then `process_geometries(..., study="simulations", run_ids=run_ids, combinations=False, workers=1)` and `process_studies(...)` with the same run ids.

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py::test_analysis_publishes_cyp2e1_exports -q`
Expected: FAIL with `FileNotFoundError` for `cyp2e1-necrosis-position.tsv`.

- [ ] **Step 3: Call the new module from `render_analysis`**

In the staging block, after `spatial-profiles.tsv` is written and beside the existing convergence branch, build the frames from `all_profiles` and `sources`, write each non-empty frame with `to_csv(..., sep="\t", index=False)`, and render the figures inside the existing `mpl.rc_context` block.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py::test_analysis_publishes_cyp2e1_exports -q`
Expected: PASS

- [ ] **Step 5: Write the failing test that a selection without CYP2E1 writes nothing**

```python
def test_selection_without_cyp2e1_writes_no_exports(study: tuple[Path, Path]) -> None:
    root, output = study
    process_analysis(
        root,
        output,
        output,
        study="simulations",
        run_ids=("sim001",),
        options=AnalysisOptions(time_points=3, bins=4),
        plots=PlotOptions(formats=("png",), dpi=45),
        workers=1,
    )
    directory = output / "analysis/simulations_lobulus"
    generation = (
        directory / json.loads((directory / "manifest.json").read_text())["generation"]
    )
    assert not list(generation.glob("cyp2e1-*"))
```

Import the `study` fixture by re-exporting it from `research/tests/test_spt_analysis.py`, or move it into `research/tests/conftest.py` if that keeps both modules simpler.

- [ ] **Step 6: Run test to verify it fails or passes, then make it pass**

Run: `uv run --no-sync pytest research/tests/test_spt_cyp2e1.py -q`
Expected: PASS once the guard from Task 1 Step 3 is honoured end to end.

- [ ] **Step 7: Record the published paths in `analysis.py`**

Add `cyp2e1_position` and `cyp2e1_protein` to the `row` dict at `research/spt/analysis.py:197-199`, each set to the generation path when the file exists and left out otherwise.

- [ ] **Step 8: Run the whole research suite**

Run: `uv run --no-sync pytest research/tests -q`
Expected: PASS, no regressions

- [ ] **Step 9: Commit**

```bash
git add research/spt/analysis_plots.py research/spt/analysis.py research/tests/test_spt_cyp2e1.py
git commit -m "feat: publish CYP2E1 necrosis exports from the analysis stage"
```

---

### Task 4: Document the exports

**Files:**
- Modify: `docs/research/spt.md` (the results layout block and the analysis section)

- [ ] **Step 1: Add both files to the results layout block**

List `cyp2e1-necrosis-position.tsv`, `cyp2e1-necrosis-protein.tsv` and their figures under `analysis/<study>/`.

- [ ] **Step 2: Describe the exports**

One short subsection: what each file holds, its columns, that both quantities are read at the final spatial time and why (the mesh deforms, so cells cross position bins and the protein profile at 0 s is not the one at 24 h), and that the files appear only when the selection holds CYP2E1 runs.

- [ ] **Step 3: Name the known bad run**

State that `sim050` in the 2026-11-24 dataset carries an imposed fluid flux about ten times the whole D1-D8 family, so its rows and its curve are not trustworthy until the run is repeated.

- [ ] **Step 4: Verify lint, types and tests**

Run: `uv run --no-sync ruff check . && uv run --no-sync ruff format --check . && uv run --no-sync ty check && uv run --no-sync pytest -q`
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add docs/research/spt.md
git commit -m "docs: describe the CYP2E1 necrosis exports"
```
