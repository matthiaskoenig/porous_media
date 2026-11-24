# CYP2E1 necrosis exports - design

## Intent and why

The CYP2E1 cohort (pattern 6, `sim049`-`sim056`, one run per substrate flux D1-D8)
needs two tabular deliverables that can be read outside this repository, plus one
figure each. Today the only position-resolved table is `spatial-profiles.tsv`,
which covers every run and every field but carries no pattern or flux column, so
isolating the CYP2E1 necrosis response means knowing by heart that `sim049`
through `sim056` are pattern 6.

## Goals

- `cyp2e1-necrosis-position.tsv`: necrosis over lobular position at 24 h, one block
  per substrate flux, with the flux value in the table.
- `cyp2e1-necrosis-protein.tsv`: necrosis at 24 h paired with the CYP2E1 protein
  level in the same position bin, one block per substrate flux.
- One figure per file.
- Both files travel inside the published analysis generation, with its manifest and
  provenance, and are reproduced by a normal pipeline run.

## Non-goals

- No other zonation pattern. The files cover pattern 6 only.
- No time resolution beyond 24 h. Position-resolved profiles exist at 0 s, 36000 s
  and 86400 s; these exports use 86400 s alone.
- No change to the existing `spatial-profiles.tsv` columns. Adding `pattern` and
  `dose` there was offered and not taken up; it stays open.
- No new simulation and no re-reading of XDMF/HDF sources.

## Anti-goals

- A second command the user must remember to run, or output that goes stale when
  the analysis is regenerated.
- Silently dropping a run, or inventing a label for a pattern this repository does
  not name.

## Constraints

- A lobule's position bins are 20 wide, centers 0.025 to 0.975, so each file holds
  8 runs x 20 bins = 160 data rows.
- Protein and necrosis must come from the same 86400 s profile. Per-cell protein
  does not change, but the mesh deforms, so cells cross position bins and the
  binned protein profile at 0 s is not the one at 86400 s. Pairing across times
  would mismatch the bins.
- The analysis stage publishes a generation directory atomically under a lock;
  anything written must go into the same staging directory before publication.

## Approaches

- **A. Family-conditional export inside `render_analysis` (chosen).** Mirrors the
  existing `convergence-means.tsv` / `convergence-profiles.tsv` branch in
  `research/spt/analysis_plots.py:613`, which already writes extra tables and
  figures for one family only. Inherits the manifest, atomic publish, provenance
  and report listing. Chosen by the user.
- **B. Standalone command module** reading the published `analysis.nc`. No pipeline
  change and easy to re-point at another results tree, but the output lands outside
  the generation directory with no manifest, and goes stale when the analysis is
  regenerated.
- **C. A general `--export-profile <pattern>:<field>` option.** More machinery than
  two deliverables earn, and each figure still needs its own special case.

## The design

Both tables are built in `render_analysis`, in the same staging directory as the
existing tables, when the family is `simulations`. The source is `all_profiles`,
already concatenated over runs, selected to `spatial_time = 86400`, fields
`rr_necrosis` and `rr_protein`, and to the runs whose `condition.pattern_key` names
the CYP2E1 pattern. Flux columns come from each source's `SimulationCondition`, which
already carries `boundary_flow_key` and `boundary_flow`.

`cyp2e1-necrosis-position.tsv`

| column | meaning |
| --- | --- |
| `run` | run id, `sim049` to `sim056` |
| `pattern` | pattern index, 6 |
| `pattern_name` | `cyp2e1` |
| `dose` | one-based flux index, 1 to 8, matching the D1-D8 labels |
| `boundary_flow` | the substrate flux of that run |
| `time` | 86400.0 |
| `position` | bin center, 0.025 to 0.975 |
| `bin_left`, `bin_right` | bin edges |
| `necrosis_mean`, `necrosis_std`, `necrosis_minimum`, `necrosis_maximum` | necrosis in that bin |
| `cell_count`, `volume` | bin occupancy, carried through from the profile |

`cyp2e1-necrosis-protein.tsv` has the same identifying columns and, in place of the
bin occupancy block, `protein_mean`, `protein_std`, `necrosis_mean`, `necrosis_std`,
so each row pairs the two quantities in one bin.

Figures, rendered in the same `mpl.rc_context` block as the existing analysis
figures and in the formats the stage already uses:

- `cyp2e1-necrosis-position`: necrosis on y against lobular position on x, one curve
  per substrate flux, colored by `simulation_color(6, dose_key)`, which is the
  `YlOrBr` shade already used for this pattern everywhere else. Periportal to
  perivenous direction and axis labels follow the existing spatial figures.
- `cyp2e1-necrosis-protein`: necrosis on y against protein on x, one series per
  substrate flux in the same colors, points ordered by position so each series reads
  as a path from periportal to perivenous.

Both figures carry a D1-D8 legend consistent with the study figures.

## Calls made here, open to correction

- All eight runs are included, `sim050` among them, although its imposed fluid flux
  is about ten times the whole D1-D8 family and its results are known bad. Dropping
  data silently is worse than showing it. The figure captions and
  `docs/research/spt.md` will name `sim050` as carrying a wrong boundary flow until
  it is re-run. No column marks it, because that would invent a concept the rest of
  the tree does not have.
- The CYP2E1 pattern is selected by `pattern_name(...) == "cyp2e1"`, the name already
  defined in `research/spt/conditions.py`, not by hard-coded run ids.
- If the cohort is ever absent from a selection, the export is skipped quietly rather
  than failing, in the same way the convergence tables simply do not appear for the
  simulations family.

## Left to the builder

- Exact figure sizing, legend placement and whether the two figures share a helper
  with the existing `_curves` renderer.
- Whether the two tables share one filtered frame internally.
- Test placement under `research/tests/`, following the existing analysis tests.
