# Agent instructions

## Scope and project layout

- Follow the contribution guide in `docs/contributing.md` and release procedure
  in `docs/release.md`. Keep instructions consistent with the actual configuration.
- Python 3.13 is the minimum supported version. Development defaults to the
  interpreter in `.python-version`; tests also cover Python 3.14. Python 3.15
  coverage is experimental until upstream dependencies provide distributions.
- Reusable package code lives in `src/porous_media/data`, `mesh`,
  `visualization`, and `febio`. Tests live under `tests/`.
- `src/porous_media/analyses` contains research workflows and historical
  scripts. Do not treat these as general public APIs or make their unpublished
  datasets, local paths, or external solvers prerequisites for the package.
- `src/porous_media/resources` contains bundled examples distributed with the
  package. Preserve their inclusion in wheels and source distributions.

## Environment and dependencies

Use uv and the committed lockfile:

```bash
uv sync --locked --group dev --group docs
uv run --no-sync pre-commit install
```

- Use `uv run --no-sync` after syncing the environment. Use Ruff for linting and
  formatting and ty for type checking. Do not add mypy dependencies or
  configuration, a second package manager, or another lockfile.
- Declare runtime dependencies and development/documentation groups in
  `pyproject.toml`. Regenerate `uv.lock` with uv when dependencies change.
  Check that dependency updates support both Python 3.13 and 3.14.
- Keep tool requirements, hook commands, tox environments, and CI consistent.
  Do not edit `uv.lock` or other generated files manually.

## Code and data safety

- Follow `.ruff.toml`: 88-character lines, double quotes, sorted imports, and
  Google-style docstrings. Add accurate type annotations for changed APIs.
- Fix lint and type errors at their source. Do not silence checks with broad
  exclusions, blanket ignores, or unneeded `Any` annotations.
- Use `pathlib.Path` for file operations. Library functions must not change the
  process working directory.
- Use `porous_media.data.xdmf_tools.XDMFTimeSeriesWriter` for HDF-backed XDMF
  output. Keep HDF sidecars beside their XDMF files with relative references so
  the pair remains relocatable. Never create or overwrite an unrelated sidecar
  in the working directory.
- Close files, readers, writers, images, and plots reliably. Preserve headless
  execution for noninteractive workflows.
- Do not modify research datasets, historical results, bundled binary fixtures,
  or generated simulation outputs unless the requested change requires it.
  Do not commit caches, environments, build outputs, or private data.
- Never manually edit `CHANGELOG.md` or files marked as generated. Write normal
  prose in documentation, comments, commits, and pull requests. Do not add
  attribution lines or use the em dash character.

## Validation

Run the checks appropriate to the change; the complete local checks are:

```bash
uv run --no-sync ruff check .
uv run --no-sync ruff format --check .
uv run --no-sync ty check
uv run --no-sync pytest -q
uv run --no-sync zensical build --clean --strict
uv run --no-sync tox run -e build
uv run --no-sync pre-commit run --all-files
```

- ty treats warnings as errors. Keep its checks over source, tests, and scripts.
- Tox provides isolated, locked environments, including installed-wheel tests:
  `uv run --no-sync tox run -e ruff,ty,py313,py314,docs,build`.
- CI tests installed wheels on Linux with Python 3.13 and 3.14, and on Windows
  and macOS with Python 3.14. Preserve package metadata and resource coverage.
  Python 3.15 has a nonblocking Linux job and an opt-in `py315` tox environment;
  do not claim full support until its installation and tests pass.
- Reproduce a bug through the user-visible workflow before changing code. Add
  regression coverage that checks behavior rather than implementation details.
- Prefer synthetic fixtures or bundled resources and pytest's `tmp_path` for
  outputs. Tests must not require private research data or FEBio installations.
- Use `MPLBACKEND=Agg` and `PYVISTA_OFF_SCREEN=true` for headless rendering.
  FFmpeg integration tests require the FFmpeg executable; preserve their
  explicit availability checks.

## Documentation and repository workflow

- Edit documentation under `docs/` and navigation in `zensical.toml`. Zensical
  renders the API reference from source docstrings through mkdocstrings.
- Keep project documentation under `docs/`, including `docs/release.md` and
  `docs/release-notes/`. The root `README.md` is the project overview;
  `AGENTS.md` and GitHub templates provide repository instructions. Root
  `LICENSE` and `CITATION.cff` provide license and citation metadata.
- Keep installation instructions, quickstarts, and contributor/release guidance
  in `docs/`. The README should contain the overview, documentation links,
  citation, license, and funding information without duplicating guides.
- Keep contribution information canonical in `docs/contributing.md`. Root
  contribution links should point readers to that guide instead of duplicating
  its instructions.
- Update affected guides, examples, and docstrings when behavior changes. Use
  bundled or synthetic data in examples, or state external requirements.
- Start work from the current `main` on a short-lived branch. Preserve existing
  user edits and deletions; do not reset, overwrite, or stage unrelated work.
- Merge changes through a pull request against `main`. Required checks are
  `tests`, `ruff`, `ty`, `docs`, and `build`; resolve review conversations and
  preserve linear history. Do not push directly to `main` or bypass protections.
- Keep `.github/workflows/ci-cd.yml` and its README badge aligned. Release tags
  validate notes under `docs/release-notes/` before publishing.
- Versioned branch and tag rules live in `.github/rulesets/`. Keep their
  definitions and documented policy consistent. Committing definitions alone
  does not apply GitHub settings. Existing release tags must not move or vanish.
- Describe resulting behavior and relevant validation in pull requests. Follow
  `docs/release.md` for version changes, release notes, and tags.
