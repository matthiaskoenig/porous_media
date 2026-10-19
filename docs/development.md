# Development

See [contributing](contributing.md) for bug reports, contribution guidelines,
and the pull request process.

## Environment

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then clone
the repository and sync the committed lockfile:

```bash
git clone https://github.com/matthiaskoenig/porous_media.git
cd porous_media
uv sync --locked --group dev --group docs
uv run --no-sync pre-commit install
```

Python 3.13 is the minimum supported version. `.python-version` selects the
default development interpreter. To check the minimum version explicitly:

```bash
uv sync --locked --group dev --group docs --python 3.13
```

`uv.lock` records the resolved environment; include intentional lockfile changes
in pull requests that update dependencies.

The `dev` dependency group contains Ruff, ty, pytest, tox, pre-commit, build
validation, and release tooling. The `docs` group contains Zensical and
mkdocstrings. `uv run --no-sync` uses the environment prepared by the sync step.

## Update dependencies

Resolve current stable releases and sync the development environment:

```bash
uv lock --upgrade
uv sync --locked --group dev --group docs
```

Update minimum versions in `pyproject.toml` when raising dependency requirements,
and commit the regenerated `uv.lock` with the change. Keep the CI uv version and
tox requirements consistent with the development environment. Verify the full
test matrix before merging a dependency update. Transitive dependencies can
remain below their newest release when an upstream package requires an older
version; resolve compatibility through the upstream package rather than
overriding its requirements.

## Required checks

Run the same checks used by CI:

```bash
uv run --no-sync ruff check .
uv run --no-sync ruff format --check .
uv run --no-sync ty check
uv run --no-sync pytest
uv run --no-sync zensical build --clean --strict
uv run --no-sync tox run -e build
```

Ruff handles linting and formatting. ty checks types. pytest exercises the
bundled workflows and regressions. The strict documentation build includes the
source-based API reference. Package checks build both distribution formats and
validate distribution metadata.

To format files, run `uv run --no-sync ruff format .`. To check all pre-commit
hooks, run `uv run --no-sync pre-commit run --all-files`.

Tox provides isolated environments:

```bash
uv run --no-sync tox run -e ty,ruff,docs,py313,py314,build
```

Install the corresponding Python interpreters or let uv provision them.
CI tests Python 3.13 and 3.14 on Linux and Python 3.14 on Windows and macOS.
The `tests` job aggregates the supported matrix into a stable required check.

An experimental Linux job also attempts the installed-wheel tests on Python
3.15. This job can fail without blocking the required checks while upstream
dependencies gain compatible distributions. The current blocker is VTK, which
has neither a Python 3.15 wheel nor a source distribution on PyPI. The package
does not yet advertise Python 3.15 support in its metadata or change the default
development interpreter. Reproduce the experimental check separately:

```bash
uv python install 3.15
uv run --no-sync tox run -e py315
```

Keep `py315` out of the default tox environments until installation and the
complete test suite pass. At that point, make the CI job required, add the
Python 3.15 classifier, and update the supported-version documentation.

## Documentation

Preview the site locally:

```bash
uv run --no-sync zensical serve
```

Edit pages under `docs/` and navigation in `zensical.toml`. The build renders API
pages from source docstrings using mkdocstrings. Keep examples executable with
bundled resources or label external data requirements explicitly.
The numerical [quickstart](quickstart.md) and nonvisual
[zonation example](zonation.md#bundled-example) do not require a display.

The documentation workflow publishes successful builds from `main` to GitHub
Pages. The configured URL is
[`https://matthiaskoenig.github.io/porous_media/`](https://matthiaskoenig.github.io/porous_media/).
Publication requires GitHub Pages to use the Actions deployment source and a
successful deployment run.

## Repository policies

The repository ruleset for `main` requires:

- a pull request and resolved review conversations;
- passing `tests`, `ruff`, `ty`, `docs`, and `build` checks;
- linear history, with no force pushes or branch deletion.

The solo-maintainer policy requires zero approving reviews so the maintainer
can merge their own pull requests after checks pass. Contributor pull requests
still need maintainer review. Tag rules prevent moving or deleting existing
tags, preserving release provenance.
GitHub's additional approval requirement for unattributed changes remains enabled.

The versioned definitions live in
[`.github/rulesets/`](https://github.com/matthiaskoenig/porous_media/tree/main/.github/rulesets).
Repository settings must also be applied on GitHub; committing a JSON definition
alone does not activate it.

Maintainers with repository admin access can apply the definitions with the
authenticated `gh-axi` wrapper:

```bash
bash .github/rulesets/apply.sh
```

The script updates existing rulesets by name and applies the merge settings, so
running it again does not create duplicate rulesets. Review policy changes in a
pull request before applying them.

## Project-specific research code

Reusable modules live under `porous_media.data`, `porous_media.mesh`,
`porous_media.visualization`, and `porous_media.febio`. The `porous_media.analyses`
tree contains research workflows and historical scripts. Some require
unpublished datasets, local paths, or external solvers. Avoid adding those
requirements to portable tests or presenting those scripts as general APIs.

## Releases

Follow [release guide](release.md)
for version changes, release notes, and immutable release tags. Do not edit
generated files or historical changelog files manually.
