# Contributing

Open an [issue](https://github.com/matthiaskoenig/porous_media/issues) for a bug
report or proposal. Include a small reproducer, expected and actual behavior,
Python and package versions, and relevant input-format details. Do not include
private research data or credentials.

Development targets `main`. Create a short-lived branch from the current `main`
and submit a pull request against `main`.

## Setup

```bash
git clone https://github.com/matthiaskoenig/porous_media.git
cd porous_media
uv sync --locked --group dev --group docs
uv run --no-sync pre-commit install
```

Python 3.13 or later and [uv](https://docs.astral.sh/uv/) are required.
Use the committed `uv.lock` for reproducible development and CI environments.

## Validate a change

```bash
uv run --no-sync ruff check .
uv run --no-sync ruff format --check .
uv run --no-sync ty check
uv run --no-sync pytest
uv run --no-sync zensical build --clean --strict
uv run --no-sync tox run -e build
```

Update guides and source docstrings when behavior changes. Add meaningful
regression coverage for bug fixes. Prefer bundled or synthetic data so checks
run without external solvers, graphics sessions, or unpublished datasets.

Pull requests require the `tests`, `ruff`, `ty`, `docs`, and `build` checks,
resolved conversations, and linear history. The solo-maintainer policy does
not require approving reviews. Existing release tags cannot be moved or deleted.
See the [development guide](docs/development.md) for the test matrix, isolated
tox environments, documentation preview, and repository policy definitions.

Describe the resulting behavior and validation in the pull request. Follow
[RELEASE.md](RELEASE.md) for releases; do not edit generated changelog files.
