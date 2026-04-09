# Release

* update release notes in `release-notes` with commit
* make sure all tests run (`tox -p`)
* check formating and linting (`ruff check`)
* test bump version (`uvx bump-my-version bump [major|minor|patch] --dry-run -vv`)
* bump version (`uvx bump-my-version bump [major|minor|patch]`)
* `git push --tags` (triggers release)
* `git push`
* test installation in virtualenv from pypi
```bash
uv venv --python 3.14
uv pip install porous_media
```

# Development
## Setup environment
```bash
uv sync --group dev
```
Setup pre-commit
```bash
pre-commit install
pre-commit run
```

## Testing
Run single target
```bash
tox r -e py314
```
Run all tests in parallel
```bash
tox run-parallel
```


