# Release

Release from `main` after the required checks pass. The release tag must identify
the reviewed version commit on `main`. Existing tags are immutable.

## Prepare the release pull request

1. Start from an up-to-date `main` and create a release branch.
2. Sync the locked environment and run validation:

   ```bash
   uv sync --locked --group dev --group docs
   uv run --no-sync ruff check .
   uv run --no-sync ruff format --check .
   uv run --no-sync ty check
   uv run --no-sync pytest
   uv run --no-sync zensical build --clean --strict
   uv run --no-sync tox run -e build
   ```

3. Preview the version change, selecting `major`, `minor`, or `patch`:

   ```bash
   uv run --no-sync bump-my-version bump patch --dry-run -vv
   ```

4. Apply the version change without creating a tag or commit:

   ```bash
   uv run --no-sync bump-my-version bump patch --no-tag --no-commit
   ```

5. Add `docs/release-notes/<version>.md` describing user-visible changes and migration
   requirements. Update `uv.lock` if the project metadata change requires it and
   verify that the locked environment remains current. Do not edit generated
   changelog files manually.
   Validate the version and nonempty release notes locally:

   ```bash
   uv run --no-sync python scripts/validate_release.py <version>
   ```

6. Commit the version files and release notes, open a pull request against
   `main`, and wait for all required checks. Merge through the protected branch.

## Preview the release workflow

Run the `CI-CD` workflow manually with `release_tag` set to the prospective
version. Select the release branch to validate it before merging, or `main`
after merging. The run validates the version and release notes, runs the
checks, and uploads distributions. Manual runs never publish to PyPI, create a
GitHub release, or deploy documentation, including runs selected on a tag.

For example, after pushing the release branch:

```bash
gh workflow run ci-cd.yml --ref <release-branch> -f release_tag=<version>
```

An empty `release_tag` runs the ordinary validation checks for supported Python
versions. Python 3.15 CI coverage is temporarily disabled.

## Tag the reviewed commit

Fetch the merged `main`, verify its version and release notes, and create the
annotated tag using the version number without a `v` prefix, for example:

```bash
git switch main
git pull --ff-only
git tag -a 0.3.1 -m "Release 0.3.1"
git push origin 0.3.1
```

Use the actual release version. Push only that tag. Do not tag a release branch
commit before merging: squash or rebase merging can change the final commit.
Never move or delete a published tag. Fix a release with a new version.

## Verify publication

The tag workflow runs validation, builds distributions, and publishes to PyPI
through the `pypi` environment using trusted publishing. The repository owner
must configure the PyPI trusted publisher and any environment approval policy
before the first automated publication.
Use `.github/workflows/ci-cd.yml` as the trusted publisher workflow filename.
Update an existing publisher registration if the workflow filename changes.
The publisher must match repository owner `matthiaskoenig`, repository
`porous_media`, workflow `ci-cd.yml`, and environment `pypi`. In particular,
replace an older registration for `main.yml`; GitHub settings alone cannot
update the PyPI registration.

If PyPI publication succeeds but the GitHub release step fails, rerun the failed
release job using its original build artifact. Publication skips files already
present on PyPI so the job can finish creating the GitHub release. Do not move
the tag, rebuild different files under the same version, or delete published
files to retry. A manual workflow run is a preview and does not publish.

Check the GitHub Actions release run, the GitHub release notes, and PyPI.
Test the published package in a fresh environment outside the checkout:

```bash
uv venv /tmp/porous-media-release-check --python 3.14
uv pip install --python /tmp/porous-media-release-check/bin/python porous_media==0.3.1
/tmp/porous-media-release-check/bin/python -c "import porous_media; print(porous_media.__version__)"
```

Replace the example version and temporary path as needed. Run the documented
quickstart with that interpreter to confirm the wheel includes its resources.

## Citation metadata

After Zenodo archives a release, update `CITATION.cff` and the README citation
from its published record. Use the released version, publication date, and
version-specific DOI. Keep the concept DOI `10.5281/zenodo.8335487` in badges
and concept identifiers. Do not invent a DOI or update a citation before the
corresponding archive exists.
