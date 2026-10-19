# Contributing

Contributions to the code, documentation, examples, and tests are welcome.
Development targets the `main` branch.

## Report a bug or propose a change

Open an [issue](https://github.com/matthiaskoenig/porous_media/issues) to report a
bug or discuss a proposal. For bug reports, include:

- A small reproducer using bundled or synthetic data when possible.
- The expected behavior and the actual behavior, including relevant error messages.
- The Python version, the `porous_media` version, and relevant dependency versions.
- Input-format details and any external programs required to reproduce the problem.

Do not include private research data or credentials. For a larger change, discuss
the intended behavior in an issue before starting implementation.

## Prepare a contribution

Create a short-lived branch from the current `main` branch. Follow the
[development environment setup](development.md#environment) to install Python,
sync the committed lockfile with uv, and install the pre-commit hooks.

Keep each change focused. Update the user guides and source docstrings when
behavior changes. For a bug fix, first reproduce the user-visible failure, then
add regression coverage that verifies meaningful behavior.

Use bundled or synthetic data for portable tests. Tests must run without external
solvers, graphics sessions, or unpublished datasets. Research scripts under
`research/` can have additional requirements; see the
[project-specific research guidance](development.md#project-specific-research-code).
Keep research data and generated simulation outputs out of commits.

Run the [required checks](development.md#required-checks) before submitting a pull
request. Include intentional `uv.lock` changes when updating dependencies.

## Submit a pull request

Submit the pull request against `main`. Describe the problem, the resulting
behavior, and the validation performed. Link the relevant issue and explain any
external data or program requirements.

Pull requests require passing `tests`, `ruff`, `ty`, `docs`, and `build` checks,
resolved review conversations, and linear history. Contributor pull requests need
maintainer review. The solo-maintainer policy does not require approving reviews,
so the maintainer can merge their own pull requests after checks pass.

See [repository policies](development.md#repository-policies) for the versioned
ruleset definitions and maintainer instructions. Existing release tags cannot be
moved or deleted.

## Releases

Follow [release guide](release.md)
for version changes, release notes, and release tags. Do not edit generated files
or historical changelog files manually.
