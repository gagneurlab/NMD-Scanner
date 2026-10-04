# Contributing to `nmd_scanner`

## Development environment

The project uses [`uv`](https://docs.astral.sh/uv/) for dependency management.

### Initial setup

```bash
# Create the conda env (provides Python + uv)
micromamba env create -f environment-dev.yml
micromamba activate nmd_scanner

# Install all dependency groups (runtime + dev + test + lint) into a uv-managed venv
uv sync --all-groups --all-extras
```

All subsequent commands assume the env is activated and `uv` is on `PATH`.

## Running tasks

The project standardises tasks through [`tox`](https://tox.wiki/) with the `tox-uv` runner. Run any environment with:

```bash
uv run tox -e <env>
```

Available environments (defined in `pyproject.toml`):

| Env             | Purpose                                  |
|-----------------|------------------------------------------|
| `format-check`  | `ruff format --check .`                  |
| `lints`         | `ruff check .`                           |
| `typecheck`     | `mypy src/nmd_scanner`                   |
| `py3.12`        | Run pytest under Python 3.12             |
| `py3.14`        | Run pytest under Python 3.14             |

Run the full matrix CI runs with:

```bash
uv run tox
```

Tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one
process, which a debugger needs and which restores per-test output order.

### Quick commands

```bash
# Format code
uv run ruff format .

# Lint with autofix
uv run ruff check --fix .

# Run tests directly (single Python)
uv run pytest

# Run a single test
uv run pytest tests/test_rules.py::test_name -x
```

## Releases

Versioning and tagging are automated by [release-please](https://github.com/googleapis/release-please) (`.github/workflows/release-please.yml`). Publishing is handled by `.github/workflows/publish.yml`:

- release-please watches commits on `main` and opens/maintains a release PR that bumps `pyproject.toml` and updates `CHANGELOG.md`.
- Merging the release PR cuts a `vX.Y.Z` tag and a GitHub release.
- `release-please.yml` then starts `publish.yml` with `workflow_dispatch` for the new tag. `publish.yml` builds sdist + wheel and uploads to PyPI via trusted publishing (OIDC). It can also be triggered by hand from the Actions tab.
- `release-please.yml` does not call `publish.yml` via `workflow_call`, because PyPI then rejects the upload's attestation: the attestation names `release-please.yml`, but the trusted publisher on PyPI is `publish.yml`.

Use [Conventional Commits](https://www.conventionalcommits.org/) on `main` so release-please can pick the next version (`fix:` → patch, `feat:` → minor, `feat!:` / `BREAKING CHANGE:` → major). Before 1.0.0, a breaking change bumps the minor version instead (`bump-minor-pre-major` in `release-please-config.json`).

