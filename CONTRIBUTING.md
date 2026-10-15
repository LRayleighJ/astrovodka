# Contributing to astrovodka

Thanks for considering a contribution. This page covers the development setup,
the conventions the project follows, and the checks a change is expected to
pass.

- [Setting up a development checkout](#setting-up-a-development-checkout)
- [Running the tests](#running-the-tests)
- [What the tests do not cover](#what-the-tests-do-not-cover)
- [Code style](#code-style)
- [Documentation](#documentation)
- [A note on bit-exactness](#a-note-on-bit-exactness)
- [Reporting a bug](#reporting-a-bug)
- [Pull requests](#pull-requests)
- [Licence and provenance](#licence-and-provenance)

---

## Setting up a development checkout

```bash
git clone https://github.com/LRayleighJ/astrovodka.git
cd astrovodka
conda env create -f environment.yml
conda activate astrovodka
```

or, without conda:

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install -e ".[dev]"
```

`[dev]` brings in h5py, matplotlib, multiprocess, pytest and pytest-cov — that
is, everything needed for the test suite, the examples and the tutorials.

**TNG data is not required for development.** The suite is entirely synthetic;
you only need a snapshot if you want to run the tutorial notebooks end to end.

## Running the tests

```bash
python -m pytest tests/ -q                     # everything, ~15 s
python -m pytest tests/test_hi_models.py -v    # one file
python -m pytest tests/ -q --cov=astrovodka    # with coverage
```

`pythonpath` is configured in `pyproject.toml`, so plain `pytest` works from a
fresh clone without installing anything.

The suite runs on Python 3.9 through 3.13 in CI. Please make sure it passes on
at least the extremes if you change anything version-sensitive.

## What the tests do not cover

Two categories of behaviour are verified outside `tests/`, and a change that
touches them needs more than a green test run:

- **The demo subhalo numbers** (HI and H₂ masses, `f_gas_in_map`, the UV field
  range) are checked by the tutorial notebook, which asserts on its own
  results as it runs. If you change the physics, expect to update
  `docs/en/verification.md` and its Chinese counterpart.
- **Bit-exactness against the reference implementation.** See
  [below](#a-note-on-bit-exactness).

## Code style

| | |
|---|---|
| Formatter / linter | `ruff`, line length 88, target Python 3.9 |
| Docstrings | numpydoc convention |
| Quotes | double |
| Type annotations | expected on public functions |
| Test files | exempt from annotation requirements |

```bash
ruff check .
ruff format --check --diff .
```

The ruff configuration is deliberately narrow — `E4`, `E7`, `E9`, `F` — because
the observation machinery (`astrovodka/martini.py`, `datacube.py`, `beams.py`,
`noise.py`, `spectral_models.py`, `sph_kernels.py`) is vendored from MARTINI and
is not restyled here. Restyling it would make future comparisons with upstream
much harder for no benefit.

## Documentation

The Markdown under `docs/en/` and `docs/zh/` is the source of truth. The HTML in
`docs/html/` is **generated**:

```bash
python docs/build_html.py           # regenerate
python docs/build_html.py --check   # fail if stale
```

If you edit or add a Markdown page you must regenerate the HTML and commit it,
or CI will fail. The two trees are kept parallel: a page that exists in one
language must exist in the other, and `tests/test_docs_sync.py` enforces that,
along with checking that relative links resolve and that no page names a
machine-specific path.

When adding a public module, add an API page for it in **both** trees —
`test_api_reference_is_complete_in_both_languages` will tell you if you forget.

## A note on bit-exactness

The HI/H₂ models in `astrovodka/himodel/` are ported from a reference
implementation, and the port is verified to agree with it **bit for bit** — see
`docs/en/verification.md`. That agreement is a property worth protecting, and it
is fragile in a specific way: it rests on the arithmetic being performed in the
same order as the original.

Two consequences for contributors:

- **Do not reorder operations** in the physics for readability. `x / h` is not
  `x * (1 / h)`; the two differ in the last bit, and the difference shows up in
  the verification table.
- **Do not replace the astropy unit chain** in `rho_code_to_physical` with a
  hard-coded constant. astropy's kpc and the module's own `kpc_to_cm` differ by
  about 8 × 10⁻⁷ relative, which is a silent change at the 10⁻⁶ level in every
  density the model computes.

If a change genuinely requires breaking exactness, say so explicitly in the pull
request and update the verification record.

## Reporting a bug

Please include:

- the astrovodka version (`python -c "import astrovodka; print(astrovodka.__version__)"`),
  the Python version, and the versions of numpy, scipy and astropy;
- what you did, what you expected, and what happened;
- the full traceback, if there is one.

For anything involving TNG data, please also say which simulation and snapshot,
and the subhalo ID. **Do not attach the snapshot files.**

If your cubes look wrong rather than crashing, `docs/en/guide.md` has a
[Pitfalls](docs/en/guide.md#pitfalls) section with the four failure modes that
account for most of these — a pruned source, dropped padding, a shared spectral
model, and a spectral axis too narrow for the line. It is worth ruling those out
first.

## Pull requests

- Keep the change focused; one topic per pull request.
- Add a test for a bug fix. The regressions fixed in 0.2.0 and 0.3.0 each have
  one, and the two most instructive are described in
  `docs/en/verification.md` — both were written specifically so that the
  *obvious* formulation of the test would not have caught the bug.
- Update `CHANGELOG.md` under `## [Unreleased]`.
- If you add a public name, export it in the relevant `__all__` and document it
  in both language trees.
- Run `ruff check .` and `python -m pytest tests/ -q` before pushing.

## Licence and provenance

astrovodka is distributed under the **GNU GPL v3 or later**, inherited from
MARTINI, of which it is a derivative work. By contributing you agree that your
contribution is licensed under the same terms.

Two rules follow from this, and both are enforced by tests:

1. **No runtime dependency on `astromartini`.** The observation machinery is
   vendored, not imported. `tests/test_selfcontained.py` parses the package with
   `ast` and fails on any absolute `import martini`.
2. **`astrovodka.himodel` must not import the observation layer.** The physics
   has to remain usable with no data cube anywhere in sight.

If you port code from elsewhere, say where it came from in the pull request and
add it to `NOTICE`. If you port physics from a paper, cite the paper in the
module docstring and add it to `CITATION.cff`.
