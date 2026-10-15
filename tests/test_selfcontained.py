"""astrovodka must not depend on the upstream `astromartini` distribution.

astrovodka is a fork of MARTINI: the mock-observation machinery (data cubes,
beams, spectral models, SPH kernels, noise, the `Martini` driver) was copied
into this package rather than imported from it. That decision is easy to undo
by accident -- an IDE auto-import, a copied snippet, a "just reuse the upstream
helper" refactor -- and the failure would be invisible until someone tried to
install astrovodka on a machine without astromartini.

These tests pin the decision in both directions: nothing imports `martini`, and
the modules that used to be upstream's are present and importable from here.
"""

import ast
import pathlib
import subprocess
import sys

import pytest

import astrovodka

PACKAGE_ROOT = pathlib.Path(astrovodka.__file__).parent

# Modules vendored from upstream MARTINI. If one of these disappears, the
# package has stopped being self-contained in the way the docs claim.
VENDORED = [
    "astrovodka.martini",
    "astrovodka.datacube",
    "astrovodka.beams",
    "astrovodka.noise",
    "astrovodka.spectral_models",
    "astrovodka.sph_kernels",
    "astrovodka.L_coords",
    "astrovodka.sources.sph_source",
    "astrovodka.sources._illustris_tools",
    "astrovodka.sources._L_align",
    "astrovodka.sources._extra_cartesian_transforms",
]

# Modules that are astrovodka's own additions, not forks of anything upstream.
OWN = [
    "astrovodka.cosmology_utils",
    "astrovodka.data_products",
    "astrovodka.sources.tng_source",
    "astrovodka.himodel",
]


def _python_files():
    return sorted(PACKAGE_ROOT.rglob("*.py"))


def test_no_module_imports_martini():
    """No `import martini` / `from martini import ...` anywhere in the package.

    Parsed rather than grepped so that mentions inside docstrings and strings
    (of which there are many, since the vendored code still documents its
    provenance) do not count.
    """
    offenders = []
    for path in _python_files():
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name == "martini" or alias.name.startswith("martini."):
                        offenders.append(f"{path.name}:{node.lineno} import {alias.name}")
            elif isinstance(node, ast.ImportFrom):
                # level > 0 is a relative import (e.g. `from .martini import ...`),
                # which is this package's own module and perfectly fine
                if node.level == 0 and node.module and (
                    node.module == "martini" or node.module.startswith("martini.")
                ):
                    offenders.append(
                        f"{path.name}:{node.lineno} from {node.module} import ..."
                    )
    assert not offenders, (
        "astrovodka must not import the upstream astromartini distribution; "
        f"found: {offenders}"
    )


@pytest.mark.parametrize("module", VENDORED + OWN)
def test_module_is_importable_from_astrovodka(module):
    __import__(module)


def test_the_public_api_does_not_leak_a_martini_dependency():
    """Importing the package must not pull `martini` into sys.modules."""
    code = (
        "import sys; import astrovodka; "
        "bad = [m for m in sys.modules if m == 'martini' or m.startswith('martini.')]; "
        "print(bad)"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=str(PACKAGE_ROOT.parent),
        check=True,
    )
    assert proc.stdout.strip() == "[]", f"importing astrovodka loaded {proc.stdout}"


def test_declared_dependencies_do_not_name_astromartini():
    """The packaging metadata must not declare astromartini as a dependency."""
    pyproject = PACKAGE_ROOT.parent / "pyproject.toml"
    if not pyproject.exists():  # pragma: no cover - source checkout only
        pytest.skip("pyproject.toml not found")
    assert "astromartini" not in pyproject.read_text().lower()


def test_himodel_is_self_contained_too():
    """The HI models must not reach back into the cube-making modules.

    `himodel` is the physics; `martini`/`datacube`/`sph_kernels` are the
    observation machinery. Keeping the dependency one-way means the physics can
    be tested and reused without a data cube anywhere in sight.
    """
    forbidden = {
        "astrovodka.martini",
        "astrovodka.datacube",
        "astrovodka.beams",
        "astrovodka.noise",
        "astrovodka.spectral_models",
        "astrovodka.sph_kernels",
    }
    himodel = PACKAGE_ROOT / "himodel"
    offenders = []
    for path in sorted(himodel.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name in forbidden:
                        offenders.append(f"{path.name}:{node.lineno} {alias.name}")
            elif isinstance(node, ast.ImportFrom) and node.module in forbidden:
                offenders.append(f"{path.name}:{node.lineno} {node.module}")
    assert not offenders, f"astrovodka.himodel imported the observation layer: {offenders}"
