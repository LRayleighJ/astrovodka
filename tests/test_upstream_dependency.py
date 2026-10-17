"""astrovodka must *call* MARTINI, not contain a copy of it.

The observation machinery (data cubes, beams, spectral models, SPH kernels,
noise, the ``Martini`` driver, ``SPHSource`` and the IllustrisTNG I/O helpers)
lives in the ``astromartini`` distribution. astrovodka imports and subclasses
it; it does not vendor it.

These tests pin that decision in both directions:

- the packaging metadata declares ``astromartini`` and importing astrovodka
  loads the real upstream package;
- the classes astrovodka re-exports either *are* MARTINI's objects, or are
  subclasses of them — nothing is a private reimplementation;
- the ``himodel`` physics layer still does not reach into the observation layer
  (or into MARTINI) at all.
"""

import ast
import pathlib
import re
import subprocess
import sys

import pytest

import astrovodka
import martini
import martini.beams
import martini.datacube
import martini.martini
import martini.noise
import martini.sources
import martini.sources._illustris_tools
import martini.spectral_models
import martini.sph_kernels

import astrovodka.beams
import astrovodka.datacube
import astrovodka.martini
import astrovodka.noise
import astrovodka.sources
import astrovodka.sources._illustris_tools
import astrovodka.spectral_models
import astrovodka.sph_kernels

PACKAGE_ROOT = pathlib.Path(astrovodka.__file__).parent

# Classes astrovodka re-exports unchanged. They must be the *same objects* as
# upstream's: a copy would silently diverge on the next MARTINI release.
PASSTHROUGH = [
    (astrovodka.beams, "GaussianBeam", martini.beams, "GaussianBeam"),
    (astrovodka.beams, "_BaseBeam", martini.beams, "_BaseBeam"),
    (astrovodka.noise, "GaussianNoise", martini.noise, "GaussianNoise"),
    (astrovodka.noise, "_BaseNoise", martini.noise, "_BaseNoise"),
    (
        astrovodka.sph_kernels,
        "CubicSplineKernel",
        martini.sph_kernels,
        "CubicSplineKernel",
    ),
    (
        astrovodka.sph_kernels,
        "DiracDeltaKernel",
        martini.sph_kernels,
        "DiracDeltaKernel",
    ),
    (astrovodka.sph_kernels, "_BaseSPHKernel", martini.sph_kernels, "_BaseSPHKernel"),
]

# Classes astrovodka specialises. They must inherit from the upstream class,
# so that MARTINI's implementation remains the implementation that runs.
SUBCLASSES = [
    (astrovodka.datacube, "DataCube", martini.datacube, "DataCube"),
    (astrovodka.martini, "Martini", martini.martini, "Martini"),
    (astrovodka.martini, "GlobalProfile", martini.martini, "GlobalProfile"),
    (astrovodka.sources, "SPHSource", martini.sources, "SPHSource"),
    (
        astrovodka.spectral_models,
        "GaussianSpectrum",
        martini.spectral_models,
        "GaussianSpectrum",
    ),
    (
        astrovodka.spectral_models,
        "DiracDeltaSpectrum",
        martini.spectral_models,
        "DiracDeltaSpectrum",
    ),
]


def _declared_dependencies():
    """Return the ``[project] dependencies`` list from pyproject.toml."""
    pyproject = PACKAGE_ROOT.parent / "pyproject.toml"
    text = pyproject.read_text(encoding="utf-8")
    try:  # Python 3.11+
        import tomllib
    except ImportError:  # pragma: no cover - exercised on Python 3.9/3.10
        block = text.split("dependencies = [", 1)[1].split("]", 1)[0]
        return re.findall(r'"([^"]+)"', block)
    return tomllib.loads(text)["project"]["dependencies"]


def test_astromartini_is_a_required_dependency():
    """``pip install astrovodka`` must install MARTINI."""
    dependencies = _declared_dependencies()
    assert any(
        re.match(r"^astromartini\s*[><=!~]", dependency) or dependency == "astromartini"
        for dependency in dependencies
    ), f"astromartini missing from project.dependencies: {dependencies}"


def test_python_floor_matches_the_supported_martini_releases():
    """The declared Python floor must be one MARTINI can be installed on."""
    pyproject = PACKAGE_ROOT.parent / "pyproject.toml"
    match = re.search(
        r'requires-python\s*=\s*">=(\d+)\.(\d+)"', pyproject.read_text(encoding="utf-8")
    )
    assert match, "requires-python not found in pyproject.toml"
    # astromartini 2.1.9 through 2.1.15 support Python 3.9; this test exists to
    # catch a floor bump that would strand the versions declared above.
    assert (int(match.group(1)), int(match.group(2))) <= (3, 9)


def test_importing_astrovodka_loads_upstream_martini():
    """``import astrovodka`` must pull in the real ``martini`` package."""
    code = (
        "import sys; import astrovodka; "
        "print(sorted(m for m in sys.modules if m == 'martini' or m.startswith('martini.')))"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=str(PACKAGE_ROOT.parent),
        check=True,
    )
    assert "martini" in proc.stdout, proc.stdout
    # Compare against the package directory, not a path component: a checkout
    # may legitimately live under a directory that is itself named astrovodka.
    martini_path = pathlib.Path(martini.__file__).resolve()
    assert PACKAGE_ROOT.resolve() not in martini_path.parents, (
        f"the imported `martini` is inside astrovodka: {martini.__file__}"
    )


@pytest.mark.parametrize("ours, name, theirs, their_name", PASSTHROUGH)
def test_re_exported_classes_are_upstream_objects(ours, name, theirs, their_name):
    assert getattr(ours, name) is getattr(theirs, their_name)


@pytest.mark.parametrize("ours, name, theirs, their_name", SUBCLASSES)
def test_specialised_classes_subclass_upstream(ours, name, theirs, their_name):
    assert issubclass(getattr(ours, name), getattr(theirs, their_name))


def test_illustris_io_is_a_re_export_not_a_copy():
    """The TNG I/O shim must contain imports only, no implementation."""
    path = PACKAGE_ROOT / "sources" / "_illustris_tools.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    definitions = [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    ]
    assert not definitions, (
        "astrovodka.sources._illustris_tools must re-export "
        "martini.sources._illustris_tools, not define its own helpers"
    )
    for name in (
        "loadHeader",
        "loadSingle",
        "loadSubset",
        "getSnapOffsets",
        "snapPath",
    ):
        assert getattr(astrovodka.sources._illustris_tools, name) is getattr(
            martini.sources._illustris_tools, name
        )


def test_no_vendored_upstream_modules_remain():
    """The modules that used to be vendored from MARTINI must be gone."""
    sources = PACKAGE_ROOT / "sources"
    for name in ("_L_align.py", "_extra_cartesian_transforms.py"):
        assert not (sources / name).exists(), f"{name} is a vendored MARTINI module"


def test_himodel_is_independent_of_the_observation_layer():
    """The HI models must not import MARTINI, directly or through astrovodka.

    ``himodel`` is the physics; the cube-making machinery is MARTINI's. Keeping
    the dependency one-way means the physics can be tested and reused without a
    data cube — or an installed MARTINI — anywhere in sight.
    """
    forbidden = {
        "martini",
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
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    root = alias.name.split(".")[0]
                    if alias.name in forbidden or root == "martini":
                        offenders.append(f"{path.name}:{node.lineno} {alias.name}")
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                if (
                    node.module in forbidden
                    or (node.module or "").split(".")[0] == "martini"
                ):
                    offenders.append(f"{path.name}:{node.lineno} {node.module}")
    assert not offenders, (
        f"astrovodka.himodel imported the observation layer: {offenders}"
    )
