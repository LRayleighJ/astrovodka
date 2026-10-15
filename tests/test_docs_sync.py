"""Tests that keep the documentation trees honest.

Two failure modes are guarded here, both of which the repository has actually
suffered from:

1. ``docs/html/`` drifting away from the Markdown it is generated from. The two
   used to be maintained by hand in parallel, and ended up disagreeing about
   both the version number and a parameter's default value.
2. The English and Chinese trees falling out of step, so that a page exists in
   one language and not the other.

No TNG data, no network, and only the standard library plus pytest.
"""

import importlib.util
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS = REPO_ROOT / "docs"
LANGS = ("en", "zh")


def _load_build_html():
    """Import ``docs/build_html.py``, which is not part of the package."""
    path = DOCS / "build_html.py"
    spec = importlib.util.spec_from_file_location("_build_html", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


build_html = _load_build_html()


def _markdown_pages(lang):
    """Every Markdown page under ``docs/<lang>/``, as paths relative to it."""
    return sorted(p.relative_to(DOCS / lang) for p in (DOCS / lang).rglob("*.md"))


# --------------------------------------------------------------------------
# docs/html/ is in sync with the Markdown
# --------------------------------------------------------------------------


def test_html_docs_are_up_to_date():
    """Regenerating must be a no-op: committed HTML matches the Markdown."""
    assert build_html.build(check_only=True) == 0, (
        "docs/html/ is stale. Run `python docs/build_html.py` and commit the result."
    )


def test_html_output_is_deterministic():
    """Two builds of the same source produce identical bytes.

    Without this, the sync check above would be flaky, and the CI job that
    regenerates and diffs would fail for no real reason.
    """
    pages = sorted((DOCS / "en").rglob("*.md"))
    assert pages, "expected at least one English Markdown page"

    for page in pages:
        source = page.read_text(encoding="utf-8")
        lines = source.splitlines()
        known = {p.resolve() for p in DOCS.glob("*/*.md")}
        known |= {p.resolve() for p in DOCS.glob("*/*/*.md")}
        link_map = build_html.make_link_map(page, "en", known)
        first = build_html.render_blocks(lines, link_map, [])
        second = build_html.render_blocks(lines, link_map, [])
        assert first == second, f"non-deterministic rendering of {page}"


def test_every_markdown_page_has_html():
    for lang in LANGS:
        for rel in _markdown_pages(lang):
            html = DOCS / "html" / lang / rel.with_suffix(".html")
            assert html.is_file(), f"missing generated page for docs/{lang}/{rel}"


def test_generated_html_is_valid_json_free_document():
    """Each generated page is a complete, standalone HTML document."""
    for lang in LANGS:
        for html in sorted((DOCS / "html" / lang).rglob("*.html")):
            text = html.read_text(encoding="utf-8")
            assert text.startswith("<!DOCTYPE html>"), html
            assert text.rstrip().endswith("</html>"), html
            assert f'<html lang="{lang}">' in text, html
            # The wrong lang attribute shipped on all 14 pages of the old
            # hand-written HTML tree; make that unrepresentable.
            assert 'lang="zh-CN"' not in text, html
            assert "<title>" in text, html


# --------------------------------------------------------------------------
# The English and Chinese trees are parallel
# --------------------------------------------------------------------------


def test_language_trees_are_parallel():
    """A page present in one language must be present in the other."""
    en, zh = set(_markdown_pages("en")), set(_markdown_pages("zh"))

    assert not en - zh, f"missing Chinese counterparts: {sorted(en - zh)}"
    assert not zh - en, f"missing English counterparts: {sorted(zh - en)}"


def test_api_reference_is_complete_in_both_languages():
    """Every public module has an API page, in both languages.

    Modules whose names begin with an underscore are internals, and whether to
    document them is the documentation's choice — ``_illustris_tools`` has a
    page because users reach for it, ``_L_align`` and
    ``_extra_cartesian_transforms`` do not. What must not happen is a *public*
    module going undocumented, or a page existing in one language only.
    """
    # Pages are per-module, plus one for the `himodel` subpackage, which is a
    # body of work in its own right.
    modules = {p.stem for p in (REPO_ROOT / "astrovodka").glob("*.py")}
    modules |= {p.stem for p in (REPO_ROOT / "astrovodka" / "sources").glob("*.py")}
    modules.add("himodel")

    # `sources` is a namespace for its two classes, both of which have their
    # own pages; it does not need one of its own.
    exempt = {"__init__", "__version__", "sources"}

    public = {
        name
        for name in modules
        if not name.startswith("_") and name not in exempt
    }

    for lang in LANGS:
        documented = {p.stem for p in (DOCS / lang / "api").glob("*.md")}
        documented.discard("index")
        missing = public - documented
        assert not missing, f"undocumented public modules in docs/{lang}/api: {sorted(missing)}"
        # A page for a private module is optional, but it has to exist in both
        # languages like any other page.
        assert not (documented - public - modules), (
            f"docs/{lang}/api has pages for modules that do not exist: "
            f"{sorted(documented - public - modules)}"
        )


def test_both_readmes_exist_and_cross_link():
    en = REPO_ROOT / "README.md"
    zh = REPO_ROOT / "README.zh-CN.md"
    assert en.is_file() and zh.is_file()

    assert "README.zh-CN.md" in en.read_text(encoding="utf-8"), (
        "the English README must link to the Chinese one"
    )
    assert "README.md" in zh.read_text(encoding="utf-8"), (
        "the Chinese README must link back to the English one"
    )


# --------------------------------------------------------------------------
# No stale cross-references
# --------------------------------------------------------------------------


def test_no_references_to_out_of_repo_report():
    """`report.md` used to live outside the repository and was cited as evidence.

    It is now docs/en/verification.md; nothing should still point at the old
    path, since a reader cannot follow it.
    """
    offenders = []
    for path in list(REPO_ROOT.glob("*.md")) + list(DOCS.rglob("*.md")):
        if "report.md" in path.read_text(encoding="utf-8"):
            offenders.append(str(path.relative_to(REPO_ROOT)))
    assert not offenders, f"stale `report.md` references in: {offenders}"


def test_no_hardcoded_machine_paths_in_docs():
    """Docs must not name one machine's directories."""
    pattern = re.compile(r"/home/[a-z]|/data/hw2")
    offenders = []
    for path in list(REPO_ROOT.glob("*.md")) + list(DOCS.rglob("*.md")):
        text = path.read_text(encoding="utf-8")
        if pattern.search(text):
            offenders.append(str(path.relative_to(REPO_ROOT)))
    assert not offenders, f"machine-specific paths in: {offenders}"


@pytest.mark.parametrize("lang", LANGS)
def test_relative_links_resolve(lang):
    """Every relative Markdown link inside a language tree points at a real file.

    Links that escape the docs tree are exempt: the HTML generator rewrites
    those to repository URLs, since a relative path cannot be resolved from the
    rendered output.
    """
    broken = []
    for rel in _markdown_pages(lang):
        page = DOCS / lang / rel
        for target in re.findall(r"\]\(([^)\s]+)\)", page.read_text(encoding="utf-8")):
            if re.match(r"^(https?:|mailto:|#)", target) or target.startswith("$"):
                continue
            # A fragment (`page.md#section`) still has to point at a real file.
            path_part = target.split("#", 1)[0]
            if not path_part:
                continue
            resolved = (page.parent / path_part).resolve()
            if not resolved.is_relative_to(DOCS):
                continue  # rendered as a repository URL instead
            if not resolved.exists():
                broken.append(f"docs/{lang}/{rel} -> {target}")
    assert not broken, "broken relative links:\n  " + "\n  ".join(broken)
