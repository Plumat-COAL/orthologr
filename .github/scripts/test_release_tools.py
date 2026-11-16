"""Tests of release_tools.py (run by release.yml and bundle.yml before they use it):
python3 -I -m unittest discover -s .github/scripts -p 'test_*.py'"""
import contextlib
import hashlib
import io
import json
import os
import shutil
import tarfile
import tempfile
import unittest

import release_tools as rt

DESC = """Package: orthologr
Type: Package
Title: Comparative Genomics with R
Version: {v}
Depends:
    R (>= 4.4.0)
Imports:
    data.table (>= 1.9.4),
    parallel,
    Biostrings (>= 2.72.0),
    pwalign
LinkingTo: Rcpp
"""


def tmpdir(test):
    d = tempfile.mkdtemp()
    test.addCleanup(shutil.rmtree, d, ignore_errors=True)
    return d


def read(path, mode="r"):
    with open(path, mode) as fh:
        return fh.read()


def write(path, text, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode) as fh:
        fh.write(text)


def r_tarball(path, package, version, namespace="export(orthologr_sitrep)\n", imports="Rcpp", package_field=None):
    desc = f"Package: {package_field or package}\nVersion: {version}\nImports: {imports}\nLinkingTo: Rcpp\n"
    with tarfile.open(path, "w:gz") as tf:
        for name, text in (("DESCRIPTION", desc), ("NAMESPACE", namespace)):
            data = text.encode()
            info = tarfile.TarInfo(f"{package}/{name}")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))


def fake_release(directory, tag, files, json_path, draft=False, prerelease=False, sums_extra=None):
    """files {name: bytes}: writes them, SHA256SUMS.txt and the API JSON; returns the JSON dict."""
    os.makedirs(directory, exist_ok=True)
    lines = []
    for name, data in sorted(files.items()):
        write(os.path.join(directory, name), data, "wb")
        lines.append(f"{hashlib.sha256(data).hexdigest()}  {name}")
    lines += sums_extra or []
    write(os.path.join(directory, rt.SUMS), "\n".join(lines) + "\n")
    assets = []
    for name in sorted(os.listdir(directory)):
        p = os.path.join(directory, name)
        assets.append({"name": name, "size": os.path.getsize(p), "digest": "sha256:" + rt.sha256(p)})
    rel = {"tag_name": tag, "draft": draft, "prerelease": prerelease, "assets": assets,
           "html_url": f"https://example.invalid/{tag}"}
    write(json_path, json.dumps(rel))
    return rel


def api_json(directory, tag, json_path, digests=True):
    """The API JSON of a release from the files of a folder as they are now (digest left out with digests=False)."""
    assets = []
    for name in sorted(os.listdir(directory)):
        p = os.path.join(directory, name)
        a = {"name": name, "size": os.path.getsize(p)}
        if digests:
            a["digest"] = "sha256:" + rt.sha256(p)
        assets.append(a)
    write(json_path, json.dumps({"tag_name": tag, "draft": False, "prerelease": False, "assets": assets}))


class Version(unittest.TestCase):
    def repo(self, version="0.6.0", heading="## `orthologr` 0.6.0", remotes=None, package="orthologr", extra="",
             sep=",\n    "):
        d = tmpdir(self)
        desc = DESC.format(v=version).replace("Package: orthologr", f"Package: {package}")
        if remotes is not None:
            desc += "Remotes:\n    " + sep.join(remotes) + "\n"
        desc += extra
        write(os.path.join(d, "DESCRIPTION"), desc)
        write(os.path.join(d, "NEWS.md"), f"{heading}\n\n### Section\n\n## `orthologr` 0.5.0\n")
        return d

    def test_tag_description_news_agree(self):
        self.assertEqual(rt.check_version("v0.6.0", self.repo()), ("0.6.0", "v0.6.0"))
        self.assertEqual(rt.check_version("v0.6.0", self.repo(heading="# orthologr 0.6.0 (2026-11-20)"))[0],
                         "0.6.0")

    def test_another_package(self):
        with self.assertRaisesRegex(rt.Failure, "is the package 'orthologrx', not orthologr"):
            rt.check_version("v0.6.0", self.repo(package="orthologrx"))
        with self.assertRaisesRegex(rt.Failure, "not orthologr"):
            rt.check_version("", self.repo(package="orthologrx"), dry_run=True)

    SHA = "0123456789abcdef0123456789abcdef01234567"
    COALR = "coalr=Plumat-COAL/COAL/packaging/r/coalr@v0.6.1"

    def refused(self, msg=None, **kw):
        """A release refuses the DESCRIPTION; a dry run only warns ("::warning::")."""
        with self.assertRaisesRegex(rt.Failure, "not pinned", msg=msg):
            rt.check_version("v0.6.0", self.repo(**kw))
        err = io.StringIO()
        with contextlib.redirect_stderr(err):   # "::warning::" in a workflow log would be an annotation
            self.assertEqual(rt.check_version("", self.repo(**kw), dry_run=True)[0], "0.6.0")
        self.assertIn("::warning::", err.getvalue())

    def test_remotes_pinned(self):
        ok = [f"drostlab/metablastr@{self.SHA}", self.COALR, f"github::owner/repo@{self.SHA}",
              "Plumat-COAL/orthologr@v0.6.0", "github::plumat-coal/COAL@v0.6.1",
              f"git::https://github.com/owner/repo.git@{self.SHA}"]
        self.assertEqual(rt.check_version("v0.6.0", self.repo(remotes=ok))[0], "0.6.0")
        self.assertEqual(rt.check_version("v0.6.0", self.repo(remotes=[]))[0], "0.6.0")

    def test_remotes_not_pinned(self):
        sha = self.SHA
        for bad in ("drostlab/metablastr", "drostlab/metablastr@main", "drostlab/metablastr@v0.6.0rc1",
                    f"drostlab/metablastr@{sha[:12]}", f"drostlab/metablastr@{sha.upper()}",
                    "drostlab/metablastr#12", f"drostlab/metablastr#{sha}", f"drostlab/metablastr@{sha}#12",
                    "coalr=Plumat-COAL/COAL@0.6.1", "coalr=Plumat-COAL/COAL#v0.6.1",
                    "git::https://github.com/owner/repo.git", "url::https://example.org/pkg.tar.gz"):
            self.refused(msg=bad, remotes=[self.COALR, bad])

    def test_tag_pins_only_for_plumat_coal(self):
        for bad in ("drostlab/metablastr@v1.0.0", "github::drostlab/metablastr@v1.0.0",
                    "coalr=someone/COAL/packaging/r/coalr@v0.6.1", "Plumat-COAL-fork/COAL@v0.6.1",
                    "git::https://github.com/Plumat-COAL/COAL.git@v0.6.1", "x/Plumat-COAL/COAL@v0.6.1"):
            self.refused(msg=bad, remotes=[bad])

    def test_missing_comma(self):
        # without the comma the two lines are one entry "a b": refused whichever comes first
        for remotes in (["drostlab/metablastr", self.COALR], [self.COALR, "drostlab/metablastr"],
                        ["drostlab/metablastr", f"owner/repo@{self.SHA}"],
                        [f"drostlab/metablastr@{self.SHA}", self.COALR]):
            self.refused(msg=str(remotes), remotes=remotes, sep="\n    ")

    def test_config_needs(self):
        ok = ("Config/Needs/website: any::pkgdown, tidyverse/tidytemplate@" + self.SHA + "\n"
              "Config/Needs/check: rcmdcheck, cran::xml2, bioc::Biostrings, knitr@>=1.6\n")
        self.assertEqual(rt.check_version("v0.6.0", self.repo(extra=ok))[0], "0.6.0")
        for bad in ("github::tidyverse/tidytemplate", "tidyverse/tidytemplate@main", "git::https://x.org/r.git",
                    "github::owner/repo@v1.0.0", "url::https://example.org/p.tar.gz", "local::."):
            self.refused(msg=bad, extra=f"Config/Needs/website: any::pkgdown,\n    {bad}\n")

    def test_additional_repositories(self):
        self.refused(extra="Additional_repositories: https://example.org/drat\n")

    def test_tag_differs_from_description(self):
        with self.assertRaisesRegex(rt.Failure, "DESCRIPTION says Version: 0.6.0"):
            rt.check_version("v0.6.1", self.repo())

    def test_news_heading_must_name_the_version(self):
        with self.assertRaisesRegex(rt.Failure, "first heading of NEWS.md"):
            rt.check_version("v0.6.0", self.repo(heading="## `orthologr` 0.5.0"))
        for heading in ("## orthologr 0.6.0.9000", "## orthologr 10.6.0", "## orthologr 0.6.01"):
            with self.assertRaises(rt.Failure, msg=heading):
                rt.check_version("v0.6.0", self.repo(heading=heading))

    def test_bad_tags(self):
        for tag in ("0.6.0", "v0.6.0rc1", "bundle-2026.11", "v0.6.0 ", "v0..6", "v0.6.0.1.2", "v0.6.0\n",
                    "v0.6.0\nversion=9.9.9"):
            with self.assertRaises(rt.Failure, msg=tag):
                rt.check_version(tag, self.repo())

    def test_development_version(self):
        d = self.repo(version="0.6.0.9000", heading="## orthologr 0.6.0.9000")
        with self.assertRaisesRegex(rt.Failure, "development version"):
            rt.check_version("v0.6.0.9000", d)
        self.assertEqual(rt.check_version("", d, dry_run=True), ("0.6.0.9000", "v0.6.0.9000"))

    def test_dry_run_without_tag_still_checks_news(self):
        self.assertEqual(rt.check_version("", self.repo(), dry_run=True), ("0.6.0", "v0.6.0"))
        with self.assertRaises(rt.Failure):
            rt.check_version("", self.repo(heading="## orthologr 0.5.0"), dry_run=True)
        with self.assertRaises(rt.Failure):
            rt.check_version("", self.repo())   # an empty tag is a dry run only


class Inputs(unittest.TestCase):
    def test_good(self):
        rt.check_inputs("v0.6.1", "v0.6.0", "bundle-2026.11")
        rt.check_inputs("v0.6.1", "v0.6.0", "bundle-2026.11.1_a")

    def test_bad(self):
        for args in (("0.6.1", "v0.6.0", "bundle-2026.11"), ("v0.6.1rc1", "v0.6.0", "bundle-2026.11"),
                     ("v0.6.1", "v0.6.0", "v2026.11"), ("v0.6.1", "v0.6.0", "2026.11"),
                     ("v0.6.1", "v0.6.0", "bundle-"), ("v0.6.1", "v0.6.0", "bundle-a..b"),
                     ("v0.6.1", "v0.6.0", "bundle-x.lock"), ("v0.6.1", "v0.6.0", "bundle-a b"),
                     ("v0.6.1", "v0.6.0", "bundle-$(id)"), ("v0.6.1", "v0.6.0", "bundle-a/b"),
                     ("v0.6.1\n", "v0.6.0", "bundle-2026.11"), ("v0.6.1", "v0.6.0\n", "bundle-2026.11"),
                     ("v0.6.1", "v0.6.0", "bundle-2026.11\n"), ("v0.6.1", "v0.6.0", "bundle-2026.11\nx=y")):
            with self.assertRaises(rt.Failure, msg=args):
                rt.check_inputs(*args)


class Verify(unittest.TestCase):
    def setUp(self):
        self.d = tmpdir(self)
        self.rel_dir = os.path.join(self.d, "rel")
        self.json = os.path.join(self.d, "rel.json")
        self.files = {"a.tar.gz": b"aaa", "b.whl": b"bbb"}

    def test_ok(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        self.assertEqual(set(rt.verify_release(self.rel_dir, self.json, "v1.0.0")), set(self.files))

    def test_changed_file(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        write(os.path.join(self.rel_dir, "a.tar.gz"), b"AAA", "wb")
        with self.assertRaisesRegex(rt.Failure, "checksum failure"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_github_digest_disagrees(self):
        rel = fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        rel["assets"][0]["digest"] = "sha256:" + "0" * 64
        write(self.json, json.dumps(rel))
        with self.assertRaisesRegex(rt.Failure, "GitHub records"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_unlisted_file(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        write(os.path.join(self.rel_dir, "extra.txt"), "x")
        rel = json.loads(read(self.json))
        rel["assets"].append({"name": "extra.txt", "size": 1})
        write(self.json, json.dumps(rel))
        with self.assertRaisesRegex(rt.Failure, "does not list"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_listed_but_absent(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json, sums_extra=["0" * 64 + "  gone.zip"])
        with self.assertRaisesRegex(rt.Failure, "does not have"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_download_incomplete(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        os.remove(os.path.join(self.rel_dir, "b.whl"))
        with self.assertRaisesRegex(rt.Failure, "exactly the files"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_draft_prerelease_wrong_tag(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json, prerelease=True)
        with self.assertRaisesRegex(rt.Failure, "prerelease"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        with self.assertRaisesRegex(rt.Failure, "not v2"):
            rt.verify_release(self.rel_dir, self.json, "v2.0.0")

    def wrong_sums_line(self):
        """SHA256SUMS.txt gives a.tar.gz a wrong hash; the file itself is as released."""
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        sums = os.path.join(self.rel_dir, rt.SUMS)
        text = read(sums).replace(hashlib.sha256(b"aaa").hexdigest(), hashlib.sha256(b"AAA").hexdigest())
        write(sums, text)

    def test_wrong_sums_line_github_digest_agrees(self):
        self.wrong_sums_line()
        api_json(self.rel_dir, "v1.0.0", self.json)   # GitHub's digests are those of the files as they are
        with self.assertRaisesRegex(rt.Failure, "a.tar.gz: SHA256SUMS.txt says"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_wrong_sums_line_no_github_digest(self):
        self.wrong_sums_line()
        api_json(self.rel_dir, "v1.0.0", self.json, digests=False)
        with self.assertRaisesRegex(rt.Failure, "a.tar.gz: SHA256SUMS.txt says"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_no_github_digest_is_fine_otherwise(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        api_json(self.rel_dir, "v1.0.0", self.json, digests=False)
        rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_size_mismatch(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json)
        rel = json.loads(read(self.json))
        rel["assets"][0]["size"] += 1
        del rel["assets"][0]["digest"]
        write(self.json, json.dumps(rel))
        with self.assertRaisesRegex(rt.Failure, "bytes, the download has"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_duplicate_sums_lines(self):
        line = f"{hashlib.sha256(b'aaa').hexdigest()}  a.tar.gz"
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json, sums_extra=[line])
        with self.assertRaisesRegex(rt.Failure, "listed twice"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")

    def test_patterns_are_anchored_at_the_very_end(self):
        # a trailing newline must not pass (it would reach $GITHUB_OUTPUT or a file name)
        h = hashlib.sha256(b"x").hexdigest()
        self.assertIsNotNone(rt.SUM_LINE_RE.fullmatch(f"{h}  a.tar.gz"))
        self.assertIsNone(rt.SUM_LINE_RE.fullmatch(f"{h}  a.tar.gz\n"))
        self.assertIsNone(rt.SUM_LINE_RE.fullmatch(f"{h}  a.tar.gz\r"))
        self.assertIsNone(rt.RELEASE_TAG_RE.fullmatch("v0.6.0\n"))
        self.assertIsNone(rt.BUNDLE_TAG_RE.fullmatch("bundle-2026.11\n"))

    def test_bad_sums_lines(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json, sums_extra=["0" * 64 + "  ../escape"])
        with self.assertRaisesRegex(rt.Failure, "not '<sha256>"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")


class Bundle(unittest.TestCase):
    def make(self, coalr=True, sitrep=True, orthologr_version="0.6.0", orthologr_package="orthologr"):
        d = tmpdir(self)
        cdir, odir = os.path.join(d, "coal"), os.path.join(d, "orthologr")
        os.makedirs(cdir)
        os.makedirs(odir)
        if coalr:
            r_tarball(os.path.join(cdir, "coalr_0.6.1.tar.gz"), "coalr", "0.6.1", "useDynLib(coalr)\n")
        r_tarball(os.path.join(odir, "orthologr_0.6.0.tar.gz"), "orthologr", orthologr_version,
                  package_field=orthologr_package,
                  namespace="export(orthologr_sitrep)\n" if sitrep else "export(dNdS)\n",
                  imports="data.table (>= 1.9.4), parallel, Biostrings,\n    pwalign")
        for x, extra in ((cdir, "coal_src.tar.gz"), (odir, "orthologr_0.6.0-check.log")):
            write(os.path.join(x, extra), "x")
            files = {n: read(os.path.join(x, n), "rb") for n in os.listdir(x)}
            fake_release(x, "v", files, os.path.join(d, "unused.json"))
        args = type("A", (), dict(coal_dir=cdir, coal_tag="v0.6.1", coal_repo="Plumat-COAL/COAL",
                                  orthologr_dir=odir, orthologr_tag="v0.6.0", orthologr_repo="Plumat-COAL/orthologr",
                                  bundle_tag="bundle-2026.11", out=os.path.join(d, "out", "files"),
                                  notes=os.path.join(d, "out", "notes.md")))
        return args

    def test_ok(self):
        a = self.make()
        names = rt.bundle(a)
        self.assertEqual(names, ["INSTALL.md", "coalr_0.6.1.tar.gz", "orthologr_0.6.0.tar.gz"])
        sums = rt.read_sums(os.path.join(a.out, rt.SUMS))
        self.assertEqual(sorted(sums), names)
        for n in names:
            self.assertEqual(sums[n], rt.sha256(os.path.join(a.out, n)))
        notes = read(a.notes)
        self.assertEqual(notes.rstrip("\n").splitlines()[-1],
                         "Authors: Chloe Newton (@chloe-newt) and Keanu Murati (@Keanu-MUR)")
        install = read(os.path.join(a.out, "INSTALL.md"))
        i_coalr = install.index('install.packages("coalr_0.6.1.tar.gz", repos = NULL, type = "source")')
        i_orthologr = install.index('install.packages("orthologr_0.6.0.tar.gz", repos = NULL, type = "source")')
        i_sitrep = install.index("library(orthologr)\northologr_sitrep()")
        self.assertLess(i_coalr, i_orthologr)
        self.assertLess(i_orthologr, i_sitrep)
        self.assertIn('BiocManager::install(c("Biostrings", "data.table", "pwalign", "Rcpp")', install)

    def test_no_coalr_in_the_coal_release(self):
        with self.assertRaisesRegex(rt.Failure, "has no coalr_0.6.1.tar.gz"):
            rt.bundle(self.make(coalr=False))

    def test_orthologr_without_sitrep(self):
        with self.assertRaisesRegex(rt.Failure, "orthologr_sitrep"):
            rt.bundle(self.make(sitrep=False))

    def test_tarball_of_another_package(self):
        with self.assertRaisesRegex(rt.Failure, "DESCRIPTION says orthologrx 0.6.0, not orthologr 0.6.0"):
            rt.bundle(self.make(orthologr_package="orthologrx"))

    def test_tarball_of_another_version(self):
        with self.assertRaisesRegex(rt.Failure, "DESCRIPTION says orthologr 0.5.9"):
            rt.bundle(self.make(orthologr_version="0.5.9"))


if __name__ == "__main__":
    unittest.main()
