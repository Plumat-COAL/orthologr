"""Tests of release_tools.py (run by release.yml and bundle.yml before they use it):
python3 -I -m unittest discover -s .github/scripts -p 'test_*.py'"""
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


def r_tarball(path, package, version, namespace="export(orthologr_sitrep)\n", imports="Rcpp"):
    desc = f"Package: {package}\nVersion: {version}\nImports: {imports}\nLinkingTo: Rcpp\n"
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


class Version(unittest.TestCase):
    def repo(self, version="0.6.0", heading="## `orthologr` 0.6.0"):
        d = tmpdir(self)
        write(os.path.join(d, "DESCRIPTION"), DESC.format(v=version))
        write(os.path.join(d, "NEWS.md"), f"{heading}\n\n### Section\n\n## `orthologr` 0.5.0\n")
        return d

    def test_tag_description_news_agree(self):
        self.assertEqual(rt.check_version("v0.6.0", self.repo()), ("0.6.0", "v0.6.0"))
        self.assertEqual(rt.check_version("v0.6.0", self.repo(heading="# orthologr 0.6.0 (2026-11-20)"))[0],
                         "0.6.0")

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
        for tag in ("0.6.0", "v0.6.0rc1", "bundle-2026.11", "v0.6.0 ", "v0..6", "v0.6.0.1.2"):
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
                     ("v0.6.1", "v0.6.0", "bundle-$(id)"), ("v0.6.1", "v0.6.0", "bundle-a/b")):
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

    def test_bad_sums_lines(self):
        fake_release(self.rel_dir, "v1.0.0", self.files, self.json, sums_extra=["0" * 64 + "  ../escape"])
        with self.assertRaisesRegex(rt.Failure, "not '<sha256>"):
            rt.verify_release(self.rel_dir, self.json, "v1.0.0")


class Bundle(unittest.TestCase):
    def make(self, coalr=True, sitrep=True, orthologr_version="0.6.0"):
        d = tmpdir(self)
        cdir, odir = os.path.join(d, "coal"), os.path.join(d, "orthologr")
        os.makedirs(cdir)
        os.makedirs(odir)
        if coalr:
            r_tarball(os.path.join(cdir, "coalr_0.6.1.tar.gz"), "coalr", "0.6.1", "useDynLib(coalr)\n")
        r_tarball(os.path.join(odir, "orthologr_0.6.0.tar.gz"), "orthologr", orthologr_version,
                  "export(orthologr_sitrep)\n" if sitrep else "export(dNdS)\n",
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

    def test_tarball_of_another_version(self):
        with self.assertRaisesRegex(rt.Failure, "DESCRIPTION says orthologr 0.5.9"):
            rt.bundle(self.make(orthologr_version="0.5.9"))


if __name__ == "__main__":
    unittest.main()
