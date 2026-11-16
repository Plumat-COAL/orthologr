"""Checks and files for the two release workflows of this fork (.github/workflows/release.yml, bundle.yml).

Standard library only; run it with `python3 -I` from the repository root. Every subcommand exits 1 with a message
starting "release_tools:" when a check fails, 0 otherwise.

  version TAG [--dry-run] [--github-output]
      release.yml: the tag must be v<Version of DESCRIPTION> and the first heading of NEWS.md must name that version.
      A development version (a component >= 9000, R's convention) is no release. With --dry-run an empty TAG means
      "the version DESCRIPTION has" and a development version is allowed; NEWS.md is still checked.
      Prints the version; --github-output writes version=, tag= to $GITHUB_OUTPUT.

  check-inputs COAL_TAG ORTHOLOGR_TAG BUNDLE_TAG
      bundle.yml: the two release tags look like v<digits>, and the bundle tag is bundle-<...>, so that it can never
      be taken for (or collide with) a release tag of either project.

  verify DIR --release-json FILE --tag TAG
      bundle.yml: DIR holds every file of a GitHub release (gh release download) and FILE is that release from the
      REST API (gh api repos/<owner>/<repo>/releases/tags/<tag>). The release must be the published, non-draft,
      non-prerelease release of TAG; DIR must hold exactly its files; SHA256SUMS.txt must list every other file,
      once, and every listed SHA-256 must be the file's; GitHub's own digest of each file (when it records one) must
      agree too.

  bundle --coal-dir D --coal-tag T --coal-repo R --orthologr-dir D --orthologr-tag T --orthologr-repo R
         --bundle-tag B --out DIR --notes FILE
      bundle.yml, after verify: copies coalr_<v>.tar.gz and orthologr_<v>.tar.gz into DIR (after checking the
      Package and Version fields of each tarball's DESCRIPTION, and that orthologr exports orthologr_sitrep(), which
      INSTALL.md runs), writes INSTALL.md and SHA256SUMS.txt (every file of DIR), and the release notes to FILE.
"""
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile

SUMS = "SHA256SUMS.txt"
AUTHORS_LINE = "Authors: Chloe Newton (@chloe-newt) and Keanu Murati (@Keanu-MUR)"
# v + an R package version: 2 to 4 numbers (R allows "-" too; a tag uses ".")
RELEASE_TAG_RE = re.compile(r"^v(\d+(?:\.\d+){1,3})$")
# bundle-<letters, digits, "." "_" "-">: never "v<digit>", so never a release tag of COAL or orthologr
BUNDLE_TAG_RE = re.compile(r"^bundle-[0-9A-Za-z](?:[0-9A-Za-z._-]{0,62}[0-9A-Za-z])?$")
SUM_LINE_RE = re.compile(r"^([0-9a-f]{64}) [ *]([^/\\]+)$")
# packages that come with R and are not installed from CRAN or Bioconductor
BASE_PACKAGES = {"R", "base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods", "parallel",
                 "splines", "stats", "stats4", "tcltk", "tools", "utils"}


class Failure(Exception):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_dcf(text):
    """The fields of a DESCRIPTION file (continuation lines joined with a space)."""
    fields, key = {}, None
    for line in text.splitlines():
        if not line.strip():
            continue
        if line[0] in " \t" and key:
            fields[key] += " " + line.strip()
            continue
        m = re.match(r"^([A-Za-z0-9@/._-]+):\s*(.*)$", line)
        if not m:
            raise Failure(f"DESCRIPTION: cannot read the line {line!r}")
        key = m.group(1)
        fields[key] = m.group(2).strip()
    return fields


def package_names(field):
    """'Rcpp (>= 0.11.3), seqinr,\n parallel' -> ['Rcpp', 'seqinr', 'parallel']."""
    out = []
    for item in (field or "").split(","):
        name = re.sub(r"\(.*?\)", "", item).strip()
        if name:
            out.append(name)
    return out


def is_development(version):
    return any(int(c) >= 9000 for c in re.split(r"[.-]", version)[2:])


def first_heading(news_text):
    for line in news_text.splitlines():
        if line.startswith("#"):
            return line
    return None


def heading_names(heading, version):
    """True when the heading has the version as a whole token ("0.6.0" but not "0.6.0.9000" or "10.6.0")."""
    return re.search(r"(?<![0-9.])" + re.escape(version) + r"(?![0-9]|\.[0-9])", heading) is not None


# ---------------------------------------------------------------------------------------------------------------
# version (release.yml)

def check_version(tag, root=".", dry_run=False):
    with open(os.path.join(root, "DESCRIPTION"), encoding="utf-8") as fh:
        desc = parse_dcf(fh.read())
    version = desc.get("Version", "")
    if desc.get("Package") != "orthologr":
        raise Failure(f"DESCRIPTION is the package {desc.get('Package')!r}, not orthologr")
    if not re.fullmatch(r"\d+(?:[.-]\d+){1,3}", version):
        raise Failure(f"DESCRIPTION has no usable Version ({version!r})")
    if tag == "" and dry_run:
        tag = "v" + version
    m = RELEASE_TAG_RE.match(tag)
    if not m:
        raise Failure(f"tag {tag!r} is not v<version>, e.g. v0.6.0")
    if m.group(1) != version:
        raise Failure(f"tag {tag} says {m.group(1)}, but DESCRIPTION says Version: {version}: set the version, "
                      f"commit, then tag that commit")
    if is_development(version) and not dry_run:
        raise Failure(f"{version} is a development version (a component >= 9000): a release sets X.Y.Z first")
    with open(os.path.join(root, "NEWS.md"), encoding="utf-8") as fh:
        heading = first_heading(fh.read())
    if heading is None:
        raise Failure("NEWS.md has no heading")
    if not heading_names(heading, version):
        raise Failure(f"the first heading of NEWS.md is {heading.strip()!r}: it must name {version}, the version "
                      f"being released")
    return version, tag


# ---------------------------------------------------------------------------------------------------------------
# check-inputs (bundle.yml)

def check_inputs(coal_tag, orthologr_tag, bundle_tag):
    for what, tag in (("coal_tag", coal_tag), ("orthologr_tag", orthologr_tag)):
        if not RELEASE_TAG_RE.match(tag):
            raise Failure(f"{what} {tag!r} is not the tag of a release, v<version> (e.g. v0.6.1); a prerelease "
                          f"(rc) cannot go into a bundle")
    if not BUNDLE_TAG_RE.match(bundle_tag) or ".." in bundle_tag or bundle_tag.endswith(".lock"):
        raise Failure(f"bundle_tag {bundle_tag!r} is not bundle-<name> (letters, digits, '.', '_', '-'; e.g. "
                      f"bundle-2026.11): a bundle tag must never look like a release tag of COAL or orthologr")


# ---------------------------------------------------------------------------------------------------------------
# verify (bundle.yml)

def read_sums(path):
    sums = {}
    with open(path, encoding="utf-8") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            m = SUM_LINE_RE.match(line)
            if not m:
                raise Failure(f"{path}:{n}: not '<sha256>  <file name>': {line!r}")
            name = m.group(2)
            if name in (".", "..") or name in sums:
                raise Failure(f"{path}:{n}: {name!r} is not a file name, or is listed twice")
            sums[name] = m.group(1)
    if not sums:
        raise Failure(f"{path} lists no files")
    return sums


def verify_release(directory, release_json, tag):
    with open(release_json, encoding="utf-8") as fh:
        rel = json.load(fh)
    where = rel.get("html_url", tag)
    if rel.get("tag_name") != tag:
        raise Failure(f"the release found is {rel.get('tag_name')!r}, not {tag}")
    if rel.get("draft") or rel.get("prerelease"):
        raise Failure(f"{where} is a draft or a prerelease: a bundle takes published releases only")
    assets = {a["name"]: a for a in rel.get("assets", [])}
    if SUMS not in assets:
        raise Failure(f"{where} has no {SUMS}")
    files = sorted(f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)))
    if set(files) != set(assets):
        raise Failure(f"{directory} does not hold exactly the files of {where}: missing "
                      f"{sorted(set(assets) - set(files))}, unexpected {sorted(set(files) - set(assets))}")
    sums = read_sums(os.path.join(directory, SUMS))
    unlisted = sorted(set(files) - set(sums) - {SUMS})
    if unlisted:
        raise Failure(f"{where}: {SUMS} does not list {unlisted}")
    missing = sorted(set(sums) - set(files))
    if missing:
        raise Failure(f"{where}: {SUMS} lists files the release does not have: {missing}")
    bad = []
    for name in files:
        got = sha256(os.path.join(directory, name))
        if name in sums and got != sums[name]:
            bad.append(f"{name}: {SUMS} says {sums[name]}, the file is {got}")
        digest = assets[name].get("digest")
        if digest and digest != "sha256:" + got:
            bad.append(f"{name}: GitHub records {digest}, the downloaded file is sha256:{got}")
        if assets[name].get("size") not in (None, os.path.getsize(os.path.join(directory, name))):
            bad.append(f"{name}: GitHub records {assets[name]['size']} bytes, the download has "
                       f"{os.path.getsize(os.path.join(directory, name))}")
    if bad:
        raise Failure(f"{where}: checksum failure:\n  " + "\n  ".join(bad))
    return sums


# ---------------------------------------------------------------------------------------------------------------
# bundle (bundle.yml)

def tarball_members(path, package, names):
    """{name: text} of <package>/<name> in an R source tarball, read without extracting anything."""
    out = {}
    with tarfile.open(path, "r:gz") as tf:
        for name in names:
            try:
                member = tf.getmember(f"{package}/{name}")
            except KeyError:
                raise Failure(f"{os.path.basename(path)} has no {package}/{name}: not an R source package "
                              f"{package}") from None
            if not member.isfile() or member.size > 1 << 20:
                raise Failure(f"{os.path.basename(path)}: {package}/{name} is not a small regular file")
            out[name] = tf.extractfile(member).read().decode("utf-8", errors="replace")
    return out


def check_tarball(path, package, version):
    files = tarball_members(path, package, ["DESCRIPTION", "NAMESPACE"])
    desc = parse_dcf(files["DESCRIPTION"])
    if desc.get("Package") != package or desc.get("Version") != version:
        raise Failure(f"{os.path.basename(path)}: DESCRIPTION says {desc.get('Package')} {desc.get('Version')}, "
                      f"not {package} {version}")
    deps = []
    for field in ("Depends", "Imports", "LinkingTo"):
        deps += [p for p in package_names(desc.get(field)) if p not in BASE_PACKAGES]
    return desc, files["NAMESPACE"], deps


def install_md(coal_tag, coal_repo, coalr_file, orthologr_tag, orthologr_repo, orthologr_file, bundle_tag, deps):
    dep_list = ", ".join(f'"{d}"' for d in deps)
    return f"""# COAL {coal_tag} + orthologr {orthologr_tag} ({bundle_tag})

| File | What | From |
|---|---|---|
| `{coalr_file}` | the R package coalr (COAL's engine for R) | https://github.com/{coal_repo}/releases/tag/{coal_tag} |
| `{orthologr_file}` | the R package orthologr | https://github.com/{orthologr_repo}/releases/tag/{orthologr_tag} |

Each tarball is byte-identical to the file of the same name in its project's release and has the SHA-256 that
release's `SHA256SUMS.txt` gives for it.

1. Download the files of this release into one folder and check them there:
   `sha256sum -c SHA256SUMS.txt` (macOS: `shasum -a 256 -c SHA256SUMS.txt`).
2. Both packages are installed from source, so a C++ compiler is needed (Windows: Rtools; macOS: the Xcode
   command line tools). In R, with that folder as the working directory:

```r
install.packages("BiocManager")
BiocManager::install(c({dep_list}), update = FALSE, ask = FALSE)
# coalr first: orthologr uses it when it is installed
install.packages("{coalr_file}", repos = NULL, type = "source")
install.packages("{orthologr_file}", repos = NULL, type = "source")
library(orthologr)
orthologr_sitrep()
```

`orthologr_sitrep()` reports what orthologr finds (coalr, DIAMOND, BLAST+ and the other external programs) and
what is missing.
"""


def bundle(args):
    coal_v = RELEASE_TAG_RE.match(args.coal_tag).group(1)
    orthologr_v = RELEASE_TAG_RE.match(args.orthologr_tag).group(1)
    coalr_file = f"coalr_{coal_v}.tar.gz"
    orthologr_file = f"orthologr_{orthologr_v}.tar.gz"
    coalr_path = os.path.join(args.coal_dir, coalr_file)
    orthologr_path = os.path.join(args.orthologr_dir, orthologr_file)
    if not os.path.isfile(coalr_path):
        raise Failure(f"the COAL release {args.coal_tag} has no {coalr_file} (COAL attaches coalr to its releases "
                      f"from 0.6.1 on): files are {sorted(os.listdir(args.coal_dir))}")
    if not os.path.isfile(orthologr_path):
        raise Failure(f"the orthologr release {args.orthologr_tag} has no {orthologr_file}: files are "
                      f"{sorted(os.listdir(args.orthologr_dir))}")
    _, _, coalr_deps = check_tarball(coalr_path, "coalr", coal_v)
    _, namespace, orthologr_deps = check_tarball(orthologr_path, "orthologr", orthologr_v)
    if not re.search(r"\borthologr_sitrep\b", namespace):
        raise Failure(f"{orthologr_file} does not export orthologr_sitrep(), which INSTALL.md tells people to run")
    deps = sorted(set(coalr_deps + orthologr_deps) - {"coalr", "orthologr"}, key=str.lower)

    os.makedirs(args.out, exist_ok=False)
    shutil.copyfile(coalr_path, os.path.join(args.out, coalr_file))
    shutil.copyfile(orthologr_path, os.path.join(args.out, orthologr_file))
    with open(os.path.join(args.out, "INSTALL.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(install_md(args.coal_tag, args.coal_repo, coalr_file, args.orthologr_tag, args.orthologr_repo,
                            orthologr_file, args.bundle_tag, deps))
    names = sorted(os.listdir(args.out))
    lines = [f"{sha256(os.path.join(args.out, n))}  {n}" for n in names]
    # the two tarballs: the very lines of their projects' SHA256SUMS.txt (verify has checked those)
    for src_dir, name in ((args.coal_dir, coalr_file), (args.orthologr_dir, orthologr_file)):
        if read_sums(os.path.join(src_dir, SUMS))[name] != sha256(os.path.join(args.out, name)):
            raise Failure(f"{name} changed while it was copied")
    with open(os.path.join(args.out, SUMS), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")

    table = "\n".join(f"| `{n}` | `{s}` |" for s, n in (line.split("  ", 1) for line in lines))
    notes = f"""The R packages of COAL {args.coal_tag} and orthologr {args.orthologr_tag} in one place: `{coalr_file}` \
from https://github.com/{args.coal_repo}/releases/tag/{args.coal_tag} and `{orthologr_file}` from \
https://github.com/{args.orthologr_repo}/releases/tag/{args.orthologr_tag}, unchanged. Each was checked against its \
release's `SHA256SUMS.txt` and GitHub's record of the file before it was copied here.

Installation: see `INSTALL.md` (install coalr, then orthologr, then run `orthologr_sitrep()`).

| File | SHA-256 |
|---|---|
{table}

{AUTHORS_LINE}
"""
    with open(args.notes, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(notes)
    with open(args.notes, encoding="utf-8") as fh:
        if fh.read().rstrip("\n").splitlines()[-1] != AUTHORS_LINE:
            raise Failure("the release notes do not end with the authors line")
    return names


# ---------------------------------------------------------------------------------------------------------------

def write_github_output(**kv):
    path = os.environ.get("GITHUB_OUTPUT")
    if not path:
        return
    with open(path, "a", encoding="utf-8") as fh:
        for k, v in kv.items():
            fh.write(f"{k}={v}\n")


def main(argv):
    import argparse
    p = argparse.ArgumentParser(prog="release_tools.py", description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    pv = sub.add_parser("version")
    pv.add_argument("tag")
    pv.add_argument("--dry-run", action="store_true")
    pv.add_argument("--github-output", action="store_true")
    pv.add_argument("--root", default=".")
    pi = sub.add_parser("check-inputs")
    pi.add_argument("coal_tag")
    pi.add_argument("orthologr_tag")
    pi.add_argument("bundle_tag")
    pr = sub.add_parser("verify")
    pr.add_argument("dir")
    pr.add_argument("--release-json", required=True)
    pr.add_argument("--tag", required=True)
    pb = sub.add_parser("bundle")
    for a in ("coal-dir", "coal-tag", "coal-repo", "orthologr-dir", "orthologr-tag", "orthologr-repo",
              "bundle-tag", "out", "notes"):
        pb.add_argument("--" + a, required=True)
    args = p.parse_args(argv)
    try:
        if args.cmd == "version":
            version, tag = check_version(args.tag, args.root, args.dry_run)
            print(version)
            if args.github_output:
                write_github_output(version=version, tag=tag)
        elif args.cmd == "check-inputs":
            check_inputs(args.coal_tag, args.orthologr_tag, args.bundle_tag)
            print(f"inputs: COAL {args.coal_tag}, orthologr {args.orthologr_tag}, bundle {args.bundle_tag}")
        elif args.cmd == "verify":
            sums = verify_release(args.dir, args.release_json, args.tag)
            print(f"{args.tag}: {len(sums)} files, each matches {SUMS} and GitHub's record")
        elif args.cmd == "bundle":
            check_inputs(args.coal_tag, args.orthologr_tag, args.bundle_tag)
            for name in bundle(args):
                print(name)
    except (Failure, OSError, tarfile.TarError, json.JSONDecodeError) as e:
        print(f"release_tools: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
