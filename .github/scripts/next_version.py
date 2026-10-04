"""Compute the next calendar version (YYYY.MAJOR.FIX) from git release tags.

Release tags look like ``v2026.1.0``.  The first release of a calendar year is
``YYYY.1.0`` whatever the release type; after that a ``major`` release bumps
MAJOR and resets FIX, and a ``bugfix`` release bumps FIX.

Usage::

    python .github/scripts/next_version.py {major,bugfix}

Prints the new version (without the leading ``v``).  Exits non-zero if HEAD
already carries a release tag, so the same commit is not released twice.
"""
import argparse
import datetime
import re
import subprocess
import sys

TAG_PATTERN = re.compile(r"^v(\d{4})\.(\d+)\.(\d+)$")


def parse_tag(tag):
    match = TAG_PATTERN.match(tag.strip())
    if match is None:
        return None
    return tuple(int(part) for part in match.groups())


def next_version(tags, release_type, year):
    releases = [v for v in (parse_tag(t) for t in tags) if v is not None]
    latest = max(releases) if releases else None
    if latest is None or latest[0] < year:
        return year, 1, 0
    if latest[0] > year:
        raise ValueError("latest release v{}.{}.{} is dated after {}".format(*latest, year))
    _, major, fix = latest
    if release_type == "major":
        return year, major + 1, 0
    if release_type == "bugfix":
        return year, major, fix + 1
    raise ValueError("release type must be 'major' or 'bugfix', not {!r}".format(release_type))


def git_lines(*args):
    output = subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout
    return [line for line in output.splitlines() if line.strip()]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("release_type", choices=["major", "bugfix"])
    args = parser.parse_args(argv)

    head_releases = [t for t in git_lines("tag", "--points-at", "HEAD") if parse_tag(t)]
    if head_releases:
        sys.exit("HEAD is already released as {}".format(", ".join(head_releases)))

    year = datetime.datetime.now(datetime.timezone.utc).year
    print("{}.{}.{}".format(*next_version(git_lines("tag", "--list"), args.release_type, year)))


if __name__ == "__main__":
    main()
