#!/usr/bin/env python3
"""
Check that all Rust source files start with the OpenModelica Run-Time System
license header and that its copyright year is the current year.

Usage:
  check_license.py [--fix] [FILE...]

  --fix   Add missing headers and update outdated copyright years.
  FILE    Rust files to check. Default: all Rust files in the repository.

Exit codes:
  0  All files pass (or all failures were fixed with --fix).
  1  One or more files fail.

Based on .CI/scripts/check_runtime_license.py from the OpenModelica repository.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

CURRENT_YEAR = datetime.now().year

HEADER = """\
/*
 * This file belongs to the OpenModelica Run-Time System
 *
 * Copyright (c) 1998-{year}, Open Source Modelica Consortium (OSMC), c/o Linköpings
 * universitet, Department of Computer and Information Science, SE-58183 Linköping, Sweden. All rights
 * reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * AGPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.8. ANY
 * USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES RECIPIENT'S
 * ACCEPTANCE OF THE BSD NEW LICENSE OR THE OSMC PUBLIC LICENSE OR THE AGPL
 * VERSION 3, ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium) Public License
 * (OSMC-PL) are obtained from OSMC, either from the above address, from the URLs:
 * http://www.openmodelica.org or https://github.com/OpenModelica/ or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica distribution. GNU
 * AGPL version 3 is obtained from: https://www.gnu.org/licenses/licenses.html#GPL. The BSD NEW
 * License is obtained from: http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY
 * SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF
 * OSMC-PL.
 *
 */
"""

# The header with any year, to tell an outdated header from a missing one.
HEADER_ANY_YEAR = re.compile(
    re.escape(HEADER.format(year="@YEAR@")).replace("@YEAR@", r"(\d{4})")
)


def rust_files(root: Path) -> list[Path]:
    """Return all Rust files tracked by Git, plus untracked ones not ignored."""
    out = subprocess.run(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "*.rs"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return sorted(root / line for line in out.splitlines() if (root / line).exists())


def check_file(path: Path, fix: bool) -> str | None:
    """Return an error message for *path*, or None if it passes or was fixed."""
    content = path.read_text(encoding="utf-8")
    expected = HEADER.format(year=CURRENT_YEAR)
    if content.startswith(expected):
        return None

    m = HEADER_ANY_YEAR.match(content)
    if m:
        error = f"copyright year out of date ({m.group(1)}, expected {CURRENT_YEAR})"
        new_content = expected + content[m.end():]
    else:
        error = "license header missing or modified"
        body = content.lstrip("\n")
        new_content = expected + ("\n" + body if body else "")

    if fix:
        path.write_text(new_content, encoding="utf-8")
        return None
    return error


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--fix",
        action="store_true",
        help="add missing headers and update outdated copyright years",
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="Rust files to check (default: all Rust files in the repository)",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    paths = [f.resolve() for f in args.files] if args.files else rust_files(root)
    failed = 0
    for path in paths:
        error = check_file(path, args.fix)
        if error:
            failed += 1
            print(f"{path.relative_to(root)}: {error}", file=sys.stderr)

    if failed:
        print(
            f"\n{failed} file(s) failed. Run '.ci/check_license.py --fix' to fix them.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
