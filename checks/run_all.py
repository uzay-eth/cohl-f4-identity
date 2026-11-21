#!/usr/bin/env python3
"""Run every check of Cohl's Appell F4 identity and print PASS or FAIL.

Run from the repository root after installing requirements.txt:

    python3 checks/run_all.py

The four programs are independent re-implementations and cross-checks;
each one re-derives or re-sums both sides of the statement in its own way:

  1. checks/check_numeric.py  -- the proof's numeric check: both defining
     series summed literally as printed, per case at 85 and 115 decimal
     digits with two independent degree cutoffs (September 2023 statement
     transcription, seeded parameter cases plus deliberate special cases).
  2. checks/check_symbolic.py  -- the proof's exact check: every coefficient
     [t^k y^s] compared as an exact polynomial over QQ in the parameters.
  3. checks/check_opus.py     -- first independent review's check: exact
     rational coefficients, a fully symbolic comparison, a 60-digit
     numeric comparison including the Jacobi generating function and the
     Horn H3 form, and negative controls that must fail.
  4. checks/verify_glm.py      -- second independent review's verification:
     exact two-variable coefficient comparisons, proof-step identities,
     60-digit numeric summations, the Horn H3 reduction, and a mutation
     test of its own sensitivity.
  5. checks/verify_gemini_exact.py  -- third independent review's exact
     check: literal right side against the literal F4 series in exact
     rational arithmetic, polynomial-in-x coefficients, the Horn H3
     form, and four planted errors that must be caught.
  6. checks/verify_gemini_numeric.py  -- third independent review's
     numeric check: 110-digit decimal-complex summation of both literal
     series, the removable S-pole case, boundary and conditional-
     convergence diagnostics, and planted errors.

A check passes only if its program exits 0, prints a PASS summary and never
prints FAIL. The script exits non-zero if any check fails. Total run time
is about three minutes on a single core.

The first four programs default to a reduced size so that the suite
finishes quickly; the saved logs numeric.txt, symbolic.txt,
opus-output.txt and glm-output.txt record their full-size runs, which
those programs regenerate with --full (for check_symbolic.py, --order 40).
The two Gemini programs already run at full size; their logs
gemini-exact-output.txt and gemini-numeric-output.txt are always identical.
"""

import subprocess
import sys
import time

CHECKS = [
    ("check_numeric (proof, mpmath)", "checks/check_numeric.py"),
    ("check_symbolic (proof, SymPy over QQ)", "checks/check_symbolic.py"),
    ("check_opus (independent review 1)", "checks/check_opus.py"),
    ("verify_glm (independent review 2)", "checks/verify_glm.py"),
    ("verify_gemini_exact (independent review 3, exact)", "checks/verify_gemini_exact.py"),
    ("verify_gemini_numeric (independent review 3, numeric)", "checks/verify_gemini_numeric.py"),
]


def run_one(name, script):
    print(f"====== {name}: {script}")
    sys.stdout.flush()
    start = time.monotonic()
    proc = subprocess.Popen(
        [sys.executable, "-u", script],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    output = []
    for line in proc.stdout:
        print(line, end="")
        output.append(line)
    proc.wait()
    elapsed = time.monotonic() - start
    text = "".join(output)
    passed = proc.returncode == 0 and "PASS" in text and "FAIL" not in text
    print(f"====== {name}: {'PASS' if passed else 'FAIL'} "
          f"({elapsed:.0f} s, exit code {proc.returncode})\n")
    sys.stdout.flush()
    return passed


def main():
    print("Running all checks from the repository root. "
          f"python {sys.version.split()[0]}\n")
    results = []
    for name, script in CHECKS:
        results.append((name, run_one(name, script)))
    print("Summary:")
    for name, passed in results:
        print(f"  {name}: {'PASS' if passed else 'FAIL'}")
    failed = [name for name, passed in results if not passed]
    if failed:
        print(f"\n{len(failed)} check(s) FAILED.")
        return 1
    print("\nAll checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
