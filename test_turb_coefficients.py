"""
ADflow Turbulence Coefficient Validation (SA + SST) via ctypes
===============================================================
Verify that SA and SST closure coefficients can be read, written, and
reset at runtime through the ctypes interface (adflow_turb_ctypes.py).

Usage (after patch + rebuild):
    python3 test_turb_coefficients.py

NOTE: This script validates the ctypes API, NOT the f2py interface.
      f2py module variable access has a known memory duplication bug
      (reads back correctly but has zero effect on computation).
"""

import sys
import os

# Ensure adflow_turb_ctypes can be imported from the repo root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from adflow_turb_ctypes import (
    set_sa_constants, set_sa_defaults, get_sa_constants,
    set_sst_constants, set_sst_defaults, get_sst_constants,
)


def test_sa_defaults():
    """Test 1: Verify SA default values after set_sa_defaults()."""
    print("\n[Test 1] SA defaults after set_sa_defaults()")
    print("-" * 60)

    set_sa_defaults()
    c = get_sa_constants()

    expected = {
        "cb1":   0.1355,
        "cb2":   0.622,
        "sigma": 0.66666666667,   # 2/3
        "kappa": 0.41,
        "cv1":   7.1,
        "cw2":   0.3,
        "cw3":   2.0,
        "ct3":   1.2,
        "ct4":   0.5,
    }

    all_ok = True
    for key, exp in expected.items():
        actual = c[key]
        ok = abs(actual - exp) < 1e-6
        if not ok:
            all_ok = False
        print(f"  {key:8s} = {actual:12.8f}  expected {exp:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Check derived cw1
    expected_cw1 = 0.1355 / (0.41 ** 2) + (1.0 + 0.622) / (2.0 / 3.0)
    actual_cw1 = c["cw1"]
    ok = abs(actual_cw1 - expected_cw1) < 1e-6
    if not ok:
        all_ok = False
    print(f"  {'cw1':8s} = {actual_cw1:12.8f}  expected {expected_cw1:12.8f}  [{'OK' if ok else 'FAIL'}]  (derived)")

    print(f"\n  SA defaults: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_defaults():
    """Test 2: Verify SST default values after set_sst_defaults()."""
    print("\n[Test 2] SST defaults after set_sst_defaults()")
    print("-" * 60)

    set_sst_defaults()
    c = get_sst_constants()

    expected = {
        "sstk":  0.41,
        "a1":    0.31,
        "betas": 0.09,
        "sigk1": 0.85,
        "sigw1": 0.5,
        "beta1": 0.075,
        "sigk2": 1.0,
        "sigw2": 0.856,
        "beta2": 0.0828,
    }

    all_ok = True
    for key, exp in expected.items():
        actual = c[key]
        ok = abs(actual - exp) < 1e-6
        if not ok:
            all_ok = False
        print(f"  {key:8s} = {actual:12.8f}  expected {exp:12.8f}  [{'OK' if ok else 'FAIL'}]")

    print(f"\n  SST defaults: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sa_setter():
    """Test 3: set_sa_constants() with 9 params + auto cw1 recompute."""
    print("\n[Test 3] set_sa_constants(9 params) + cw1 auto-recompute")
    print("-" * 60)

    cb1, cb2, sigma, kappa = 0.14, 0.65, 0.7, 0.40
    cv1, cw2, cw3 = 7.0, 0.055, 2.5
    ct3, ct4 = 1.0, 0.4

    set_sa_constants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)
    c = get_sa_constants()

    checks = {
        "cb1":   cb1,
        "cb2":   cb2,
        "sigma": sigma,
        "kappa": kappa,
        "cv1":   cv1,
        "cw2":   cw2,
        "cw3":   cw3,
        "ct3":   ct3,
        "ct4":   ct4,
    }

    all_ok = True
    for key, exp in checks.items():
        actual = c[key]
        ok = abs(actual - exp) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {key:8s} = {actual:12.8f}  expected {exp:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Verify cw1 derived
    expected_cw1 = cb1 / (kappa ** 2) + (1.0 + cb2) / sigma
    actual_cw1 = c["cw1"]
    ok = abs(actual_cw1 - expected_cw1) < 1e-6
    if not ok:
        all_ok = False
    print(f"  {'cw1':8s} = {actual_cw1:12.8f}  expected {expected_cw1:12.8f}  [{'OK' if ok else 'FAIL'}]  (derived)")

    # Restore
    set_sa_defaults()

    print(f"\n  SA setter: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_setter():
    """Test 4: set_sst_constants() with 9 params."""
    print("\n[Test 4] set_sst_constants(9 params)")
    print("-" * 60)

    # Perturbed values
    params = {
        "sstk": 0.38, "a1": 0.28, "betas": 0.08,
        "sigk1": 0.75, "sigw1": 0.45, "beta1": 0.06,
        "sigk2": 0.9, "sigw2": 0.8, "beta2": 0.07,
    }

    set_sst_constants(**params)
    c = get_sst_constants()

    all_ok = True
    for key, exp in params.items():
        actual = c[key]
        ok = abs(actual - exp) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {key:8s} = {actual:12.8f}  expected {exp:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Restore
    set_sst_defaults()

    print(f"\n  SST setter: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sa_reset():
    """Test 5: Verify set_sa_defaults() restores after modification."""
    print("\n[Test 5] SA reset after modification")
    print("-" * 60)

    # Modify
    set_sa_constants(0.5, 0.8, 0.9, 0.35, 6.0, 0.4, 3.0, 1.5, 0.8)
    c = get_sa_constants()
    modified_cb1 = c["cb1"]

    # Reset
    set_sa_defaults()
    c = get_sa_constants()
    reset_cb1 = c["cb1"]

    ok_modified = abs(modified_cb1 - 0.5) < 1e-10
    ok_reset = abs(reset_cb1 - 0.1355) < 1e-6

    print(f"  After set:   cb1 = {modified_cb1:.8f}  (expect 0.5)       [{'OK' if ok_modified else 'FAIL'}]")
    print(f"  After reset: cb1 = {reset_cb1:.8f}  (expect 0.1355)    [{'OK' if ok_reset else 'FAIL'}]")

    all_ok = ok_modified and ok_reset
    print(f"\n  SA reset: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_reset():
    """Test 6: Verify set_sst_defaults() restores after modification."""
    print("\n[Test 6] SST reset after modification")
    print("-" * 60)

    # Modify
    set_sst_constants(0.35, 0.25, 0.07, 0.7, 0.4, 0.06, 0.8, 0.7, 0.06)
    c = get_sst_constants()
    modified_betas = c["betas"]

    # Reset
    set_sst_defaults()
    c = get_sst_constants()
    reset_betas = c["betas"]

    ok_modified = abs(modified_betas - 0.07) < 1e-10
    ok_reset = abs(reset_betas - 0.09) < 1e-6

    print(f"  After set:   betas = {modified_betas:.8f}  (expect 0.07)    [{'OK' if ok_modified else 'FAIL'}]")
    print(f"  After reset: betas = {reset_betas:.8f}  (expect 0.09)    [{'OK' if ok_reset else 'FAIL'}]")

    all_ok = ok_modified and ok_reset
    print(f"\n  SST reset: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def main():
    print("=" * 70)
    print("ADflow SA + SST Turbulence Coefficient Validation (ctypes API)")
    print("=" * 70)

    results = {}
    results["sa_defaults"] = test_sa_defaults()
    results["sst_defaults"] = test_sst_defaults()
    results["sa_setter"] = test_sa_setter()
    results["sst_setter"] = test_sst_setter()
    results["sa_reset"] = test_sa_reset()
    results["sst_reset"] = test_sst_reset()

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    all_pass = True
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"  {name:20s}: {status}")

    print()
    if all_pass:
        print("  ALL TESTS PASSED (ctypes API)")
    else:
        print("  SOME TESTS FAILED")
        sys.exit(1)

    print("=" * 70)


if __name__ == "__main__":
    main()
