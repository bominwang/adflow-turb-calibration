"""
ADflow Turbulence Coefficient Validation (SA + SST)
====================================================
Verify that SA and SST closure coefficients can be modified at runtime
after applying the paramTurb patch.

Usage (inside MDO Lab Docker, after patch + rebuild):
    python3 test_turb_coefficients.py
"""

import sys


def test_sa_defaults(pt):
    """Test 1: Verify SA default values after setSADefaults()."""
    print("\n[Test 1] SA defaults after setSADefaults()")
    print("-" * 60)

    pt.setsadefaults()

    sa_expected = {
        "rsak":   0.41,
        "rsacb1": 0.1355,
        "rsacb2": 0.622,
        "rsacb3": 0.66666666667,   # sigma = 2/3
        "rsacv1": 7.1,
        "rsacw2": 0.3,
        "rsacw3": 2.0,
        "rsact1": 1.0,
        "rsact2": 2.0,
        "rsact3": 1.2,
        "rsact4": 0.5,
        "rsacrot": 2.0,
    }

    all_ok = True
    for var, expected in sa_expected.items():
        actual = getattr(pt, var)
        ok = abs(actual - expected) < 1e-6
        if not ok:
            all_ok = False
        print(f"  {var:12s} = {actual:12.8f}  expected {expected:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Check derived cw1
    expected_cw1 = 0.1355 / (0.41 ** 2) + (1.0 + 0.622) / (2.0 / 3.0)
    actual_cw1 = pt.rsacw1
    ok = abs(actual_cw1 - expected_cw1) < 1e-6
    if not ok:
        all_ok = False
    print(f"  {'rsacw1':12s} = {actual_cw1:12.8f}  expected {expected_cw1:12.8f}  [{'OK' if ok else 'FAIL'}]  (derived)")

    print(f"\n  SA defaults: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_defaults(pt):
    """Test 2: Verify SST default values after setSSTDefaults()."""
    print("\n[Test 2] SST defaults after setSSTDefaults()")
    print("-" * 60)

    pt.setsstdefaults()

    sst_expected = {
        "rsstk":     0.41,
        "rssta1":    0.31,
        "rsstbetas": 0.09,
        "rsstsigk1": 0.85,
        "rsstsigw1": 0.5,
        "rsstbeta1": 0.075,
        "rsstsigk2": 1.0,
        "rsstsigw2": 0.856,
        "rsstbeta2": 0.0828,
    }

    all_ok = True
    for var, expected in sst_expected.items():
        actual = getattr(pt, var)
        ok = abs(actual - expected) < 1e-6
        if not ok:
            all_ok = False
        print(f"  {var:12s} = {actual:12.8f}  expected {expected:12.8f}  [{'OK' if ok else 'FAIL'}]")

    print(f"\n  SST defaults: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sa_readwrite(pt):
    """Test 3: SA individual variable read/write."""
    print("\n[Test 3] SA individual variable read/write")
    print("-" * 60)

    pt.setsadefaults()

    sa_vars = [
        "rsak", "rsacb1", "rsacb2", "rsacb3", "rsacv1",
        "rsacw1", "rsacw2", "rsacw3",
        "rsact1", "rsact2", "rsact3", "rsact4", "rsacrot",
    ]

    all_ok = True
    for var in sa_vars:
        original = getattr(pt, var)
        test_val = original * 1.5 + 0.01  # some different value
        setattr(pt, var, test_val)
        readback = getattr(pt, var)
        ok = abs(readback - test_val) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {var:12s}: set {test_val:.8f} -> read {readback:.8f}  [{'OK' if ok else 'FAIL'}]")
        # Restore
        setattr(pt, var, original)

    print(f"\n  SA read/write: {'ALL OK (13/13)' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_readwrite(pt):
    """Test 4: SST individual variable read/write."""
    print("\n[Test 4] SST individual variable read/write")
    print("-" * 60)

    pt.setsstdefaults()

    sst_vars = [
        "rsstk", "rssta1", "rsstbetas",
        "rsstsigk1", "rsstsigw1", "rsstbeta1",
        "rsstsigk2", "rsstsigw2", "rsstbeta2",
    ]

    all_ok = True
    for var in sst_vars:
        original = getattr(pt, var)
        test_val = original * 0.8 + 0.005
        setattr(pt, var, test_val)
        readback = getattr(pt, var)
        ok = abs(readback - test_val) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {var:12s}: set {test_val:.8f} -> read {readback:.8f}  [{'OK' if ok else 'FAIL'}]")
        setattr(pt, var, original)

    print(f"\n  SST read/write: {'ALL OK (9/9)' if all_ok else 'FAILED'}")
    return all_ok


def test_sa_setter(pt):
    """Test 5: setSAConstants() with 9 params + auto cw1 recompute."""
    print("\n[Test 5] setSAConstants(9 params) + cw1 auto-recompute")
    print("-" * 60)

    # Set non-default values
    cb1, cb2, sigma, kappa = 0.14, 0.65, 0.7, 0.40
    cv1, cw2, cw3 = 7.0, 0.055, 2.5
    ct3, ct4 = 1.0, 0.4

    pt.setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)

    checks = {
        "rsacb1": cb1,
        "rsacb2": cb2,
        "rsacb3": sigma,
        "rsak":   kappa,
        "rsacv1": cv1,
        "rsacw2": cw2,
        "rsacw3": cw3,
        "rsact3": ct3,
        "rsact4": ct4,
    }

    all_ok = True
    for var, expected in checks.items():
        actual = getattr(pt, var)
        ok = abs(actual - expected) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {var:12s} = {actual:12.8f}  expected {expected:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Verify cw1 derived
    expected_cw1 = cb1 / (kappa ** 2) + (1.0 + cb2) / sigma
    actual_cw1 = pt.rsacw1
    ok = abs(actual_cw1 - expected_cw1) < 1e-6
    if not ok:
        all_ok = False
    print(f"  {'rsacw1':12s} = {actual_cw1:12.8f}  expected {expected_cw1:12.8f}  [{'OK' if ok else 'FAIL'}]  (derived)")

    # Verify ct1, ct2, crot unchanged (not in setter)
    pt.setsadefaults()
    pt.setsaconstants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4)
    for var, expected in [("rsact1", 1.0), ("rsact2", 2.0), ("rsacrot", 2.0)]:
        actual = getattr(pt, var)
        ok = abs(actual - expected) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {var:12s} = {actual:12.8f}  expected {expected:12.8f}  [{'OK' if ok else 'FAIL'}]  (unchanged by setter)")

    # Restore
    pt.setsadefaults()

    print(f"\n  SA setter: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_sst_setter(pt):
    """Test 6: setSSTConstants() with 9 params."""
    print("\n[Test 6] setSSTConstants(9 params)")
    print("-" * 60)

    # Set non-default values (ONERA M6 midpoints)
    sstk  = 0.41
    a1    = 0.31
    betas = 0.09
    sigk1 = 0.85
    sigw1 = 0.5
    beta1 = 0.075
    sigk2 = 1.0
    sigw2 = 0.856
    beta2 = 0.0828

    # Use slightly perturbed values
    sstk_p  = 0.38
    a1_p    = 0.28
    betas_p = 0.08
    sigk1_p = 0.75
    sigw1_p = 0.45
    beta1_p = 0.06
    sigk2_p = 0.9
    sigw2_p = 0.8
    beta2_p = 0.07

    pt.setsstconstants(sstk_p, a1_p, betas_p,
                       sigk1_p, sigw1_p, beta1_p,
                       sigk2_p, sigw2_p, beta2_p)

    checks = {
        "rsstk":     sstk_p,
        "rssta1":    a1_p,
        "rsstbetas": betas_p,
        "rsstsigk1": sigk1_p,
        "rsstsigw1": sigw1_p,
        "rsstbeta1": beta1_p,
        "rsstsigk2": sigk2_p,
        "rsstsigw2": sigw2_p,
        "rsstbeta2": beta2_p,
    }

    all_ok = True
    for var, expected in checks.items():
        actual = getattr(pt, var)
        ok = abs(actual - expected) < 1e-10
        if not ok:
            all_ok = False
        print(f"  {var:12s} = {actual:12.8f}  expected {expected:12.8f}  [{'OK' if ok else 'FAIL'}]")

    # Restore
    pt.setsstdefaults()

    print(f"\n  SST setter: {'ALL OK' if all_ok else 'FAILED'}")
    return all_ok


def test_discover_all(pt):
    """Test 0: Discover all paramturb variables and subroutines."""
    print("\n[Test 0] Discover paramturb module contents")
    print("-" * 60)

    all_items = sorted([v for v in dir(pt) if not v.startswith("_")])
    vars_list = []
    subs_list = []
    for item in all_items:
        if callable(getattr(pt, item)):
            subs_list.append(item)
        else:
            vars_list.append(item)

    print(f"  Variables ({len(vars_list)}):")
    for v in vars_list:
        val = getattr(pt, v)
        print(f"    {v:15s} = {val}")
    print(f"  Subroutines ({len(subs_list)}):")
    for s in subs_list:
        print(f"    {s}")

    # Expect 22 variables (13 SA + 9 SST) and 4 subroutines
    ok_vars = len(vars_list) >= 22
    ok_subs = len(subs_list) >= 4
    print(f"\n  Variables >= 22: {'OK' if ok_vars else f'FAIL ({len(vars_list)})'}")
    print(f"  Subroutines >= 4: {'OK' if ok_subs else f'FAIL ({len(subs_list)})'}")
    return ok_vars and ok_subs


def main():
    print("=" * 70)
    print("ADflow SA + SST Turbulence Coefficient Validation")
    print("=" * 70)

    from adflow import libadflow
    pt = libadflow.paramturb
    print("  libadflow.paramturb loaded OK")

    results = {}
    results["discover"] = test_discover_all(pt)
    results["sa_defaults"] = test_sa_defaults(pt)
    results["sst_defaults"] = test_sst_defaults(pt)
    results["sa_readwrite"] = test_sa_readwrite(pt)
    results["sst_readwrite"] = test_sst_readwrite(pt)
    results["sa_setter"] = test_sa_setter(pt)
    results["sst_setter"] = test_sst_setter(pt)

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
        print("  ALL TESTS PASSED")
    else:
        print("  SOME TESTS FAILED")
        sys.exit(1)

    print("=" * 70)


if __name__ == "__main__":
    main()
