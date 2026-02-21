"""
ADflow Turbulence Coefficient Interface via ctypes.

Bypasses f2py module variable access (which has a known memory duplication bug)
by directly writing to the Fortran paramTurb module symbols through ctypes.

Usage:
    from adflow_turb_ctypes import set_sa_constants, set_sst_constants, get_sa_constants, get_sst_constants

    # Set SA coefficients
    set_sa_constants(cb1=0.5, cb2=0.622, sigma=2/3, kappa=0.41, cv1=7.1, cw2=0.3, cw3=2.0, ct3=1.2, ct4=0.5)

    # Read back
    print(get_sa_constants())

    # Set SST coefficients
    set_sst_constants(sstk=0.41, a1=0.31, betas=0.09, sigk1=0.85, sigw1=0.5, beta1=0.075, sigk2=1.0, sigw2=0.856, beta2=0.0828)

Verified on Paracloud HPC (GCC 12.2 + OpenMPI 4.1.5), Job 36920788.
"""
import ctypes

# --- Module-level state ---
_dll = None
_sa_vars = {}
_sst_vars = {}

# gfortran symbol mangling: __<module>_MOD_<variable> (all lowercase)
_SA_SYMBOL_MAP = {
    'rsacb1': '__paramturb_MOD_rsacb1',
    'rsacb2': '__paramturb_MOD_rsacb2',
    'rsacb3': '__paramturb_MOD_rsacb3',   # sigma
    'rsak':   '__paramturb_MOD_rsak',     # kappa
    'rsacv1': '__paramturb_MOD_rsacv1',
    'rsacw1': '__paramturb_MOD_rsacw1',   # derived = cb1/k^2 + (1+cb2)/sigma
    'rsacw2': '__paramturb_MOD_rsacw2',
    'rsacw3': '__paramturb_MOD_rsacw3',
    'rsact1': '__paramturb_MOD_rsact1',
    'rsact2': '__paramturb_MOD_rsact2',
    'rsact3': '__paramturb_MOD_rsact3',
    'rsact4': '__paramturb_MOD_rsact4',
    'rsacrot':'__paramturb_MOD_rsacrot',
}

_SST_SYMBOL_MAP = {
    'rsstk':     '__paramturb_MOD_rsstk',
    'rssta1':    '__paramturb_MOD_rssta1',
    'rsstbetas': '__paramturb_MOD_rsstbetas',
    'rsstsigk1': '__paramturb_MOD_rsstsigk1',
    'rsstsigw1': '__paramturb_MOD_rsstsigw1',
    'rsstbeta1': '__paramturb_MOD_rsstbeta1',
    'rsstsigk2': '__paramturb_MOD_rsstsigk2',
    'rsstsigw2': '__paramturb_MOD_rsstsigw2',
    'rsstbeta2': '__paramturb_MOD_rsstbeta2',
}


def _ensure_init():
    """Load the libadflow.so and cache ctypes references."""
    global _dll, _sa_vars, _sst_vars
    if _dll is not None:
        return

    import adflow.libadflow as _lib
    _dll = ctypes.CDLL(_lib.__file__, mode=ctypes.RTLD_GLOBAL)

    for name, symbol in _SA_SYMBOL_MAP.items():
        try:
            _sa_vars[name] = ctypes.c_double.in_dll(_dll, symbol)
        except ValueError:
            pass  # symbol not found (e.g. if not patched)

    for name, symbol in _SST_SYMBOL_MAP.items():
        try:
            _sst_vars[name] = ctypes.c_double.in_dll(_dll, symbol)
        except ValueError:
            pass


# ============================================================
# SA Interface
# ============================================================

def set_sa_constants(cb1, cb2, sigma, kappa, cv1, cw2, cw3, ct3, ct4):
    """Set SA turbulence model closure coefficients.

    Parameters (9 calibration params):
        cb1   : production coefficient (default 0.1355)
        cb2   : diffusion coefficient (default 0.622)
        sigma : diffusion ratio, stored as rsaCb3 (default 2/3)
        kappa : von Karman constant (default 0.41)
        cv1   : wall damping coefficient (default 7.1)
        cw2   : destruction coefficient (default 0.3)
        cw3   : destruction coefficient (default 2.0)
        ct3   : transition coefficient (default 1.2)
        ct4   : transition coefficient (default 0.5)

    Note: cw1 is automatically recomputed as cb1/kappa^2 + (1+cb2)/sigma.
    """
    _ensure_init()
    _sa_vars['rsacb1'].value = cb1
    _sa_vars['rsacb2'].value = cb2
    _sa_vars['rsacb3'].value = sigma
    _sa_vars['rsak'].value = kappa
    _sa_vars['rsacv1'].value = cv1
    _sa_vars['rsacw2'].value = cw2
    _sa_vars['rsacw3'].value = cw3
    _sa_vars['rsact3'].value = ct3
    _sa_vars['rsact4'].value = ct4
    # Derived constant
    _sa_vars['rsacw1'].value = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma


def set_sa_defaults():
    """Reset SA coefficients to standard defaults."""
    set_sa_constants(
        cb1=0.1355, cb2=0.622, sigma=2.0/3.0, kappa=0.41,
        cv1=7.1, cw2=0.3, cw3=2.0, ct3=1.2, ct4=0.5
    )


def get_sa_constants():
    """Read current SA coefficients.

    Returns:
        dict with keys: cb1, cb2, sigma, kappa, cv1, cw1 (derived), cw2, cw3, ct3, ct4
    """
    _ensure_init()
    return {
        'cb1':   _sa_vars['rsacb1'].value,
        'cb2':   _sa_vars['rsacb2'].value,
        'sigma': _sa_vars['rsacb3'].value,
        'kappa': _sa_vars['rsak'].value,
        'cv1':   _sa_vars['rsacv1'].value,
        'cw1':   _sa_vars['rsacw1'].value,
        'cw2':   _sa_vars['rsacw2'].value,
        'cw3':   _sa_vars['rsacw3'].value,
        'ct3':   _sa_vars['rsact3'].value,
        'ct4':   _sa_vars['rsact4'].value,
    }


# ============================================================
# SST Interface
# ============================================================

def set_sst_constants(sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2):
    """Set SST turbulence model closure coefficients.

    Parameters (9 calibration params):
        sstk  : von Karman constant kappa (default 0.41)
        a1    : limiter constant (default 0.31)
        betas : beta* (default 0.09)
        sigk1 : sigma_k1 (default 0.85)
        sigw1 : sigma_omega1 (default 0.5)
        beta1 : beta_1 (default 0.075)
        sigk2 : sigma_k2 (default 1.0)
        sigw2 : sigma_omega2 (default 0.856)
        beta2 : beta_2 (default 0.0828)
    """
    _ensure_init()
    _sst_vars['rsstk'].value = sstk
    _sst_vars['rssta1'].value = a1
    _sst_vars['rsstbetas'].value = betas
    _sst_vars['rsstsigk1'].value = sigk1
    _sst_vars['rsstsigw1'].value = sigw1
    _sst_vars['rsstbeta1'].value = beta1
    _sst_vars['rsstsigk2'].value = sigk2
    _sst_vars['rsstsigw2'].value = sigw2
    _sst_vars['rsstbeta2'].value = beta2


def set_sst_defaults():
    """Reset SST coefficients to standard defaults."""
    set_sst_constants(
        sstk=0.41, a1=0.31, betas=0.09,
        sigk1=0.85, sigw1=0.5, beta1=0.075,
        sigk2=1.0, sigw2=0.856, beta2=0.0828
    )


def get_sst_constants():
    """Read current SST coefficients.

    Returns:
        dict with keys: sstk, a1, betas, sigk1, sigw1, beta1, sigk2, sigw2, beta2
    """
    _ensure_init()
    return {
        'sstk':  _sst_vars['rsstk'].value,
        'a1':    _sst_vars['rssta1'].value,
        'betas': _sst_vars['rsstbetas'].value,
        'sigk1': _sst_vars['rsstsigk1'].value,
        'sigw1': _sst_vars['rsstsigw1'].value,
        'beta1': _sst_vars['rsstbeta1'].value,
        'sigk2': _sst_vars['rsstsigk2'].value,
        'sigw2': _sst_vars['rsstsigw2'].value,
        'beta2': _sst_vars['rsstbeta2'].value,
    }
