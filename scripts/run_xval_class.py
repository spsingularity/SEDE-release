#!/usr/bin/env python3
"""
CLASS-vs-CAMB background cross-check: does SEDE's dark-energy w(a) give the SAME
distances and sound horizon r_drag in two INDEPENDENT Boltzmann codes?  Agreement to
<~0.1% confirms the headline geometry is not a CAMB code artifact.

SEDE's w(z) (canonical Barrow λ=0.5) is projected onto its CPL form (w0,wa) — the w(a)
parametrisation BOTH codes support identically (PPF) — and fed to CAMB and CLASS with a
matched cosmology (massless ν on both sides to avoid neutrino-convention mismatch). We
compare r_drag, D_M(z)/r_d, D_H(z)/r_d, H(z) over z∈[0.3,2.3] (the DESI range).

CLASS (`classy`) is optional. If absent, `compare()` returns skipped=True (the suite stays
green); install with `pip install classy` to activate the live cross-check.

Run:  python run_xval_class.py
"""
import numpy as np
C = 299792.458

def sede_cpl(Om=0.30, gamma=1.4964, zmax=3.0, n=60):
    """CPL (w0,wa) projection of the canonical SEDE-H w(z)."""
    from sede import friedmann as fr
    z = np.linspace(0, zmax, n); E = fr.E_SEDE_lambda(z, Om, gamma, 0.5)
    rho = np.maximum(E**2 - Om*(1+z)**3 - 9e-5*(1+z)**4, 1e-12)
    w = -1 + (1/3.)*(1+z)*np.gradient(np.log(rho), z)
    A = np.column_stack([np.ones_like(z), z/(1+z)]); c, *_ = np.linalg.lstsq(A, w, rcond=None)
    return float(c[0]), float(c[1])

def camb_background(Om, H0, ombh2, w0, wa, zg):
    import camb
    from camb import dark_energy
    h = H0/100.; pa = camb.CAMBparams()
    pa.set_cosmology(H0=H0, ombh2=ombh2, omch2=Om*h**2 - ombh2, mnu=0.0, num_massive_neutrinos=0)
    de = dark_energy.DarkEnergyPPF(); de.set_params(w=w0, wa=wa); pa.DarkEnergy = de
    pa.WantCls = False; pa.WantTransfer = False
    bg = camb.get_background(pa); d = bg.get_derived_params()
    DM = np.array([bg.comoving_radial_distance(z) for z in zg])
    DH = np.array([C/bg.hubble_parameter(z) for z in zg])
    Hz = np.array([bg.hubble_parameter(z) for z in zg])
    return d['rdrag'], DM, DH, Hz

def class_background(Om, H0, ombh2, w0, wa, zg):
    from classy import Class
    h = H0/100.
    cosmo = Class()
    cosmo.set({
        'h': h, 'omega_b': ombh2, 'omega_cdm': Om*h**2 - ombh2,
        'N_ur': 3.044, 'N_ncdm': 0,                       # massless neutrinos (match CAMB mnu=0)
        'Omega_Lambda': 0.0, 'w0_fld': w0, 'wa_fld': wa, 'use_ppf': 'yes',
        'output': '', 'background_verbose': 0,
    })
    cosmo.compute()
    rd = cosmo.get_current_derived_parameters(['rs_d'])['rs_d']
    DM = np.array([(1+z)*cosmo.angular_distance(z) for z in zg])     # comoving (flat)
    DH = np.array([1.0/cosmo.Hubble(z) for z in zg])                  # c/H, H in 1/Mpc -> Mpc
    Hz = np.array([cosmo.Hubble(z)*C for z in zg])                    # km/s/Mpc
    cosmo.struct_cleanup(); cosmo.empty()
    return rd, DM, DH, Hz

def compare(Om=0.30, H0=67.5, ombh2=0.02237, verbose=False):
    """Returns dict(skipped, ok, rd_camb, rd_class, max_frac_dev, ...)."""
    try:
        import classy  # noqa: F401
    except Exception:
        return {'skipped': True, 'ok': True, 'reason': 'classy not installed (pip install classy)'}
    w0, wa = sede_cpl(Om)
    zg = np.array([0.3, 0.5, 0.7, 0.93, 1.32, 1.49, 2.33])
    rdC, DMC, DHC, HC = camb_background(Om, H0, ombh2, w0, wa, zg)
    rdL, DML, DHL, HL = class_background(Om, H0, ombh2, w0, wa, zg)
    dev = {
        'r_drag': abs(rdC - rdL)/rdL,
        'DM/rd':  float(np.max(np.abs((DMC/rdC)/(DML/rdL) - 1))),
        'DH/rd':  float(np.max(np.abs((DHC/rdC)/(DHL/rdL) - 1))),
        'H(z)':   float(np.max(np.abs(HC/HL - 1))),
    }
    maxdev = max(dev.values()); ok = maxdev < 2e-3
    if verbose:
        print(f"  SEDE CPL: w0={w0:+.3f}, wa={wa:+.3f}")
        print(f"  r_drag: CAMB={rdC:.3f} Mpc  CLASS={rdL:.3f} Mpc  dev={dev['r_drag']:.2e}")
        for k in ('DM/rd', 'DH/rd', 'H(z)'):
            print(f"  {k:6s} max frac dev (CAMB vs CLASS): {dev[k]:.2e}")
    return {'skipped': False, 'ok': ok, 'rd_camb': rdC, 'rd_class': rdL,
            'max_frac_dev': maxdev, 'dev': dev, 'w0': w0, 'wa': wa}

if __name__ == "__main__":
    print("="*72); print("CLASS-vs-CAMB BACKGROUND CROSS-CHECK — SEDE w(a) in two Boltzmann codes"); print("="*72)
    r = compare(verbose=True)
    if r['skipped']:
        print(f"\n  SKIPPED: {r['reason']}")
        print("  -> install classy to activate; the test is wired into the suite (V33) and")
        print("     will verify CAMB and CLASS agree on r_drag and the DESI distances to <0.2%.")
    else:
        print(f"\n  max fractional deviation (CAMB vs CLASS) = {r['max_frac_dev']:.2e}  "
              f"-> {'AGREE (<0.2%, not a code artifact)' if r['ok'] else 'DISAGREE'}")
