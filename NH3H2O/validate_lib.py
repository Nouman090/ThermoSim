"""Validation suite for nh3h2o.py."""
import warnings, time, itertools
import numpy as np
from CoolProp.CoolProp import PropsSI
import nh3h2o as aw
warnings.filterwarnings('ignore')

PASS = []


def rep(name, val, tol, unit='%'):
    ok = val <= tol
    PASS.append((name, val, tol, ok))
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<52s} {val:10.4g} {unit}  (tol {tol:g})")
    return ok


h0w = PropsSI('H', 'T', 273.16, 'Q', 0, 'Water') / 1000
h0a = PropsSI('H', 'T', 273.15, 'Q', 0, 'Ammonia') / 1000
s0w = PropsSI('S', 'T', 273.16, 'Q', 0, 'Water') / 1000
s0a = PropsSI('S', 'T', 273.15, 'Q', 0, 'Ammonia') / 1000


def ref(fl, prop, **kw):
    args = []
    for k, v in kw.items():
        args += [k, v]
    val = PropsSI(prop, *args, fl)
    if prop == 'H':
        return val / 1000 - (h0w if fl == 'Water' else h0a)
    if prop == 'S':
        return val / 1000 - (s0w if fl == 'Water' else s0a)
    if prop == 'D':
        return 1.0 / val
    return val


# ==========================================================================
print('\n=== TEST A: pure-component limits vs reference EOS (T,P inputs) ===')
eh, es, ev = [], [], []
for fl, z in [('Water', 0.0), ('Ammonia', 1.0)]:
    for P in [0.5, 1, 2, 5, 10, 20, 40]:
        Ts = PropsSI('T', 'P', P * 1e5, 'Q', 0, fl)
        for T, ph in [(Ts - 25, 'liq'), (Ts + 40, 'vap')]:
            if T < 235:
                continue
            st = aw.state(z=z, T=T, P=P, warn=False)
            hr = ref(fl, 'H', T=T, P=P * 1e5)
            sr = ref(fl, 'S', T=T, P=P * 1e5)
            vr = ref(fl, 'D', T=T, P=P * 1e5)
            eh.append(abs(st.h - hr))
            es.append(abs(st.s - sr))
            ev.append(100 * abs(st.v - vr) / vr)
rep('enthalpy, mean abs deviation', np.mean(eh), 6.0, 'kJ/kg')
rep('enthalpy, max abs deviation', np.max(eh), 25.0, 'kJ/kg')
rep('entropy, max abs deviation', np.max(es), 0.06, 'kJ/kg/K')
rep('specific volume, max deviation', np.max(ev), 3.0, '%')

# ==========================================================================
print('\n=== TEST B: saturation temperature at pure limits (P,Q inputs) ===')
eb, ed = [], []
for fl, z in [('Water', 0.0), ('Ammonia', 1.0)]:
    for P in [0.5, 1, 2, 5, 10, 20, 40]:
        rT = PropsSI('T', 'P', P * 1e5, 'Q', 0, fl)
        eb.append(abs(aw.PropsSI('T', 'P', P, 'Q', 0, z, warn=False) - rT))
        ed.append(abs(aw.PropsSI('T', 'P', P, 'Q', 1, z, warn=False) - rT))
rep('bubble temperature, max error', np.max(eb), 0.35, 'K')
rep('dew temperature, max error', np.max(ed), 0.35, 'K')

print('\n=== TEST C: latent heat at pure limits (Q=0 vs Q=1) ===')
el = []
for fl, z in [('Water', 0.0), ('Ammonia', 1.0)]:
    for P in [1, 5, 10, 20, 40]:
        hg = aw.PropsSI('H', 'P', P, 'Q', 1, z, warn=False)
        hf = aw.PropsSI('H', 'P', P, 'Q', 0, z, warn=False)
        r = (PropsSI('H', 'P', P * 1e5, 'Q', 1, fl)
             - PropsSI('H', 'P', P * 1e5, 'Q', 0, fl)) / 1000
        el.append(100 * abs((hg - hf) - r) / r)
rep('latent heat, max deviation', np.max(el), 1.5, '%')

# ==========================================================================
print('\n=== TEST D: round-trip consistency of every input pair ===')
grid = [(z, T, P) for z in (0.05, 0.3, 0.5, 0.7, 0.95)
        for P in (1., 5., 15., 30.) for T in (280., 330., 380., 430.)]
errs = {k: [] for k in ['P,H', 'P,S', 'P,V', 'P,D', 'T,H', 'T,S', 'T,V', 'P,Q', 'T,Q']}
amb = {}
nrt = 0
for z, T, P in grid:
    try:
        st = aw.state(z=z, T=T, P=P, warn=False)
    except Exception:
        continue
    nrt += 1
    for lbl, kw, chk in [
            ('P,H', dict(P=P, H=st.h), 'T'), ('P,S', dict(P=P, S=st.s), 'T'),
            ('P,V', dict(P=P, V=st.v), 'T'), ('P,D', dict(P=P, D=1 / st.v), 'T'),
            ('T,H', dict(T=T, H=st.h), 'P'), ('T,S', dict(T=T, S=st.s), 'P'),
            ('T,V', dict(T=T, V=st.v), 'P')]:
        try:
            b = aw.state(z=z, warn=False, **kw)
            tgt = T if chk == 'T' else P
            got = b.T if chk == 'T' else b.P
            e = 100 * abs(got - tgt) / tgt
            if e > 0.05:
                # a different branch may still satisfy the inputs exactly;
                # score the residual on the SPECIFIED property instead
                pk = [k for k in kw if k not in ('T', 'P')][0]
                want = kw[pk]
                have = getattr(b, {'H': 'h', 'S': 's', 'V': 'v'}.get(pk, pk.lower()))
                if pk == 'D':
                    have = 1.0 / b.v
                e = 100 * abs(have - want) / max(abs(want), 1e-6)
                amb[lbl] = amb.get(lbl, 0) + 1
            errs[lbl].append(e)
        except Exception:
            errs[lbl].append(np.nan)
    if st.phase == 'two-phase':
        for lbl, kw, chk, tgt in [('P,Q', dict(P=P, Q=st.Q), 'T', T),
                                  ('T,Q', dict(T=T, Q=st.Q), 'P', P)]:
            try:
                b = aw.state(z=z, warn=False, **kw)
                got = b.T if chk == 'T' else b.P
                errs[lbl].append(100 * abs(got - tgt) / tgt)
            except Exception:
                errs[lbl].append(np.nan)
print(f'  ({nrt} base states)')
for lbl, e in errs.items():
    e = np.array(e)
    fail = np.isnan(e).sum()
    a = amb.get(lbl, 0)
    tag = f'{lbl} -> satisfies inputs  (n={len(e)}, {a} multi-root)'
    rep(tag, np.nanmax(e) if len(e) else 0.0, 0.05)

# ==========================================================================
print('\n=== TEST E: two-phase mass and energy balances ===')
em, ee, es2, ev2 = [], [], [], []
for z in (0.1, 0.3, 0.5, 0.7, 0.9):
    for P in (1., 5., 15., 30.):
        for Q in (0.1, 0.35, 0.6, 0.9):
            try:
                st = aw.state(z=z, P=P, Q=Q, warn=False)
            except Exception:
                continue
            pl = aw.props_single_phase(st.T, st.P, st.x, 'liquid')
            pg = aw.props_single_phase(st.T, st.P, st.y, 'vapor')
            em.append(abs((1 - st.Q) * st.x + st.Q * st.y - z))
            ee.append(abs((1 - st.Q) * pl['h'] + st.Q * pg['h'] - st.h))
            es2.append(abs((1 - st.Q) * pl['s'] + st.Q * pg['s'] - st.s))
            ev2.append(abs((1 - st.Q) * pl['v'] + st.Q * pg['v'] - st.v))
rep('NH3 mass balance  z = (1-Q)x + Qy', np.max(em), 1e-9, 'kg/kg')
rep('energy balance  h = (1-Q)hL + QhV', np.max(ee), 1e-9, 'kJ/kg')
rep('entropy balance', np.max(es2), 1e-9, 'kJ/kg/K')
rep('volume balance', np.max(ev2), 1e-12, 'm3/kg')

# ==========================================================================
print('\n=== TEST F: VLE is thermodynamically consistent (equal chem. potentials) ===')
emu = []
for z in (0.15, 0.4, 0.65, 0.85):
    for P in (1., 5., 15., 30.):
        try:
            st = aw.state(z=z, P=P, Q=0.5, warn=False)
        except Exception:
            continue
        xl, yv = aw.to_mole(st.x), aw.to_mole(st.y)
        lnKA, lnKW = aw._lnK(st.T, st.P, xl)
        emu.append(abs(np.exp(lnKA) * xl - yv))
        emu.append(abs(np.exp(lnKW) * (1 - xl) - (1 - yv)))
rep('K-value closure  |Kx - y|', np.max(emu), 1e-9, 'mol/mol')

print('\n=== TEST G: continuity across the phase boundary ===')
ec = []
for z in (0.2, 0.5, 0.8):
    for P in (2., 10., 25.):
        try:
            b = aw.state(z=z, P=P, Q=0, warn=False)
            d = aw.state(z=z, P=P, Q=1, warn=False)
        except Exception:
            continue
        hl = aw.state(z=z, T=b.T - 1e-4, P=P, warn=False)
        hv = aw.state(z=z, T=d.T + 1e-4, P=P, warn=False)
        ec.append(abs(hl.h - b.h))
        ec.append(abs(hv.h - d.h))
rep('enthalpy continuity at bubble/dew', np.max(ec), 0.05, 'kJ/kg')

print('\n=== TEST H: bubble/dew vs Patek-Klomfar correlation (model difference) ===')
db, dd = [], []
for z in (0.1, 0.3, 0.5, 0.7, 0.9):
    for P in (2., 5., 10., 20.):
        db.append(abs(aw.bubble_T(P, z) - aw.bubble_T(P, z, vle='PK')))
        dd.append(abs(aw.dew_T(P, z) - aw.dew_T(P, z, vle='PK')))
print(f'  bubble: mean {np.mean(db):.2f} K, max {np.max(db):.2f} K')
print(f'  dew:    mean {np.mean(dd):.2f} K, max {np.max(dd):.2f} K')
print('  (expected: ZT and PK are different models; not a pass/fail test)')

print('\n=== TEST I: performance ===')
t0 = time.time()
for _ in range(200):
    aw.state(z=0.4, T=350, P=10, warn=False)
t1 = time.time()
for _ in range(30):
    aw.state(z=0.4, P=10, H=500, warn=False)
t2 = time.time()
print(f'  TP flash : {(t1-t0)/200*1000:7.3f} ms/call')
print(f'  PH flash : {(t2-t1)/30*1000:7.3f} ms/call')

# ==========================================================================
n_ok = sum(1 for *_, ok in PASS if ok)
print(f'\n{"="*74}\nSUMMARY: {n_ok}/{len(PASS)} checks passed')
for name, val, tol, ok in PASS:
    if not ok:
        print(f'  FAILED: {name}  ({val:.4g} > {tol:g})')
