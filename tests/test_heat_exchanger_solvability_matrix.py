"""
test_heat_exchanger_solvability_matrix.py
-----------------------------------------
Case-to-test matrix for the pinch-constrained heat-exchanger solver.

Every specification of the six exchanger quantities
    {H_hi, H_ho, H_ci, H_co, mh, mc}
with two or three of them unknown is exercised against a reference state that
was built FORWARD (all six quantities known and mutually consistent, with the
minimum approach equal to the specified pinch).  The unknowns are then hidden
and the solver must recover them.

    2 unknowns: C(6,2) = 15 specifications -> 12 solved, 3 refused
    3 unknowns: C(6,3) = 20 specifications ->  8 solved, 12 refused
                (for both HEX_type='Evaporator' and 'Condenser')

Acceptance criteria
    * solved   : every unknown recovered to a relative error below
                 REL_TOL_2 (two unknowns) or REL_TOL_3 (three unknowns);
                 minimum approach within DT_TOL of the specified pinch;
                 no warning raised.
    * refused  : a ValueError whose message names the reason, and the model
                 left EXACTLY as it was (every attribute of every state
                 point, including mass flows).

The module also covers terminal temperature crossings, an internal pinch
violation (warn by default, raise with strict_pinch=True), mid-solve failures
and rollback between anchors.

Run with:   pytest tests/test_heat_exchanger_solvability_matrix.py -v
"""
import sys
import os
import itertools
import warnings

import pytest
from CoolProp.CoolProp import PropsSI as PS
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ThermoSim import ThermodynamicModel, HeatExchanger

NAMES = ['H_hi', 'H_ho', 'H_ci', 'H_co', 'mh', 'mc']
PPT = 8.0
DIV_N = 20                 # the recovered unknowns are N-independent (see S9)
REL_TOL_2 = 5e-4           # observed maximum 1.1e-4 (brentq rtol = 1e-4)
REL_TOL_3 = 1e-4           # observed maximum 7.4e-6
DT_TOL = 5e-3              # K


# --------------------------------------------------------------------------
#  Reference states (built forward, then hidden)
# --------------------------------------------------------------------------
def _orc_evaporator_interior_pinch():
    """R245fa evaporated at 15 bar by 130 C brine (5 bar, 10 kg/s).  The cold
    stream leaves as saturated vapour, so the 8 K minimum approach sits at the
    boiling onset inside the exchanger (used for the 2-unknown specifications)."""
    pev, pb = 15e5, 5e5
    h3 = PS('H', 'T', 313.15, 'Q', 0, 'R245fa')
    s3 = PS('S', 'T', 313.15, 'Q', 0, 'R245fa')
    h_ci = h3 + (PS('H', 'P', pev, 'S', s3, 'R245fa') - h3) / 0.75
    h_co = PS('H', 'P', pev, 'Q', 1, 'R245fa')
    h_hi = PS('H', 'T', 403.15, 'P', pb, 'Water')
    mc, mh = 4.7347, 10.0
    h_ho = h_hi - mc * (h_co - h_ci) / mh
    return dict(truth=dict(H_hi=h_hi, H_ho=h_ho, H_ci=h_ci, H_co=h_co, mh=mh, mc=mc),
                fluids=('Water', pb, 'R245fa', pev), hex='Evaporator', heat_added=False)


def _orc_evaporator_hot_end_anchor():
    """Same exchanger with the cold stream leaving at T_hi - 8 K (122 C), so the
    hot-end anchor T_co = T_hi - PPT is consistent (used for 3 unknowns)."""
    pev, pb = 15e5, 5e5
    h3 = PS('H', 'T', 313.15, 'Q', 0, 'R245fa')
    s3 = PS('S', 'T', 313.15, 'Q', 0, 'R245fa')
    h_ci = h3 + (PS('H', 'P', pev, 'S', s3, 'R245fa') - h3) / 0.75
    h_hi = PS('H', 'T', 403.15, 'P', pb, 'Water')
    h_co = PS('H', 'T', 395.15, 'P', pev, 'R245fa')
    mc, mh = 4.1321, 10.0
    h_ho = h_hi - mc * (h_co - h_ci) / mh
    return dict(truth=dict(H_hi=h_hi, H_ho=h_ho, H_ci=h_ci, H_co=h_co, mh=mh, mc=mc),
                fluids=('Water', pb, 'R245fa', pev), hex='Evaporator', heat_added=False)


def _orc_condenser_cold_end_anchor():
    """R245fa condensed at 3 bar by water (2 bar).  Cold inlet is T_ho - 8 K and
    the water flow is chosen so the approach at the dew point is also 8 K."""
    ph, pc = 3e5, 2e5
    h_hi = PS('H', 'T', 343.15, 'P', ph, 'R245fa')
    h_ho = PS('H', 'T', 308.15, 'P', ph, 'R245fa')
    h_ci = PS('H', 'T', 300.15, 'P', pc, 'Water')
    mh = 2.0
    q = mh * (h_hi - h_ho)
    h_dew = PS('H', 'P', ph, 'Q', 1, 'R245fa')
    t_dew = PS('T', 'P', ph, 'Q', 1, 'R245fa')
    q_dew = mh * (h_dew - h_ho)
    mc = brentq(lambda m: t_dew - PS('T', 'H', h_ci + q_dew / m, 'P', pc, 'Water') - 8.0, 5, 200)
    return dict(truth=dict(H_hi=h_hi, H_ho=h_ho, H_ci=h_ci, H_co=h_ci + q / mc, mh=mh, mc=mc),
                fluids=('R245fa', ph, 'Water', pc), hex='Condenser', heat_added=None)


REF_2 = _orc_evaporator_interior_pinch()
REF_3 = {'Evaporator': _orc_evaporator_hot_end_anchor(),
         'Condenser': _orc_condenser_cold_end_anchor()}

#: specifications the solver must refuse, with the reason it must give
REFUSED_2 = {
    frozenset({'H_hi', 'H_co'}): 'Unsolvable 2-unknown',
    frozenset({'H_ho', 'H_ci'}): 'Unsolvable 2-unknown',
    frozenset({'mh', 'mc'}):     'Unsolvable 2-unknown',
}
SOLVED_3 = {
    frozenset(s) for s in (
        ('H_hi', 'H_ho', 'mh'), ('H_hi', 'H_ho', 'mc'),
        ('H_hi', 'H_ci', 'mh'), ('H_hi', 'H_ci', 'mc'),
        ('H_ho', 'H_co', 'mh'), ('H_ho', 'H_co', 'mc'),
        ('H_ci', 'H_co', 'mh'), ('H_ci', 'H_co', 'mc'),
    )
}


# --------------------------------------------------------------------------
#  Helpers
# --------------------------------------------------------------------------
def build(ref, unknown, ppt=PPT, strict=False, truth=None):
    truth = truth or ref['truth']
    fh, ph, fc, pc = ref['fluids']
    M = ThermodynamicModel()

    def pt(name, fluid, p, hname, mname):
        d = dict(P=p)
        if hname not in unknown:
            d['H'] = truth[hname]
        if mname not in unknown:
            d['Mass_flowrate'] = truth[mname]
        M.add_point(fluid, name, **d)

    pt('hi', fh, ph, 'H_hi', 'mh')
    pt('ho', fh, ph, 'H_ho', 'mh')
    pt('ci', fc, pc, 'H_ci', 'mc')
    pt('co', fc, pc, 'H_co', 'mc')
    hx = HeatExchanger(M, 'X', ppt, ref['hex'], ref['heat_added'],
                       'hi', 'ho', 'ci', 'co', div_N=DIV_N, strict_pinch=strict)
    M.Component['X'] = hx
    return M, hx


def recovered(M):
    return dict(H_hi=M.Point['hi'].H, H_ho=M.Point['ho'].H,
                H_ci=M.Point['ci'].H, H_co=M.Point['co'].H,
                mh=M.Point['hi'].Mass_flowrate, mc=M.Point['ci'].Mass_flowrate)


def contents(M):
    """Every attribute of every state point, for an exact before/after test."""
    return {k: dict(vars(p)) for k, p in M.Point.items()}


def solve_quietly(hx):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        hx.Cal()
    return w


# --------------------------------------------------------------------------
#  Two unknowns: 15 specifications
# --------------------------------------------------------------------------
@pytest.mark.parametrize('unknown', list(itertools.combinations(NAMES, 2)),
                         ids=lambda u: '{' + ','.join(u) + '}')
def test_two_unknown_matrix(unknown):
    unknown = set(unknown)
    M, hx = build(REF_2, unknown)
    if frozenset(unknown) in REFUSED_2:
        before = contents(M)
        with pytest.raises(ValueError, match=REFUSED_2[frozenset(unknown)]):
            hx.Cal()
        assert contents(M) == before, 'refusal must leave the model untouched'
        return
    w = solve_quietly(hx)
    got, truth = recovered(M), REF_2['truth']
    for n in unknown:
        assert got[n] == pytest.approx(truth[n], rel=REL_TOL_2), n
    assert hx.min_dT == pytest.approx(PPT, abs=DT_TOL)
    assert not w, [str(x.message) for x in w]


# --------------------------------------------------------------------------
#  Three unknowns: 20 specifications x {Evaporator, Condenser}
# --------------------------------------------------------------------------
@pytest.mark.parametrize('hex_type', ['Evaporator', 'Condenser'])
@pytest.mark.parametrize('unknown', list(itertools.combinations(NAMES, 3)),
                         ids=lambda u: '{' + ','.join(u) + '}')
def test_three_unknown_matrix(hex_type, unknown):
    ref, unknown = REF_3[hex_type], set(unknown)
    M, hx = build(ref, unknown)
    if frozenset(unknown) not in SOLVED_3:
        before = contents(M)
        with pytest.raises(ValueError, match='Unsolvable 3-unknown|No valid anchor'):
            hx.Cal()
        assert contents(M) == before, 'refusal must leave the model untouched'
        return
    w = solve_quietly(hx)
    got, truth = recovered(M), ref['truth']
    for n in unknown:
        assert got[n] == pytest.approx(truth[n], rel=REL_TOL_3), n
    assert hx.min_dT == pytest.approx(PPT, abs=DT_TOL)
    assert not w, [str(x.message) for x in w]


def test_matrix_sizes():
    """12 of 15 and 8 of 20 specifications are admissible."""
    assert len(list(itertools.combinations(NAMES, 2))) - len(REFUSED_2) == 12
    assert len(SOLVED_3) == 8 and len(list(itertools.combinations(NAMES, 3))) == 20


# --------------------------------------------------------------------------
#  Terminal temperature crossings (specification contradicts the pinch)
# --------------------------------------------------------------------------
@pytest.mark.parametrize('t_out_c, dT', [(135.0, -5.0), (128.0, 2.0)])
def test_terminal_crossing_is_refused_and_rolled_back(t_out_c, dT):
    """Cold outlet specified above (or within less than PPT of) the hot inlet."""
    ref = REF_2
    M = ThermodynamicModel()
    fh, ph, fc, pc = ref['fluids']
    M.add_point('Water', 'hi', Mass_flowrate=10.0, T=403.15, P=ph)
    M.add_point('Water', 'ho', P=ph)
    M.add_point('R245fa', 'ci', P=pc, H=ref['truth']['H_ci'])
    M.add_point('R245fa', 'co', P=pc, T=t_out_c + 273.15)
    hx = HeatExchanger(M, 'X', PPT, 'Evaporator', True, 'hi', 'ho', 'ci', 'co', div_N=DIV_N)
    M.Component['X'] = hx
    before = contents(M)
    with pytest.raises(ValueError, match='PPT violation at hot-inlet/cold-outlet') as e:
        hx.Cal()
    assert f'{dT:.3f}' in str(e.value)
    assert contents(M) == before, 'mass-flow propagation must be undone as well'


# --------------------------------------------------------------------------
#  Internal pinch violation: warn by default, raise with strict_pinch=True
# --------------------------------------------------------------------------
def _overloaded_evaporator(strict):
    """All six quantities supplied, consistent in energy, but with more cold
    flow than the 8 K pinch allows: the terminal approaches are fine, the
    boiling-onset approach is below 8 K."""
    truth = dict(REF_2['truth'])
    truth['mc'] = 5.3
    truth['H_ho'] = (truth['H_hi'] - truth['mc'] * (truth['H_co'] - truth['H_ci']) / truth['mh'])
    return build(REF_2, set(), strict=strict, truth=truth)


def test_internal_pinch_violation_warns_but_returns():
    M, hx = _overloaded_evaporator(strict=False)
    w = solve_quietly(hx)
    assert hx.min_dT < PPT - 0.5
    assert any('pinch' in str(x.message).lower() or 'PPT' in str(x.message) for x in w)


def test_internal_pinch_violation_raises_when_strict():
    M, hx = _overloaded_evaporator(strict=True)
    before = contents(M)
    with pytest.raises(ValueError, match='PPT|pinch'):
        hx.Cal()
    assert contents(M) == before


# --------------------------------------------------------------------------
#  Rollback
# --------------------------------------------------------------------------
def test_midsolve_refusal_leaves_model_untouched():
    """{H_hi,H_ho,mc} with an interior-pinch reference state is not consistent
    with the hot-end anchor; the solver gets as far as the energy balance and
    refuses ("heat flows the wrong way").  Up to v3.4.0 the cold-side mass flow
    and exergy of that attempt stayed in Model.Point."""
    M, hx = build(REF_2, {'H_hi', 'H_ho', 'mc'})
    before = contents(M)
    with pytest.raises(ValueError):
        hx.Cal()
    assert contents(M) == before


def test_state_restored_between_anchors():
    """With HEX_type='double_pipe' both the hot-end and the cold-end anchor are
    admissible for {H_ci,H_co,mh}.  If the first attempt fails part-way, the
    second must start from the state the exchanger had on entry."""
    ref = dict(REF_3['Evaporator'], hex='double_pipe', heat_added=None)
    M, hx = build(ref, {'H_ci', 'H_co', 'mh'})
    entry = contents(M)
    original = hx._solve_2_unknown
    seen = []

    def flaky(*args, **kwargs):
        seen.append(contents(M))
        if len(seen) == 1:                        # corrupt the model, then fail
            M.Point['hi'].H = 1.0
            M.Point['ci'].Mass_flowrate = -1.0
            raise RuntimeError('injected failure')
        return original(*args, **kwargs)

    hx._solve_2_unknown = flaky
    try:
        solve_quietly(hx)
    except ValueError:
        pass                                      # the second anchor may legitimately refuse
    assert len(seen) == 2, 'both anchors should have been attempted'
    assert seen[0] == entry
    assert seen[1] == entry, 'second anchor started from the failed attempt'


def test_exception_in_post_solve_check_restores_model():
    """A failure after the solver has converged (strict pinch check) must undo
    the solve, not leave a converged but rejected state in the model."""
    M, hx = _overloaded_evaporator(strict=True)
    before = contents(M)
    with pytest.raises(ValueError):
        hx.Cal()
    assert contents(M) == before
    assert not hx.Solution_Status


# --------------------------------------------------------------------------
#  Width of the window excluded next to the anchored end (3-unknown solves)
# --------------------------------------------------------------------------
@pytest.mark.parametrize('frac', [0.02, 0.05, 0.1, 0.2, 0.3])
def test_anchor_exclusion_fraction_does_not_move_the_answer(frac):
    """The profile points nearest the anchored end are left out of the pinch
    search because the anchor already fixes the approach there.  The width of
    that window is a parameter (default 0.1); the recovered flow must be
    insensitive to it over a wide range."""
    ref = REF_3['Evaporator']
    M, hx = build(ref, {'H_hi', 'H_ho', 'mh'})
    hx.anchor_exclude_frac = frac
    w = solve_quietly(hx)
    assert recovered(M)['mh'] == pytest.approx(ref['truth']['mh'], rel=2e-4)
    assert hx.min_dT == pytest.approx(PPT, abs=5e-3)
    assert not w
