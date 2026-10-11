#!/usr/bin/env python3
"""
Dorabella "three" battery: transpositions built on 3, and keys that repeat every 3 symbols.

Needs Dorabella_hillclimb.py, Dorabella_tests.py, Dorabella_marathon.py,
Dorabella_transcription.csv and english.txt in the same folder, plus numpy.

  python Dorabella_three.py --selftest          # ~10-20 min: checks + time estimate
  python Dorabella_three.py                     # start, or resume where it left off
  python Dorabella_three.py --hours 72          # stop by itself after 72 hours
  python Dorabella_three.py --report-only       # rewrite the report from saved state
  python Dorabella_three.py --only T3 PC        # work only on these families
  python Dorabella_three.py --fresh             # discard saved state and start over

Stop any time with Ctrl+C; the same command resumes. Progress is saved after every
batch to three_out/state.json and three_out/report.txt is rewritten each time.
When every target is met it keeps adding scrambled copies until stopped.

THE IDEA
Elgar's alphabet is built on 3 (one, two or three loops), and the message is 87 =
3 x 29 symbols in lines of 29, 31 and 27. Two kinds of "3" scheme are tested.

PRE-REGISTERED FAMILIES (frozen before any real run; 8 families, so the
Bonferroni threshold for a single family is p < 0.05/8 = 0.00625)

  T3  Transpositions built on 3, then ANY one-symbol-one-letter key (annealing):
      read the three lines as a grid down the columns (Elgar's ragged lines, and
      a regular 3 x 29 grid), every row order, columns left-right or right-left,
      snake columns, diagonals; every 3rd symbol; every 3rd symbol within each
      line; the lines in any order, each forwards or backwards; 3-rail fence;
      width-3 columnar with all 6 column orders. Each read both ways
      (as decryption and as encryption). The 5 reading orders already tested
      are left out.
  TB  Broad transpositions beyond 3 (closes the gap in Tim Roberts' widths 2-30,
      which had no scrambled-copy comparison): columnar widths 2-30 (plain and
      reversed column order), all column orders for widths 4 and 5, rail fences
      of 2-8 rails at every offset, every-d-th-symbol decimation for all d,
      spirals of the 3 x 29 grid. Members of T3 are left out.
  TE  Elgar's own published key (no search) read through every T3 and TB order.
  PA  Elgar's published key, then a letter shift that repeats every 3 symbols
      (a Vigenere-style key of length 3): all 24^3 shifts.
  PB  Elgar's published key, but each symbol is first turned (8 rotations),
      flipped (mirror), and/or its loop count stepped (+0/+1/+2 mod 3), with a
      different turn for each of the 3 positions in a cycle: all 48^3 settings.
  PC  ANY one-symbol-one-letter key, then a period-3 letter shift. Every one of
      the 576 shift pairs gets its own key search (a joint search could not find
      the shifts on planted English). Expensive, so forward reading only.
  PD  ANY key, with period-3 symbol turns as in PB. All 2,304 turn settings are
      ranked without a key (the right one makes the symbols look like a simple
      substitution), and the top 30 get a full key search.
  PE  Three completely independent alphabets used in turn (most general).

  PA, PB, PD and PE are run over 10 configurations: the 5 reading orders x the
  3-cycle either running through the whole text or restarting at each line.
  PC uses the 2 forward ones.

METHOD (unchanged): each family's statistic is its best quadgram score over all
its members. The real cipher's statistic is the MEDIAN over several search seeds
(the lesson from the held-out check). Each scrambled copy (plain shuffle, and
pair-preserving shuffle that keeps every adjacent symbol pair) goes through the
IDENTICAL search once; this is conservative (single runs have a wider spread
than a median). p = (1+k)/(n+1). Planted English, enciphered by a random member
of each family, is run through the identical search to show the family can
detect a real message ("power"). Word score is not optimised and gives an
independent check.
"""
import argparse
import csv
import datetime
import itertools
import json
import math
import multiprocessing as mp
import os
import random
import signal
import statistics
import sys
import time

import numpy as np

import Dorabella_hillclimb as H
import Dorabella_tests as T
import Dorabella_marathon as M

VERSION = 1
NL = 24
FAMILIES = ['T3', 'TB', 'TE', 'PA', 'PB', 'PC', 'PD', 'PE']
FAMILY_DESC = {
    'T3': 'transpositions built on 3 + any key',
    'TB': 'broad transpositions (widths 2-30 etc.) + any key',
    'TE': "Elgar's key through every transposition",
    'PA': "Elgar's key + period-3 letter shift",
    'PB': "Elgar's key + period-3 symbol turns",
    'PC': 'any key + period-3 letter shift',
    'PD': 'any key + period-3 symbol turns',
    'PE': 'three independent alphabets in turn',
}
SA_FAMS = {'T3', 'TB', 'PC', 'PD', 'PE'}
ENUM_FAMS = {'TE', 'PA', 'PB'}
MODEL = {'T3': 'masc', 'TB': 'masc', 'PC': 'shift3', 'PD': 'turn3', 'PE': 'alpha3'}
BONF = 0.05 / len(FAMILIES)
KNOWN_ORDERS = ['forward', 'fully reversed', 'each line reversed', 'boustrophedon A', 'boustrophedon B']
LINE_LENS = [29, 31, 27]
LINE_START = [0, 29, 60]

# Elgar's published key (allkeys check: ABCvd3 DEF^d5 GHIvd1 JKL^d7 MNO^d8 PQRvd6 STU^d4 WXYvd2;
# ^ = 1 loop is the first letter, v = 1 loop is the third letter)
ELGAR_SETS = [('ABC', 3, 'v'), ('DEF', 5, '^'), ('GHI', 1, 'v'), ('JKL', 7, '^'),
              ('MNO', 8, '^'), ('PQR', 6, 'v'), ('STU', 4, '^'), ('WXY', 2, 'v')]
STANDARD_DECRYPT = ("BLTACEIARWUNISNFNNELLHSYWYDUOINIEYARQATNNTEDMINUNEHOMSYRRYUO"
                    "TOEHOTSHGDOTNEHMOSALDOEADYA")

QL = None          # flat quadgram table (list) for the annealer
QT = None          # same as numpy array for the enumerations
SEQ = None         # real cipher codes (from the CSV)
FAM = None         # family definitions, built once per process
KEEP_AWAKE_OK = False


# ============================================================== symbols
def slot_of(code):
    return (code // 10 - 1) * 8 + (code % 10 - 1)


def code_of_slot(s):
    return (s // 8 + 1) * 10 + (s % 8 + 1)


def build_elgar():
    E = [None] * 24
    for letters, d, ud in ELGAR_SETS:
        for loop in (1, 2, 3):
            ch = letters[loop - 1] if ud == '^' else letters[3 - loop]
            E[(loop - 1) * 8 + (d - 1)] = H.LIDX[H.fold_base(ch)]
    return E


ELGAR = build_elgar()
ELGAR_INV = [[s for s in range(24) if ELGAR[s] == L] for L in range(24)]


def build_turns():
    """48 symbol transforms: mirror (s=+-1) x rotation (8) x loop step (3). Index 0 = identity."""
    trs, labels = [], []
    for s in (1, -1):
        for r in range(8):
            for m in range(3):
                t = []
                for slot in range(24):
                    y, x = divmod(slot, 8)
                    t.append(((y + m) % 3) * 8 + (s * x + r) % 8)
                trs.append(t)
                labels.append(f"{'mirror+' if s < 0 else ''}rot{r}{'+loop' + str(m) if m else ''}")
    inv = []
    for t in trs:
        iv = [0] * 24
        for a, b in enumerate(t):
            iv[b] = a
        inv.append(iv)
    return trs, inv, labels


TURN, TURN_INV, TURN_LABEL = build_turns()
TURN_LABEL[0] = 'none'


def load_seq(path='Dorabella_transcription.csv'):
    if os.path.exists(path):
        with open(path, newline='', encoding='utf-8') as fh:
            rows = list(csv.DictReader(fh))
        seq = [int(r['code']) for r in sorted(rows, key=lambda r: int(r['position']))]
        assert len(seq) == 87
        return seq, 'Dorabella_transcription.csv'
    seq = list(H.BASE_SEQ)
    seq[31] = 18
    return seq, 'embedded first transcription with 32=18 (CSV not found)'


# ============================================================== reading orders
def inverse(perm):
    iv = [0] * len(perm)
    for i, p in enumerate(perm):
        iv[p] = i
    return iv


def known_order(name):
    return T.apply_order(list(range(87)), name)


def lines_pos():
    return [list(range(a, a + n)) for a, n in zip(LINE_START, LINE_LENS)]


def grid_rows(kind):
    return lines_pos() if kind == 'ragged' else [list(range(r * 29, r * 29 + 29)) for r in range(3)]


def column_read(rows, rev_cols=False, snake=False):
    ncol = max(len(r) for r in rows)
    cols = range(ncol - 1, -1, -1) if rev_cols else range(ncol)
    out = []
    for k, c in enumerate(cols):
        col = [r[c] for r in rows if c < len(r)]
        if snake and k % 2:
            col = col[::-1]
        out.extend(col)
    return out


def diagonal_read(rows, anti, up):
    cells = [(ri, ci) for ri, r in enumerate(rows) for ci in range(len(r))]
    key = (lambda rc: (rc[0] + rc[1], rc[0])) if anti else (lambda rc: (rc[1] - rc[0], rc[0]))
    cells.sort(key=lambda rc: (key(rc)[0], -key(rc)[1] if up else key(rc)[1]))
    return [rows[r][c] for r, c in cells]


def every_kth(seq_positions, k):
    return [seq_positions[i] for off in range(k) for i in range(off, len(seq_positions), k)]


def rail_fence_enc(n, rails, offset, rail_rev=False):
    cyc = 2 * rails - 2
    rail = []
    for i in range(n):
        t = (i + offset) % cyc
        rail.append(t if t < rails else cyc - t)
    order = sorted(range(n), key=lambda i: ((rails - 1 - rail[i]) if rail_rev else rail[i], i))
    return order


def columnar_enc(n, w, colorder):
    """Plain written in rows of width w; cipher = columns in colorder, top to bottom.
    Returns enc with cipher[k] = plain[enc[k]]."""
    out = []
    for c in colorder:
        out.extend(range(c, n, w))
    return out


def spiral(rows, clockwise, corner):
    R, C = len(rows), len(rows[0])
    grid = [[rows[r][c] for c in range(C)] for r in range(R)]
    if corner in (1, 3):
        grid = [g[::-1] for g in grid]
    if corner in (2, 3):
        grid = grid[::-1]
    if not clockwise:
        grid = [list(x) for x in zip(*grid)]
    out = []
    while grid:
        out.extend(grid.pop(0))
        grid = [list(x) for x in zip(*grid)][::-1] if grid else []
    return out


def add_both(store, label, perm):
    """A reading order R (plain[i] = cipher[R[i]]) and its inverse (the same scheme
    used for encryption instead of decryption)."""
    assert sorted(perm) == list(range(87)), label
    store.append((label, list(perm)))
    store.append((label + ' (inverse)', inverse(perm)))


def build_t3():
    out = []
    for kind in ('ragged', 'regular'):
        rows0 = grid_rows(kind)
        for ro in itertools.permutations(range(3)):
            rows = [rows0[i] for i in ro]
            rlab = ''.join(str(i + 1) for i in ro)
            for rev in (False, True):
                for snake in (False, True):
                    add_both(out, f"{kind} grid, rows {rlab}, columns {'R-L' if rev else 'L-R'}"
                             f"{', snake' if snake else ''}", column_read(rows, rev, snake))
        for anti in (False, True):
            for up in (False, True):
                add_both(out, f"{kind} grid, {'anti-' if anti else ''}diagonals "
                         f"{'upward' if up else 'downward'}", diagonal_read(rows0, anti, up))
    allpos = list(range(87))
    for base in ('forward', 'fully reversed'):
        b = known_order(base)
        add_both(out, f"every 3rd symbol ({base})", every_kth(b, 3))
    for rev in (False, True):
        lp = [l[::-1] if rev else l for l in lines_pos()]
        add_both(out, f"every 3rd symbol within each line{' (lines reversed)' if rev else ''}",
                 [p for l in lp for p in every_kth(l, 3)])
    lp = lines_pos()
    for ro in itertools.permutations(range(3)):
        for flips in itertools.product((0, 1), repeat=3):
            perm = [p for i, f in zip(ro, flips) for p in (lp[i][::-1] if f else lp[i])]
            lab = ' '.join(f"{i + 1}{'r' if f else ''}" for i, f in zip(ro, flips))
            add_both(out, f"lines in order {lab}", perm)
    for off in range(4):
        for rr in (False, True):
            for base in ('forward', 'fully reversed'):
                enc = rail_fence_enc(87, 3, off, rr)
                b = known_order(base)
                add_both(out, f"3-rail fence, offset {off}{', rails bottom-up' if rr else ''} ({base})",
                         [b[i] for i in enc])
    for co in itertools.permutations(range(3)):
        add_both(out, f"width-3 columnar, column order {''.join(str(c + 1) for c in co)}",
                 columnar_enc(87, 3, co))
    return dedupe(out, exclude=set())


def build_tb(t3_perms):
    out = []
    for w in range(2, 31):
        for rev in (False, True):
            co = list(range(w))[::-1] if rev else list(range(w))
            add_both(out, f"width-{w} columnar, columns {'reversed' if rev else 'in order'}",
                     columnar_enc(87, w, co))
    for w in (4, 5):
        for co in itertools.permutations(range(w)):
            add_both(out, f"width-{w} columnar, column order {''.join(str(c + 1) for c in co)}",
                     columnar_enc(87, w, co))
    for rails in range(2, 9):
        for off in range(2 * rails - 2):
            for rr in (False, True):
                add_both(out, f"{rails}-rail fence, offset {off}{', rails bottom-up' if rr else ''}",
                         rail_fence_enc(87, rails, off, rr))
    for d in range(2, 86):
        if math.gcd(d, 87) == 1:
            out.append((f"every {d}th symbol, wrapping round", [(i * d) % 87 for i in range(87)]))
    rows = grid_rows('regular')
    for cw in (True, False):
        for corner in range(4):
            add_both(out, f"3x29 spiral, {'clockwise' if cw else 'anticlockwise'}, corner {corner}",
                     spiral(rows, cw, corner))
    return dedupe(out, exclude=t3_perms)


def dedupe(members, exclude):
    known = {tuple(known_order(o)) for o in KNOWN_ORDERS}
    seen, out = set(), []
    for lab, perm in members:
        t = tuple(perm)
        if t in seen or t in known or t in exclude:
            continue
        seen.add(t)
        out.append((lab, perm))
    return out


def p_configs():
    """(label, reading order R, phase list) for the 10 period-3 configurations."""
    out = []
    line_of = [0] * 29 + [1] * 31 + [2] * 27
    for o in KNOWN_ORDERS:
        R = known_order(o)
        out.append((f"{o}, cycle runs through", R, [i % 3 for i in range(87)]))
        cnt = [0, 0, 0]
        ph = []
        for p in R:
            ph.append(cnt[line_of[p]] % 3)
            cnt[line_of[p]] += 1
        out.append((f"{o}, cycle restarts each line", R, ph))
    return out


def build_families():
    t3 = build_t3()
    tb = build_tb({tuple(p) for _, p in t3})
    pc = p_configs()
    return {'T3': t3, 'TB': tb, 'TE': t3 + tb, 'PA': pc, 'PB': pc, 'PC': pc[:2], 'PD': pc, 'PE': pc}


def get_fam():
    global FAM
    if FAM is None:
        FAM = build_families()
    return FAM


# ============================================================== scoring
def init_tables(corpus):
    global QL, QT
    H.init_worker(corpus, None)      # quadgrams + word tables
    QL = [H.QFLOOR] * (NL ** 4)
    for k, v in H.QUAD.items():
        QL[k] = v
    QT = np.array(QL, dtype=np.float64)


def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    init_tables(corpus)
    get_fam()


def score_letters(p):
    Q = QL
    t = 0.0
    for s in range(len(p) - 3):
        t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
    return t / (len(p) - 3)


def np_scores(P):
    idx = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
    return QT[idx].sum(axis=1) / (P.shape[1] - 3)


def text_of(p):
    return ''.join(H.LETTERS[x] for x in p)


def words_of(text):
    return H.word_score(text.translate(H.TOL))


# ============================================================== annealer
SA_T0, SA_T1 = 10.0, 0.3
PSCAN = 400       # steps between parameter scans in the period-3 searches


def _index(e, n, S):
    pos = [[] for _ in range(S)]
    for q, c in enumerate(e):
        pos[c].append(q)
    starts = [tuple(sorted({s for q in pos[c] for s in range(max(0, q - 3), min(q, n - 4) + 1)}))
              for c in range(S)]
    used = [c for c in range(S) if pos[c]]
    return pos, starts, used, {}


def anneal(model, c, ph, restarts, iters, seed, fixed=None):
    """c: symbol slots (0-23) in reading order; ph: position in the 3-cycle (0,1,2).
    model: masc | shift3 | turn3 | alpha3. Returns (avg score, letters, param label)."""
    Q = QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(c)
    B = 3 if model == 'alpha3' else 1
    S = 24 * B
    has_p = model in ('shift3', 'turn3') and fixed is None
    fac = (SA_T1 / SA_T0) ** (1.0 / max(1, iters))
    allst = range(n - 3)
    best_all = (-1e18, None, None)

    def make_e(g):
        if model == 'turn3':
            return [TURN[g[ph[q]]][c[q]] for q in range(n)]
        if model == 'alpha3':
            return [ph[q] * 24 + c[q] for q in range(n)]
        return list(c)

    def full(p):
        t = 0.0
        for s in allst:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    for _ in range(restarts):
        if fixed is not None:
            prm = list(fixed)
        else:
            prm = [0, rr(24), rr(24)] if model == 'shift3' else ([0, rr(48), rr(48)] if model == 'turn3' else None)
        e = make_e(prm)
        sh = [prm[ph[q]] for q in range(n)] if model == 'shift3' else None
        pos, starts, used, pc = _index(e, n, S)
        nu = len(used)
        key = []
        for _b in range(B):
            k = list(range(24))
            rng.shuffle(k)
            key.extend(k)
        p = [(key[e[q]] + sh[q]) % 24 for q in range(n)] if sh else [key[x] for x in e]
        cur = full(p)
        best, bp, bprm = cur, p[:], (prm[:] if prm else None)
        Tm = SA_T0
        for it in range(iters):
            if has_p and it % PSCAN == 0:
                # scan every value of one cycle parameter with the key held fixed,
                # and pick one by the annealing (Boltzmann) rule at this temperature
                k = 1 + (it // PSCAN) % 2
                nv = 24 if model == 'shift3' else 48
                cands = []
                for v in range(nv):
                    prm[k] = v
                    if model == 'shift3':
                        sh2 = [prm[ph[q]] for q in range(n)]
                        p2 = [(key[e[q]] + sh2[q]) % 24 for q in range(n)]
                        cands.append((full(p2), v, p2, sh2))
                    else:
                        e2 = make_e(prm)
                        p2 = [key[x] for x in e2]
                        cands.append((full(p2), v, p2, e2))
                top = max(x[0] for x in cands)
                w = [exp((x[0] - top) / Tm) for x in cands]
                r = rand() * sum(w)
                for x, wx in zip(cands, w):
                    r -= wx
                    if r <= 0:
                        break
                s2, v, p2, aux = x
                prm[k] = v
                cur, p = s2, p2
                if model == 'shift3':
                    sh = aux
                else:
                    e = aux
                    pos, starts, used, pc = _index(e, n, S)
                    nu = len(used)
                if cur > best:
                    best, bp, bprm = cur, p[:], prm[:]
                Tm *= fac
                continue
            i = used[rr(nu)]
            blk = i - i % 24
            j = blk + rr(23)
            if j >= i:
                j += 1
            a, b = (i, j) if i < j else (j, i)
            ss = pc.get(a * S + b)
            if ss is None:
                ss = tuple(sorted(set(starts[a]) | set(starts[b])))
                pc[a * S + b] = ss
            o = 0.0
            for s in ss:
                o += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
            li, lj = key[i], key[j]
            if sh:
                for q in pos[i]:
                    p[q] = (lj + sh[q]) % 24
                for q in pos[j]:
                    p[q] = (li + sh[q]) % 24
            else:
                for q in pos[i]:
                    p[q] = lj
                for q in pos[j]:
                    p[q] = li
            nw = 0.0
            for s in ss:
                nw += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
            d = nw - o
            if d >= 0 or rand() < exp(d / Tm):
                key[i], key[j] = lj, li
                cur += d
                if cur > best:
                    best, bp, bprm = cur, p[:], (prm[:] if prm else None)
            else:
                if sh:
                    for q in pos[i]:
                        p[q] = (li + sh[q]) % 24
                    for q in pos[j]:
                        p[q] = (lj + sh[q]) % 24
                else:
                    for q in pos[i]:
                        p[q] = li
                    for q in pos[j]:
                        p[q] = lj
            Tm *= fac
        assert abs(best - full(bp)) < 1e-6, "annealer drifted"
        if best > best_all[0]:
            best_all = (best, bp, bprm)
    s, bp, bprm = best_all
    return s / (n - 3), bp, param_label(model, bprm)


def param_label(model, prm):
    if prm is None:
        return ''
    if model in ('shift3', 'PA'):
        return 'shifts ' + ','.join(str(x) for x in prm)
    return 'turns ' + ' | '.join(TURN_LABEL[g] for g in prm)


# ============================================================== one family search
def member_inputs(fam, seq, mi):
    """Symbol slots in reading order and phases for member mi of a family."""
    F = get_fam()[fam]
    if fam in ('T3', 'TB', 'TE'):
        R = F[mi][1]
        return [slot_of(seq[r]) for r in R], [0] * 87
    _, R, ph = F[mi]
    return [slot_of(seq[r]) for r in R], ph


def sa_job(fam, seq, mi, cfg, seed):
    c, ph = member_inputs(fam, seq, mi)
    if fam == 'PC':
        s, p, lab = pc_search(c, ph, cfg, seed)
    elif fam == 'PD':
        s, p, lab = pd_search(c, ph, cfg, seed)
    else:
        s, p, lab = anneal(MODEL[fam], c, ph, cfg[0], cfg[1], seed)
    return s, text_of(p), lab


def pc_search(c, ph, cfg, seed):
    """Any key + period-3 letter shift. A joint search cannot find the shifts (tested on
    planted English), so every one of the 24 x 24 shift pairs gets its own key search;
    the best few are then searched again more deeply."""
    r1, i1, keep, r2, i2 = cfg
    res = []
    for a in range(24):
        for b in range(24):
            s, _, _ = anneal('shift3', c, ph, r1, i1, derive_seed(seed, a, b), fixed=[0, a, b])
            res.append((s, a, b))
    res.sort(reverse=True)
    best = (-1e18, None, None)
    for s0, a, b in res[:keep]:
        s, p, lab = anneal('shift3', c, ph, r2, i2, derive_seed(seed, 'polish', a, b), fixed=[0, a, b])
        if s > best[0]:
            best = (s, p, lab)
    return best


def screen_turns(c, ph):
    """Rank all 48 x 48 turn settings without any key: the right setting makes the
    symbols behave like one simple substitution (higher index of coincidence, more
    repeated pairs and triples). Tested: the true setting ranks in the top 11."""
    out = []
    n = len(c)
    for g1 in range(48):
        for g2 in range(48):
            g = (0, g1, g2)
            e = [TURN[g[ph[q]]][c[q]] for q in range(n)]
            cnt = [0] * 24
            for x in e:
                cnt[x] += 1
            ic = sum(v * (v - 1) for v in cnt) / (n * (n - 1))
            seen2, seen3, rep_ = set(), set(), 0
            for q in range(n - 1):
                t2 = e[q] * 24 + e[q + 1]
                if t2 in seen2:
                    rep_ += 1
                seen2.add(t2)
                if q < n - 2:
                    t3 = t2 * 24 + e[q + 2]
                    if t3 in seen3:
                        rep_ += 2
                    seen3.add(t3)
            out.append((ic * 100 + rep_ * 0.3, g1, g2))
    out.sort(reverse=True)
    return out


def pd_search(c, ph, cfg, seed):
    keep, r2, i2 = cfg
    best = (-1e18, None, None)
    for _, g1, g2 in screen_turns(c, ph)[:keep]:
        s, p, lab = anneal('turn3', c, ph, r2, i2, derive_seed(seed, g1, g2), fixed=[0, g1, g2])
        if s > best[0]:
            best = (s, p, lab)
    return best


def enum_job(fam, seq):
    """Exhaustive families. Returns per-member list of (score, letters text, param label)."""
    F = get_fam()[fam]
    E = np.array(ELGAR, dtype=np.int64)
    out = []
    if fam == 'TE':
        P = np.array([[ELGAR[slot_of(seq[r])] for r in R] for _, R in F], dtype=np.int64)
        sc = np_scores(P)
        return [(float(sc[i]), text_of(P[i].tolist()), '') for i in range(len(F))]
    if fam == 'PA':
        S3 = np.array(list(itertools.product(range(24), repeat=3)), dtype=np.int64)
        for mi in range(len(F)):
            c, ph = member_inputs(fam, seq, mi)
            L = E[np.array(c)]
            P = (L[None, :] + S3[:, np.array(ph)]) % 24
            sc = np_scores(P)
            b = int(sc.argmax())
            out.append((float(sc[b]), text_of(P[b].tolist()), param_label('PA', S3[b].tolist())))
        return out
    if fam == 'PB':
        LT = E[np.array(TURN, dtype=np.int64)]            # LT[g, slot] = letter
        G12 = np.array(list(itertools.product(range(48), repeat=2)), dtype=np.int64)
        for mi in range(len(F)):
            c, ph = member_inputs(fam, seq, mi)
            c = np.array(c)
            ph = np.array(ph)
            cols = [np.where(ph == k)[0] for k in range(3)]
            best = (-1e18, None, None)
            P = np.empty((len(G12), 87), dtype=np.int64)
            for g0 in range(48):
                P[:, cols[0]] = LT[g0, c[cols[0]]][None, :]
                P[:, cols[1]] = LT[G12[:, 0]][:, c[cols[1]]]
                P[:, cols[2]] = LT[G12[:, 1]][:, c[cols[2]]]
                sc = np_scores(P)
                b = int(sc.argmax())
                if sc[b] > best[0]:
                    best = (float(sc[b]), P[b].tolist(), [g0] + G12[b].tolist())
            out.append((best[0], text_of(best[1]), param_label('turn3', best[2])))
        return out
    raise ValueError(fam)


# ============================================================== controls
def plant(fam, idx, noise, seed):
    """Known English enciphered by a random member/setting of this family.
    Returns (codes in written order, plaintext, true member index)."""
    rng = random.Random(seed)
    texts = M.CONTROL_TEXTS
    raw = texts[idx % len(texts)] + ' ' + texts[(idx + 7) % len(texts)]
    plain = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))[:87]
    L = [H.LIDX[ch] for ch in plain]
    F = get_fam()[fam]
    mi = rng.randrange(len(F))
    if fam in ('T3', 'TB', 'TE'):
        R = F[mi][1]
        ph = [0] * 87
    else:
        _, R, ph = F[mi]

    def rand_key():
        k = list(range(24))
        rng.shuffle(k)
        return k, inverse(k)

    slots = []
    prm = None
    if fam in ('T3', 'TB'):
        _, kinv = rand_key()
        slots = [kinv[x] for x in L]
    elif fam == 'TE':
        slots = [rng.choice(ELGAR_INV[x]) if ELGAR_INV[x] else rng.randrange(24) for x in L]
    elif fam in ('PA', 'PC'):
        s = [0 if fam == 'PC' else rng.randrange(24), rng.randrange(24), rng.randrange(24)]
        prm = s
        if fam == 'PA':
            slots = [rng.choice(ELGAR_INV[(x - s[k]) % 24]) if ELGAR_INV[(x - s[k]) % 24]
                     else rng.randrange(24) for x, k in zip(L, ph)]
        else:
            _, kinv = rand_key()
            slots = [kinv[(x - s[k]) % 24] for x, k in zip(L, ph)]
    elif fam in ('PB', 'PD'):
        g = [0 if fam == 'PD' else rng.randrange(48), rng.randrange(48), rng.randrange(48)]
        prm = g
        if fam == 'PB':
            slots = [TURN_INV[g[k]][rng.choice(ELGAR_INV[x])] if ELGAR_INV[x] else rng.randrange(24)
                     for x, k in zip(L, ph)]
        else:
            _, kinv = rand_key()
            slots = [TURN_INV[g[k]][kinv[x]] for x, k in zip(L, ph)]
    elif fam == 'PE':
        kinvs = [rand_key()[1] for _ in range(3)]
        slots = [kinvs[k][x] for x, k in zip(L, ph)]
    seq = [0] * 87
    for i, r in enumerate(R):
        seq[r] = code_of_slot(slots[i])
    for q in rng.sample(range(87), noise):            # misread a symbol by one direction
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, mi, prm


# ============================================================== utilities
def now():
    return datetime.datetime.now()


def stamp():
    return f"[{now():%Y-%m-%d %H:%M:%S}]"


def derive_seed(*parts):
    return random.Random('|'.join(str(x) for x in parts)).getrandbits(48)


def pval(real, nulls):
    if not nulls:
        return None
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


def keep_awake():
    global KEEP_AWAKE_OK
    if sys.platform == 'win32':
        try:
            import ctypes
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
            KEEP_AWAKE_OK = True
        except Exception:
            pass


def run_jobs(pool, fn, argsets):
    if pool is None:
        return [fn(*a) for a in argsets]
    res = pool.starmap_async(fn, argsets, chunksize=1)
    while True:
        try:
            return res.get(60)
        except mp.TimeoutError:
            continue


def save_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def write_text(path, text):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        fh.write(text)
    os.replace(tmp, path)


# ============================================================== running families
def family_search(pool, items, P):
    """items: list of (fam, seq, seed). Runs the full family search on each and
    returns, per item, the per-member list of (score, text, param)."""
    jobs, owners = [], []
    for k, (fam, seq, seed) in enumerate(items):
        if fam in ENUM_FAMS:
            jobs.append((enum_job, (fam, seq)))
            owners.append((k, None))
        else:
            cfg = P['sa'][fam if fam in ('PC', 'PD') else MODEL[fam]]
            for mi in range(len(get_fam()[fam])):
                jobs.append((sa_job, (fam, seq, mi, cfg, derive_seed(seed, fam, mi))))
                owners.append((k, mi))
    results = [[] for _ in items]
    # group by function so starmap gets one callable
    for fn in (enum_job, sa_job):
        idx = [i for i, (f, _) in enumerate(jobs) if f is fn]
        if not idx:
            continue
        out = run_jobs(pool, fn, [jobs[i][1] for i in idx])
        for i, r in zip(idx, out):
            k, mi = owners[i]
            if mi is None:
                results[k] = r
            else:
                results[k].append((mi, r))
    final = []
    for k, (fam, _, _) in enumerate(items):
        r = results[k]
        if fam in SA_FAMS:
            r = [x for _, x in sorted(r, key=lambda t: t[0])]
        final.append(r)
    return final


def summarize_members(members):
    scores = [m[0] for m in members]
    b = max(range(len(scores)), key=lambda i: scores[i])
    return dict(stat=scores[b], best=b, text=members[b][1], param=members[b][2],
                words=words_of(members[b][1]))


# ============================================================== state
DEFAULT_TARGETS = dict(
    T3=dict(shuffle=1000, doublet=1000, controls=40),
    TB=dict(shuffle=300, doublet=300, controls=40),
    TE=dict(shuffle=2000, doublet=2000, controls=40),
    PA=dict(shuffle=2000, doublet=2000, controls=40),
    PB=dict(shuffle=1000, doublet=1000, controls=40),
    PC=dict(shuffle=250, doublet=250, controls=40),
    PD=dict(shuffle=500, doublet=500, controls=40),
    PE=dict(shuffle=1000, doublet=1000, controls=40),
)


def blank_state(args):
    fam = get_fam()
    return dict(version=VERSION, created=str(now()), elapsed=0.0, sessions=0,
                params=dict(seed=args.seed, real_seeds=args.real_seeds,
                            sa={'masc': [args.masc_restarts, args.masc_iters],
                                'PC': [args.pc_screen_restarts, args.pc_screen_iters, args.keep,
                                       args.masc_restarts, args.masc_iters],
                                'PD': [args.pd_keep, args.masc_restarts, args.masc_iters],
                                'alpha3': [args.e_restarts, args.e_iters]}),
                sizes={f: len(fam[f]) for f in FAMILIES},
                targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                fam={f: dict(real=[], shuffle=[], doublet=[], controls=[]) for f in FAMILIES},
                unit_secs={}, log=[])


def stage_list(state, only):
    fams = only or FAMILIES
    out = []
    for f in fams:
        st = state['fam'][f]
        tg = state['targets'][f]
        nreal = 1 if f in ENUM_FAMS else state['params']['real_seeds']
        if len(st['real']) < nreal:
            out.append((0, -1.0, f, 'real'))
            continue
        # early controls first (power check), then nulls and controls by fraction done
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], f, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            frac = len(st[kind]) / tg[kind]
            if frac < 1:
                out.append((2, frac, f, kind))
    return sorted(out)


def extend_targets(state, only):
    for f in (only or FAMILIES):
        for kind in ('shuffle', 'doublet'):
            state['targets'][f][kind] = int(state['targets'][f][kind] * 1.5)



_BIGSEG = []


def show_words(text):
    """Words in a decrypt, for DISPLAY ONLY (scores never use this): the big ~11,700-word list
    built from english.txt (as in the styles test), longer words preferred; falls back to the
    small built-in list if that cannot be built."""
    if not _BIGSEG:
        try:
            import Dorabella_styles as S
            d = S.StyleDict(S.base_dictionary('english.txt'), 'S0')
            _BIGSEG.append(lambda t: d.segment(S.fold(t)))
        except Exception:
            _BIGSEG.append(lambda t: H.segment(t.translate(H.TOL)))
    return _BIGSEG[0](text)

# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def pct(vals, q):
    if not vals:
        return None
    v = sorted(vals)
    return v[min(len(v) - 1, int(q * len(v)))]


def real_summary(st):
    stats = [r['stat'] for r in st['real']]
    if not stats:
        return None, None
    med = statistics.median(stats)
    best = max(st['real'], key=lambda r: r['stat'])
    return med, dict(seed_stats=stats, best=best)


def verdict(p_d, p_s, power):
    if p_d is None:
        return 'waiting for copies'
    if p_d <= BONF and (p_s or 1) <= BONF:
        if power is not None and power < 0.5:
            return 'BEATS both nulls, but controls show low power: check by eye'
        return 'BEATS both nulls past the Bonferroni line: read the decrypts by eye'
    if p_d <= 0.05:
        return 'below 0.05 but not past Bonferroni (expected by chance in 8 families)'
    if power is not None and power < 0.5:
        return 'chance level, but controls show LOW POWER: inconclusive, not a refutation'
    return 'chance level'


def build_report(state):
    fam = get_fam()
    L = []
    L.append("DORABELLA 'THREE' REPORT")
    L.append(f"written {now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
             f"sessions {state['sessions']}")
    sa = state['params']['sa']
    L.append(f"settings: real seeds {state['params']['real_seeds']} (median used), annealing "
             f"restarts x steps: any-key {sa['masc'][0]}x{sa['masc'][1]}, PC screen "
             f"{sa['PC'][0]}x{sa['PC'][1]} per shift pair (best {sa['PC'][2]} re-searched), PD top "
             f"{sa['PD'][0]} turn settings, three alphabets {sa['alpha3'][0]}x{sa['alpha3'][1]}")
    L.append(f"cipher: {state.get('source', '')}")
    L.append(f"Bonferroni line for one of 8 families: p < {BONF:.5f}")
    L.append("")
    L.append("=" * 100)
    L.append("HEADLINE  (statistic = best quadgram score in the family; real = median over seeds;")
    L.append("           power = share of planted-English controls above the 95th percentile of the")
    L.append("           pair-preserving copies; p floor = 1/(n+1))")
    L.append("=" * 100)
    L.append(f"{'fam':<4} {'members':>7} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} "
             f"{'p_shuf':>7} {'p_words':>7} {'power':>6} {'ctrl med':>8}  verdict")
    rows_detail = []
    for f in FAMILIES:
        st = state['fam'][f]
        med, rs = real_summary(st)
        dnull = [x['stat'] for x in st['doublet']]
        snull = [x['stat'] for x in st['shuffle']]
        p_d = pval(med, dnull) if med is not None else None
        p_s = pval(med, snull) if med is not None else None
        wd = [x['words'] for x in st['doublet']]
        real_words = rs['best']['words'] if rs else None
        p_w = pval(real_words, wd) if (rs and wd) else None
        ctrl = [c['stat'] for c in st['controls']]
        q95 = pct(dnull, 0.95)
        power = (sum(1 for c in ctrl if c > q95) / len(ctrl)) if (ctrl and q95 is not None) else None
        L.append(f"{f:<4} {len(fam[f]):>7} {fmt(med, 3):>8} {len(dnull):>6} {fmt(p_d, 4):>7} "
                 f"{len(snull):>6} {fmt(p_s, 4):>7} {fmt(p_w, 3):>7} "
                 f"{('-' if power is None else f'{power:.0%}'):>6} "
                 f"{fmt(statistics.median(ctrl) if ctrl else None, 2):>8}  {verdict(p_d, p_s, power)}")
        rows_detail.append((f, st, med, rs, dnull, snull, ctrl, q95))
    L.append("")
    L.append("Reading it: a family only matters if the real cipher beats the PAIR-PRESERVING copies")
    L.append("past the Bonferroni line, the controls show the family can find a real message (power")
    L.append("near 100%), and the winning decrypt reads as English. Beating plain shuffles alone is")
    L.append("expected (the Dorabella has strong pair habits).")
    for f, st, med, rs, dnull, snull, ctrl, q95 in rows_detail:
        L.append("")
        L.append("-" * 100)
        L.append(f"{f}: {FAMILY_DESC[f]}   ({len(fam[f])} members)")
        L.append("-" * 100)
        if rs:
            ss = rs['seed_stats']
            L.append(f"real statistic per seed: {', '.join(f'{x:.3f}' for x in ss)}   median {med:.3f}")
            for kind, nul in (('pair-preserving', dnull), ('plain shuffle', snull)):
                if nul:
                    L.append(f"  vs {kind:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, "
                             f"95th pct {pct(nul, 0.95):.3f}, max {max(nul):.3f}, p = {pval(med, nul):.4f}")
            b = rs['best']
            L.append(f"best real decrypt (seed with the highest statistic): score {b['stat']:.3f}, "
                     f"words {b['words']}")
            for t in b.get('top', []):
                lab = member_label(f, t['member'])
                L.append(f"  [{t['score']:.3f}] {lab}{('  ' + t['param']) if t['param'] else ''}")
                for a in range(0, 87, 29):
                    L.append("      " + t['text'][a:a + 29])
                L.append("      segmented: " + show_words(t['text']))
            if dnull:
                cnt = {}
                for x in st['doublet']:
                    cnt[x['best']] = cnt.get(x['best'], 0) + 1
                top = sorted(cnt.items(), key=lambda kv: -kv[1])[:3]
                L.append("members most often winning in pair-preserving copies: " +
                         '; '.join(f"{member_label(f, m)} ({k}x)" for m, k in top))
        else:
            L.append("real run not done yet")
        if ctrl:
            for nz in sorted({c['noise'] for c in st['controls']}):
                cs = [c for c in st['controls'] if c['noise'] == nz]
                above = sum(1 for c in cs if q95 is not None and c['stat'] > q95)
                L.append(f"controls, {nz} misread symbols: n={len(cs)}, median stat "
                         f"{statistics.median(c['stat'] for c in cs):.3f}, true member found "
                         f"{sum(c['found'] for c in cs)}/{len(cs)}, median letters right "
                         f"{statistics.median(c['acc'] for c in cs):.0%}, above null 95th pct "
                         f"{above}/{len(cs) if q95 is not None else '-'}")
            ex = max(st['controls'], key=lambda c: c['acc'])
            L.append(f"  example control decrypt ({ex['acc']:.0%} right): {ex['text'][:60]}...")
    L.append("")
    L.append("Run log (last 20):")
    L.extend("  " + x for x in state['log'][-20:])
    return '\n'.join(L) + '\n'


def member_label(f, mi):
    F = get_fam()[f]
    return F[mi][0] if mi is not None and mi < len(F) else '?'


# ============================================================== units of work
def do_batch(pool, state, P, f, kind, count):
    st = state['fam'][f]
    seed0 = state['params']['seed']
    items, meta = [], []
    for k in range(count):
        if kind == 'real':
            n = len(st['real']) + k
            items.append((f, SEQ, derive_seed(seed0, 'real', f, n)))
            meta.append(n)
        elif kind in ('shuffle', 'doublet'):
            n = len(st[kind]) + k
            rng = random.Random(derive_seed(seed0, 'null', f, kind, n))
            if kind == 'shuffle':
                s = list(SEQ)
                rng.shuffle(s)
            else:
                s = M.doublet_shuffle(list(SEQ), rng)
            items.append((f, s, derive_seed(seed0, 'nullsearch', f, kind, n)))
            meta.append(n)
        else:
            n = len(st['controls']) + k
            noise = 0 if n % 2 == 0 else 2
            seq, plain, mi, _prm = plant(f, n, noise, derive_seed(seed0, 'plant', f, n))
            items.append((f, seq, derive_seed(seed0, 'ctrlsearch', f, n)))
            meta.append((n, noise, plain, mi))
    t0 = time.time()
    res = family_search(pool, items, P)
    dt = time.time() - t0
    for members, m in zip(res, meta):
        sm = summarize_members(members)
        if kind == 'real':
            order = sorted(range(len(members)), key=lambda i: -members[i][0])[:5]
            sm['top'] = [dict(member=i, score=members[i][0], text=members[i][1], param=members[i][2])
                         for i in order]
            sm['all'] = [round(x[0], 4) for x in members]
            st['real'].append(sm)
        elif kind in ('shuffle', 'doublet'):
            st[kind].append(dict(stat=sm['stat'], best=sm['best'], words=sm['words']))
        else:
            n, noise, plain, mi = m
            acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
            st['controls'].append(dict(stat=sm['stat'], noise=noise, found=int(sm['best'] == mi),
                                       acc=acc, words=sm['words'], text=sm['text']))
    key = f"{f}:{'real' if kind == 'real' else 'copy'}"
    state['unit_secs'][key] = dt / count
    return dt


def batch_size(state, f, kind, workers):
    if kind == 'real':
        need = (1 if f in ENUM_FAMS else state['params']['real_seeds']) - len(state['fam'][f]['real'])
        return max(1, need)
    if f in ENUM_FAMS:
        n = workers * 2
    else:
        n = max(1, math.ceil(3 * workers / len(get_fam()[f])))
    tg = state['targets'][f][kind]
    return max(1, min(n, tg - len(state['fam'][f][kind])))


# ============================================================== selftest
def selftest(args, pool, P, workers):
    print(f"\n{stamp()} SELF-TEST")
    pub = [H.LIDX[c] for c in H.fold_base(STANDARD_DECRYPT)]
    mine = [ELGAR[slot_of(c)] for c in SEQ]
    agree = sum(a == b for a, b in zip(mine, pub))
    print(f"Elgar's key on this transcription agrees with the published decrypt at {agree}/87 positions")
    print("  " + text_of(mine))
    assert agree >= 75, "Elgar key check failed"
    fam = get_fam()
    for f in FAMILIES:
        print(f"  {f}: {len(fam[f]):>4} members  ({FAMILY_DESC[f]})")
    # transformation sanity
    for t, iv in zip(TURN, TURN_INV):
        assert sorted(t) == list(range(24)) and all(iv[t[s]] == s for s in range(24))
    d = M.doublet_shuffle(list(SEQ), random.Random(1))
    print(f"pair-preserving shuffle check passed ({''.join(str(x) for x in d[:6])}...)")
    print("\nOne planted-English control per family (no misreads), through the full family search:")
    est = {}
    for f in (args.only or FAMILIES):
        seq, plain, mi, _prm = plant(f, 3, 0, derive_seed('selftest', f))
        t0 = time.time()
        members = family_search(pool, [(f, seq, derive_seed('selftest-search', f))], P)[0]
        dt = time.time() - t0
        sm = summarize_members(members)
        acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
        # a null copy for scale
        s = M.doublet_shuffle(list(SEQ), random.Random(derive_seed('selftest-null', f)))
        nm = summarize_members(family_search(pool, [(f, s, 5)], P)[0])
        est[f] = dt
        print(f"  {f}: control stat {sm['stat']:.3f} (one pair-preserving copy: {nm['stat']:.3f}), "
              f"true member found: {'yes' if sm['best'] == mi else 'NO'}, letters right {acc:.0%}, "
              f"{H.fmt_secs(dt)} per copy")
        print(f"      planted: {plain[:58]}")
        print(f"      decoded: {sm['text'][:58]}")
    tot = 0.0
    print("\nEstimated time to reach the targets on this machine:")
    for f in (args.only or FAMILIES):
        tg = DEFAULT_TARGETS[f]
        n = tg['shuffle'] + tg['doublet'] + tg['controls'] + (1 if f in ENUM_FAMS else args.real_seeds)
        jobs = 1 if f in ENUM_FAMS else len(fam[f])
        per = est[f] * min(workers, jobs) / workers   # small families run copies side by side
        t = n * per
        tot += t
        print(f"  {f}: {n} searches x {H.fmt_secs(per)} = {H.fmt_secs(t)}")
    print(f"  total about {H.fmt_secs(tot)} (it keeps going after that until you stop it)")


# ============================================================== main
def main():
    global SEQ
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='three_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0, help='stop after this many hours (0 = never)')
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--real-seeds', type=int, default=5, help='search seeds on the real cipher (median used)')
    ap.add_argument('--masc-restarts', type=int, default=8)
    ap.add_argument('--masc-iters', type=int, default=15000)
    ap.add_argument('--pc-screen-restarts', type=int, default=4)
    ap.add_argument('--pc-screen-iters', type=int, default=10000)
    ap.add_argument('--keep', type=int, default=24, help='PC: shift pairs searched again deeply')
    ap.add_argument('--pd-keep', type=int, default=30, help='PD: turn settings searched after the screen')
    ap.add_argument('--e-restarts', type=int, default=12)
    ap.add_argument('--e-iters', type=int, default=40000)
    ap.add_argument('--only', nargs='+', choices=FAMILIES)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    SEQ, source = load_seq()
    print(f"cipher read from: {source}")
    if 'not found' in source:
        print("WARNING: Dorabella_transcription.csv was not found in this folder, so the older "
              "embedded transcription is being used (position 85 = 14). Put the CSV next to this "
              "script, named exactly Dorabella_transcription.csv.")
    os.makedirs(args.out, exist_ok=True)
    state_path = os.path.join(args.out, 'state.json')
    report_path = os.path.join(args.out, 'report.txt')
    log_path = os.path.join(args.out, 'log.txt')

    def log(msg):
        line = f"{stamp()} {msg}"
        print(line, flush=True)
        with open(log_path, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    init_tables(args.corpus)
    get_fam()

    if os.path.exists(state_path) and not args.fresh:
        with open(state_path, encoding='utf-8') as fh:
            state = json.load(fh)
        if state.get('sizes') != {f: len(FAM[f]) for f in FAMILIES}:
            raise SystemExit("saved state was made with different families; use --fresh or another --out")
        if state.get('seq') is not None and state['seq'] != SEQ:
            raise SystemExit("saved state was made with a DIFFERENT transcription than the one read now "
                             f"({source}). Use --fresh (or another --out) to start over.")
        if state.get('seq') is None and state['sessions'] > 0 and not args.report_only:
            raise SystemExit("saved state does not record which transcription it used (made by an older "
                             "version of this script). Use --fresh to start over.")
    else:
        state = blank_state(args)
    state['source'] = source
    state['seq'] = SEQ
    P = state['params']

    if args.report_only:
        write_text(report_path, build_report(state))
        print(f"report written to {report_path}")
        return

    workers = max(1, args.workers)
    keep_awake()
    pool = mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, P, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}, cipher from {source}"
            f"{', keep-awake on' if KEEP_AWAKE_OK else ''}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stages = stage_list(state, args.only)
            if not stages:
                extend_targets(state, args.only)
                log("all targets met; raising the null targets by 50% and carrying on")
                continue
            _, _, f, kind = stages[0]
            n = batch_size(state, f, kind, workers)
            dt = do_batch(pool, state, P, f, kind, n)
            now_t = time.time()
            state['elapsed'] += now_t - last
            last = now_t
            st = state['fam'][f]
            msg = f"{f} {kind} x{n} in {H.fmt_secs(dt)}"
            if kind == 'real':
                msg += f"  (real stat {st['real'][-1]['stat']:.3f})"
            else:
                msg += f"  ({len(st[kind])}/{state['targets'][f][kind]})"
            state['log'].append(msg)
            state['log'] = state['log'][-200:]
            print(f"{stamp()} {msg}", flush=True)
            save_json(state_path, state)
            write_text(report_path, build_report(state))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            save_json(state_path, state)
            write_text(report_path, build_report(state))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    mp.freeze_support()
    main()
