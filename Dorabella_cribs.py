#!/usr/bin/env python3
"""
Dorabella crib test: do words Elgar would likely have written to Dora FIT the cipher,
and does fixing them make the REST of the cipher read as English?

Needs Dorabella_styles.py, Dorabella_allkeys.py, Dorabella_hillclimb.py,
Dorabella_marathon.py, Dorabella_tests.py, Dorabella_transcription.csv, english.txt and
numpy in the same folder.

  python Dorabella_cribs.py --selftest        # ~20-30 min: eligible cribs, power, timing
  python Dorabella_cribs.py                   # start, or resume where it left off
  python Dorabella_cribs.py --hours 48        # stop by itself after 48 hours
  python Dorabella_cribs.py --report-only     # rewrite the report from saved state
  python Dorabella_cribs.py --fresh           # discard saved state and start over

Ctrl+C any time; the same command resumes. Output: cribs_out/report.txt, log.txt,
state.json. When every target is met it raises the copy targets by 50% and carries on.

THE IDEA
Under any one-symbol-one-letter key a word leaves a fingerprint: its pattern of
repeated letters. DEAR DORA can only sit where the cipher repeats symbols in the same
pattern. So each crib is placed wherever it fits, its letters are fixed, the rest of
the key is searched, and the result is scored. The same is done to scrambled copies,
so the length of the crib list is automatically accounted for.

PRE-REGISTERED TESTS (frozen before any run)
  C1  ANY KEY. Statistic: the best English score (quadgrams, standard spelling) of the
      cipher completed around ANY eligible crib in ANY of the 5 reading orders, with
      the 9 uncertain symbols free. Each copy gets the identical search.
  C2  ELGAR-STYLE GROUPED KEYS (the cribber's spreadsheet approach): the key is one
      of the 10,321,920 "sets of three" keys (each direction holds one set of three
      letters, read up or down; both letter groupings). Statistic: the best score
      among grouped keys that spell an eligible crib somewhere. Main readings only.
  C3  LOOP-GROUP KEYS (the cribber's "SortLoops" approach): the loop count picks a run
      of 8 letters (A-H, I-Q, R-Z) and the direction picks the letter; every grouping
      order, reversal and rotation (24,576 keys). Statistic as for C2.
  H1  "MY DEAR DORA" - FOUND BY LOOKING, reported separately and labelled as such.
      Statistic: the best C1-style completion around DEARDORA or MYDEARDORA. The
      real cipher has one fit (start of line 3, read backwards); about 1 copy in 40
      has a fit somewhere. Copies with no fit count as failures.

ELIGIBILITY RULE (fixed before any run, decided on scrambled copies only, never on the
real cipher): a crib is used if it has at least 5 letters and fits, on average, no more
than 20 places in a pair-preserving copy (5 reading orders). Short or patternless
cribs (DORA, ELGAR, ALICE...) fit everywhere and cannot narrow anything down.

METHOD: real = median over several search seeds; plain shuffles and pair-preserving
shuffles get the identical search; planted English containing a crib measures power.
p = (1+k)/(n+1).
"""
import argparse
import datetime
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

import Dorabella_styles as S
import Dorabella_allkeys as AK
import Dorabella_hillclimb as H
import Dorabella_marathon as M

VERSION = 1
LET = H.LETTERS
LIDX = H.LIDX
ORDERS = S.ORDERS
ELIG_MIN_LEN, ELIG_MAX_FITS = 5, 20
H1_CRIBS = ('DEARDORA', 'MYDEARDORA')

# ---------------------------------------------------------------- the crib list
CRIBS = {
    'greetings and sign-offs': """MY DEAR DORA | DEAR DORA | DEAREST DORA | DEAR MISS PENNY |
        MISS PENNY | MY DEAR | DEAREST | YOURS VERY TRULY | YOURS SINCERELY | YOURS EVER |
        EVER YOURS | AFFECTIONATELY | BELIEVE ME | KIND REGARDS | GOODBYE | THANK YOU |
        THANK YOU SO MUCH | MANY THANKS""",
    'the visit': """DELIGHTFUL | TOMORROW | YESTERDAY | BICYCLE | VISIT | PUZZLE | CIPHER | CYPHER |
        WOLVERHAMPTON | RECTORY | THE RECTORY | HIGHFIELD | HASFIELD | FORLI | MALVERN |
        WORCESTER | BROADHEATH | LONDON | JULY""",
    'the Penny family': """DORA | DORA PENNY | DORA MARY PENNY | DORABELLA | PENNY | ALFRED |
        ALFRED PENNY | MRS PENNY | MARY FRANCES | HEALE | BAKER | MEATH BAKER |
        WILLIAM MEATH BAKER | TOWNSHEND""",
    'the Elgar family': """EDWARD | EDWARD ELGAR | ELGAR | EDOO | ALICE | MRS ELGAR | CAROLINE |
        CAROLINE ALICE | ROBERTS | CARICE | WILLIAM HENRY | ANN GREENING | GREENING |
        HARRY | LUCY | POLLIE | FRANK | HELEN AGNES | FREDERICK""",
    'pets': """MARCO | MEG | MINA | DAN""",
    "Elgar's circle": """JAEGER | NIMROD | AUGUST | SINCLAIR | TROYTE | GRIFFITH | YSOBEL |
        ISABEL | FITTON | WINIFRED | NORBURY | NEVINSON | STEUART POWELL | ARNOLD | LYGON""",
    'the 1920 notebook': """MARCO ELGAR | A VERY OLD CYPHER | VERY OLD | DO YOU GO TO LONDON TOMORROW |
        LONDON TOMORROW | GO TO LONDON""",
}


def crib_list():
    out, seen = [], set()
    for group, txt in CRIBS.items():
        for item in txt.split('|'):
            raw = ' '.join(item.split())
            w = S.fold(''.join(ch for ch in raw.upper() if ch.isalpha()))
            if w and w not in seen:
                seen.add(w)
                out.append((w, raw, group))
    return out


ALL_CRIBS = crib_list()

# ---------------------------------------------------------------- placements
ORDER_PERM = S.ORDER_PERM


def placements(word, seq, alts, use_alts=True):
    """All places a crib fits under SOME one-to-one key: (order, start, slots chosen).
    Uncertain symbols inside the window may take either reading."""
    L = len(word)
    out = []
    for o in ORDERS:
        R = ORDER_PERM[o]
        for i in range(87 - L + 1):
            win = R[i:i + L]
            opts = [[S.slot_of(seq[r])] + ([S.slot_of(alts[r])] if use_alts and r in alts else []) for r in win]
            choice = [None] * L

            def bt(k, m, rv):
                if k == L:
                    out.append((o, i, tuple(choice)))
                    return
                ch = word[k]
                for s in opts[k]:
                    if (ch in m and m[ch] != s) or (s in rv and rv[s] != ch):
                        continue
                    new = ch not in m
                    m[ch], rv[s] = s, ch
                    choice[k] = s
                    bt(k + 1, m, rv)
                    if new:
                        del m[ch], rv[s]
            bt(0, {}, {})
    return out


def eligibility(seq, alts, n=100, seed=4242):
    rng = random.Random(seed)
    copies = []
    for k in range(n):
        s, a = S.null_tokens('doublet', seq, alts, rng.getrandbits(48))
        copies.append((s, a))
    elig, table = [], []
    for w, raw, group in ALL_CRIBS:
        mean_fits = statistics.mean(len(placements(w, s, a)) for s, a in copies)
        ok = len(w) >= ELIG_MIN_LEN and mean_fits <= ELIG_MAX_FITS
        table.append((w, raw, group, mean_fits, ok))
        if ok:
            elig.append(w)
    return elig, table


# ---------------------------------------------------------------- constrained annealer
def anneal_fixed(c, A, fixed, frozen, restarts, iters, seed):
    """Any key with some symbol->letter pairs fixed (the crib). c: slots in reading order
    (window readings already applied); A: alternate slot or -1; frozen: reading positions
    that may not flip."""
    Q = S.QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(c)
    amb = [q for q in range(n) if A[q] >= 0 and q not in frozen]
    fixed_letters = set(fixed.values())
    free_slots = [s for s in range(24) if s not in fixed]
    free_letters = [l for l in range(24) if l not in fixed_letters]
    isfree = [s not in fixed for s in range(24)]
    fac = (S.SA_T1 / S.SA_T0) ** (1.0 / iters)
    allst = range(n - 3)

    def qs(p, ss):
        t = 0.0
        for s in ss:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    def starts_of(lst):
        return tuple(sorted({s for q in lst for s in range(max(0, q - 3), min(q, n - 4) + 1)}))

    best_all = (-1e18, None)
    nf = len(free_slots)
    for _ in range(restarts):
        e = list(c)
        for q in amb:
            if rand() < 0.5:
                e[q] = A[q]
        pos = [[] for _ in range(24)]
        for q, x in enumerate(e):
            pos[x].append(q)
        starts = [starts_of(pos[x]) for x in range(24)]
        pc = {}
        key = [0] * 24
        for s, l in fixed.items():
            key[s] = l
        fl = list(free_letters)
        rng.shuffle(fl)
        for s, l in zip(free_slots, fl):
            key[s] = l
        p = [key[x] for x in e]
        cur = qs(p, allst)
        best, bp = cur, p[:]
        Tm = S.SA_T0
        for _ in range(iters):
            if amb and rand() < S.PFLIP:
                q = amb[rr(len(amb))]
                old = e[q]
                new = c[q] if old == A[q] else A[q]
                ss = range(max(0, q - 3), min(q, n - 4) + 1)
                o = qs(p, ss)
                p[q] = key[new]
                d = qs(p, ss) - o
                if d >= 0 or rand() < exp(d / Tm):
                    e[q] = new
                    pos[old].remove(q)
                    pos[new].append(q)
                    starts[old] = starts_of(pos[old])
                    starts[new] = starts_of(pos[new])
                    pc.clear()
                    cur += d
                    if cur > best:
                        best, bp = cur, p[:]
                else:
                    p[q] = key[old]
                Tm *= fac
                continue
            i = e[rr(n)]
            if not isfree[i]:
                i = free_slots[rr(nf)]
            j = free_slots[rr(nf)]
            if j == i:
                Tm *= fac
                continue
            a, b = (i, j) if i < j else (j, i)
            ss = pc.get(a * 24 + b)
            if ss is None:
                ss = tuple(sorted(set(starts[a]) | set(starts[b])))
                pc[a * 24 + b] = ss
            o = qs(p, ss)
            li, lj = key[i], key[j]
            for q in pos[i]:
                p[q] = lj
            for q in pos[j]:
                p[q] = li
            d = qs(p, ss) - o
            if d >= 0 or rand() < exp(d / Tm):
                key[i], key[j] = lj, li
                cur += d
                if cur > best:
                    best, bp = cur, p[:]
            else:
                for q in pos[i]:
                    p[q] = li
                for q in pos[j]:
                    p[q] = lj
            Tm *= fac
        assert abs(best - qs(bp, allst)) < 1e-6
        if best > best_all[0]:
            best_all = (best, bp)
    return best_all[0] / (n - 3), ''.join(LET[x] for x in best_all[1])


def complete_job(corpus, seq, alts, items, restarts, iters, seed):
    """items: list of (crib, order, start, slots). Returns [(score, text)] per item."""
    S.use_style(corpus, 'S0')
    out = []
    for k, (w, o, i, slots) in enumerate(items):
        R = ORDER_PERM[o]
        c = [S.slot_of(seq[r]) for r in R]
        A = [S.slot_of(alts[r]) if r in alts else -1 for r in R]
        fixed = {}
        frozen = set()
        for k2, s in enumerate(slots):
            q = i + k2
            if c[q] != s:                 # window used the alternate reading here
                A[q], c[q] = c[q], s
            frozen.add(q)
            fixed[s] = LIDX[w[k2]]
        out.append(anneal_fixed(c, A, fixed, frozen, restarts, iters, S.derive_seed(seed, k)))
    return out


# ---------------------------------------------------------------- C2: grouped keys
QT = None


def group_job(corpus, seq, cribs):
    """Best score among the 10,321,920 grouped keys (x2 groupings) that spell a crib."""
    global QT
    if QT is None:
        S.use_style(corpus, 'S0')
        QT = np.array(S.QL, dtype=np.float64)
    best = (-99.0, '', '', '')
    hits = 0
    for gname, sets in AK.GROUPINGS.items():
        sidx = AK.letters_of(sets)
        for o in ORDERS:
            s_ord = [seq[r] for r in ORDER_PERM[o]]
            ci = AK.code_index(s_ord)
            plist = [(w, i) for w in cribs for (oo, i, _) in placements(w, seq, {}, False) if oo == o]
            if not plist:
                continue
            for flip in range(256):
                KL = AK.key_tables(sidx, flip)
                D = None
                for w, i in plist:
                    L = len(w)
                    cols = ci[i:i + L]
                    target = np.array([LIDX[ch] for ch in w])
                    mask = np.all(KL[:, cols] == target, axis=1)
                    if not mask.any():
                        continue
                    if D is None:
                        D = KL[:, ci]
                    rows = np.nonzero(mask)[0]
                    P = D[rows]
                    q = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
                    sc = QT[q].sum(axis=1) / 84
                    hits += len(rows)
                    b = int(sc.argmax())
                    if sc[b] > best[0]:
                        best = (float(sc[b]), w, ''.join(LET[x] for x in P[b]),
                                f"{gname} grouping, {o}, {AK.describe(sets, int(rows[b]), flip)}")
    return best, hits


LGK = None
LG_GROUPS = ['ABCDEFGH', 'IKLMNOPQ', 'RSTUWXYZ']   # I = I/J, U = U/V (same as Dorabella_loopgroups.py)


class LG:
    """The 24,576 loop-group keys: groups to loop counts (6) x each group forwards or
    backwards (8) x each group rotated to any starting direction (512)."""
    _cache = None

    @staticmethod
    def build_keys():
        if LG._cache is None:
            import itertools
            keys, rows = [], []
            for perm in itertools.permutations(range(3)):
                for rev in itertools.product((0, 1), repeat=3):
                    for rots in itertools.product(range(8), repeat=3):
                        row = [0] * 24
                        for loop in range(3):
                            g = list(LG_GROUPS[perm[loop]])
                            g = g[::-1] if rev[loop] else g
                            g = g[rots[loop]:] + g[:rots[loop]]
                            for d in range(8):
                                row[loop * 8 + d] = LIDX[g[d]]
                        rows.append(row)
                        keys.append((perm, rev, rots))
            LG._cache = (keys, np.array(rows, dtype=np.int64))
        return LG._cache


def loopgroup_job(corpus, seq, cribs):
    """Best score among the 24,576 loop-group keys that spell a crib somewhere."""
    global QT, LGK
    if QT is None:
        S.use_style(corpus, 'S0')
        QT = np.array(S.QL, dtype=np.float64)
    if LGK is None:
        LGK = LG.build_keys()
    keys, KL = LGK
    best = (-99.0, '', '', '')
    hits = 0
    for o in ORDERS:
        ci = np.array([S.slot_of(seq[r]) for r in ORDER_PERM[o]])
        D = KL[:, ci]
        for w in cribs:
            for (oo, i, _) in placements(w, seq, {}, False):
                if oo != o:
                    continue
                target = np.array([LIDX[ch] for ch in w])
                mask = np.all(D[:, i:i + len(w)] == target, axis=1)
                if not mask.any():
                    continue
                rows = np.nonzero(mask)[0]
                P = D[rows]
                q = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
                sc = QT[q].sum(axis=1) / 84
                hits += len(rows)
                b = int(sc.argmax())
                if sc[b] > best[0]:
                    perm, rev, rots = keys[int(rows[b])]
                    best = (float(sc[b]), w, ''.join(LET[x] for x in P[b]),
                            f"{o}, loop groups {perm}, reversed {rev}, rotations {rots}")
    return best, hits


def plant_c3(cribs, k, seed):
    rng = random.Random(seed)
    crib = rng.choice([w for w in cribs if len(w) >= 7] or cribs)
    words = S.words_of_text(M.CONTROL_TEXTS[k % len(M.CONTROL_TEXTS)] + ' ' +
                            M.CONTROL_TEXTS[(k + 7) % len(M.CONTROL_TEXTS)])
    body = S.fold(''.join(words))
    at = rng.randrange(0, 87 - len(crib))
    plain = (body[:at] + crib + body[at:])[:87]
    keys, KL = LG.build_keys()
    row = KL[rng.randrange(len(KL))]
    inv = {int(l): s for s, l in enumerate(row)}
    order = rng.choice(ORDERS)
    R = ORDER_PERM[order]
    seq = [0] * 87
    for i, ch in enumerate(plain):
        seq[R[i]] = S.code_of_slot(inv[LIDX[ch]])
    return seq, plain, crib, order


# ---------------------------------------------------------------- one full C1 search
def c1_search(pool, corpus, seq, alts, cribs, P, seed):
    plist = [(w, o, i, sl) for w in cribs for (o, i, sl) in placements(w, seq, alts)]
    if not plist:
        return dict(stat=-99.0, crib='', text='', order='', h1=-99.0, n=0)
    chunk = max(1, math.ceil(len(plist) / (4 * P['workers'])))
    chunks = [plist[k:k + chunk] for k in range(0, len(plist), chunk)]
    r1 = run(pool, complete_job, [(corpus, seq, alts, ch, P['r1'], P['i1'], S.derive_seed(seed, 'a', k))
                                  for k, ch in enumerate(chunks)])
    scored = [(res[0], it) for ch, rr_ in zip(chunks, r1) for it, res in zip(ch, rr_)]
    scored.sort(key=lambda x: -x[0])
    top = [it for _, it in scored[:P['keep']]]
    h1_items = [it for _, it in scored if it[0] in H1_CRIBS][:4]
    deep = top + [it for it in h1_items if it not in top]
    r2 = run(pool, complete_job, [(corpus, seq, alts, [it], P['r2'], P['i2'], S.derive_seed(seed, 'b', k))
                                  for k, it in enumerate(deep)])
    res = [(r[0][0], r[0][1], it) for r, it in zip(r2, deep)]
    best = max(res, key=lambda x: x[0])
    h1 = max([x[0] for x in res if x[2][0] in H1_CRIBS], default=-99.0)
    h1_text = max([x for x in res if x[2][0] in H1_CRIBS], key=lambda x: x[0], default=None)
    return dict(stat=best[0], crib=best[2][0], text=best[1], order=best[2][1], start=best[2][2],
                h1=h1, h1_text=(h1_text[1] if h1_text else ''), n=len(plist),
                top=[dict(score=x[0], crib=x[2][0], order=x[2][1], text=x[1])
                     for x in sorted(res, key=lambda x: -x[0])[:5]])


def run(pool, fn, argsets):
    if pool is None:
        return [fn(*a) for a in argsets]
    res = pool.starmap_async(fn, argsets, chunksize=1)
    while True:
        try:
            return res.get(60)
        except mp.TimeoutError:
            continue


def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    S.corpus_words(corpus)


# ---------------------------------------------------------------- controls
def plant_c1(cribs, k, noise, seed):
    rng = random.Random(seed)
    crib = rng.choice([w for w in cribs if len(w) >= 7] or cribs)
    words = S.words_of_text(M.CONTROL_TEXTS[k % len(M.CONTROL_TEXTS)] + ' ' +
                            M.CONTROL_TEXTS[(k + 7) % len(M.CONTROL_TEXTS)])
    body = S.fold(''.join(words))
    at = rng.randrange(0, 87 - len(crib))
    plain = (body[:at] + crib + body[at:])[:87]
    slots = list(range(24))
    rng.shuffle(slots)
    order = rng.choice(ORDERS)
    R = ORDER_PERM[order]
    seq = [0] * 87
    for i, ch in enumerate(plain):
        seq[R[i]] = S.code_of_slot(slots[LIDX[ch]])
    for q in rng.sample([i for i in range(87) if i not in S.ALTS], noise):
        seq[q] = S.neighbour(seq[q], rng)
    alts = {}
    for q in sorted(S.ALTS):
        wrong = S.neighbour(seq[q], rng)
        if rng.random() < 0.5:
            alts[q] = wrong
        else:
            alts[q], seq[q] = seq[q], wrong
    return seq, alts, plain, crib, order


def plant_c2(cribs, k, seed):
    """Planted English containing a crib, enciphered with a random GROUPED key."""
    rng = random.Random(seed)
    crib = rng.choice([w for w in cribs if len(w) >= 7] or cribs)
    words = S.words_of_text(M.CONTROL_TEXTS[k % len(M.CONTROL_TEXTS)] + ' ' +
                            M.CONTROL_TEXTS[(k + 7) % len(M.CONTROL_TEXTS)])
    body = S.fold(''.join(words))
    at = rng.randrange(0, 87 - len(crib))
    plain = (body[:at] + crib + body[at:])[:87]
    gname = rng.choice(list(AK.GROUPINGS))
    KL = AK.key_tables(AK.letters_of(AK.GROUPINGS[gname]), rng.randrange(256))
    row = KL[rng.randrange(len(KL))]
    inv = {}
    for col, l in enumerate(row):
        inv.setdefault(int(l), []).append(col)
    order = rng.choice(ORDERS)
    R = ORDER_PERM[order]
    seq = [0] * 87
    for i, ch in enumerate(plain):
        l = LIDX[ch]
        if l not in inv:                       # letter missing from this grouping (e.g. Z)
            l = LIDX['S']
        col = rng.choice(inv[l])
        d, loop = divmod(col, 3)
        seq[R[i]] = (loop + 1) * 10 + d + 1
    return seq, plain, crib, order


# ---------------------------------------------------------------- state and batches
def now():
    return datetime.datetime.now()


def stamp():
    return f"[{now():%Y-%m-%d %H:%M:%S}]"


def blank_state(args, elig, table):
    return dict(version=VERSION, created=str(now()), elapsed=0.0, sessions=0,
                seq=S.SEQ, alts={str(k): v for k, v in S.ALTS.items()},
                eligible=elig, elig_table=table,
                params=dict(seed=args.seed, real_seeds=args.real_seeds, r1=args.r1, i1=args.i1,
                            keep=args.keep, r2=args.r2, i2=args.i2, workers=args.workers),
                targets=dict(C1=dict(shuffle=args.copies, doublet=args.copies, controls=args.controls),
                             C2=dict(shuffle=args.copies_c2, doublet=args.copies_c2, controls=args.controls),
                             C3=dict(shuffle=args.copies_c3, doublet=args.copies_c3, controls=args.controls)),
                st=dict(C1=dict(real=[], shuffle=[], doublet=[], controls=[]),
                        C2=dict(real=[], shuffle=[], doublet=[], controls=[]),
                        C3=dict(real=[], shuffle=[], doublet=[], controls=[])),
                log=[])


def stages(state):
    out = []
    for t in ('C1', 'C2', 'C3'):
        st, tg = state['st'][t], state['targets'][t]
        need_real = state['params']['real_seeds'] if t == 'C1' else 1
        if len(st['real']) < need_real:
            out.append((0, 0.0, t, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], t, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            fr = len(st[kind]) / tg[kind]
            if fr < 1:
                out.append((2, fr, t, kind))
    return sorted(out)


def do_unit(pool, state, corpus, test, kind):
    P = state['params']
    st = state['st'][test]
    elig = state['eligible']
    n = len(st[kind])
    seed = S.derive_seed(P['seed'], test, kind, n)
    if test == 'C1':
        if kind == 'real':
            r = c1_search(pool, corpus, S.SEQ, S.ALTS, elig, P, seed)
            st['real'].append(r)
        elif kind in ('shuffle', 'doublet'):
            s, a = S.null_tokens(kind, S.SEQ, S.ALTS, S.derive_seed(P['seed'], 'null', kind, n))
            r = c1_search(pool, corpus, s, a, elig, P, seed)
            st[kind].append(dict(stat=round(r['stat'], 5), h1=round(r['h1'], 5), crib=r['crib'], n=r['n']))
        else:
            noise = 0 if n % 2 == 0 else 2
            s, a, plain, crib, order = plant_c1(elig, n, noise, S.derive_seed(P['seed'], 'plant', n))
            r = c1_search(pool, corpus, s, a, elig, P, seed)
            acc = sum(x == y for x, y in zip(r['text'], plain)) / 87 if r['order'] == order else 0.0
            st['controls'].append(dict(stat=round(r['stat'], 5), noise=noise, planted=crib,
                                       found=r['crib'], acc=round(acc, 3), plain=plain, text=r['text']))
        return
    # C2 / C3: one search per job, a batch of jobs in parallel
    fn = group_job if test == 'C2' else loopgroup_job
    planter = plant_c2 if test == 'C2' else plant_c3
    count = P['workers'] * (1 if test == 'C2' else 4) if kind != 'real' else 1
    count = max(1, min(count, state['targets'][test].get(kind, 1) - n)) if kind != 'real' else 1
    jobs, meta = [], []
    for k in range(count):
        m = n + k
        if kind == 'real':
            s = S.SEQ
            meta.append(None)
        elif kind in ('shuffle', 'doublet'):
            s, _ = S.null_tokens(kind, S.SEQ, S.ALTS, S.derive_seed(P['seed'], 'null' + test, kind, m))
            meta.append(None)
        else:
            s, plain, crib, order = planter(elig, m, S.derive_seed(P['seed'], 'plant' + test, m))
            meta.append((plain, crib, order))
        jobs.append((corpus, s, elig))
    res = run(pool, fn, jobs)
    for (best, hits), mt in zip(res, meta):
        if kind == 'real':
            st['real'].append(dict(stat=best[0], crib=best[1], text=best[2], key=best[3], hits=hits))
        elif kind in ('shuffle', 'doublet'):
            st[kind].append(dict(stat=round(best[0], 5), crib=best[1], hits=hits))
        else:
            plain, crib, order = mt
            acc = sum(x == y for x, y in zip(best[2], plain)) / 87 if order in best[3] else 0.0
            st['controls'].append(dict(stat=round(best[0], 5), planted=crib, found=best[1],
                                       acc=round(acc, 3), plain=plain, text=best[2]))


# ---------------------------------------------------------------- report
def pval(real, nulls):
    if real is None or not nulls:
        return None
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))] if v else None


def f(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def build_report(state, sd):
    L = ["DORABELLA CRIB TEST REPORT",
         f"written {now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"eligible cribs ({len(state['eligible'])}): " + ', '.join(state['eligible']), ""]
    L.append("=" * 104)
    L.append("HEADLINE (p vs PAIR-PRESERVING copies is the one that matters)")
    L.append("=" * 104)
    rows = []
    for t, name in (('C1', 'C1 any key, crib fixed'), ('C2', 'C2 Elgar-style grouped keys'),
                    ('C3', 'C3 loop-group keys')):
        st = state['st'][t]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        cs = [c['stat'] for c in st['controls']]
        q95 = pct(dn, .95)
        power = sum(1 for c in cs if q95 is not None and c > q95) / len(cs) if (cs and q95 is not None) else None
        rows.append((t, name, med, dn, sn, cs, power))
        L.append(f"{name:<32} real {f(med):>8}  pair n={len(dn):>4} p={f(pval(med, dn), 4):>7}  "
                 f"shuffle n={len(sn):>4} p={f(pval(med, sn), 4):>7}  power "
                 f"{'-' if power is None else f'{power:.0%}'}")
    st = state['st']['C1']
    if st['real']:
        h1 = statistics.median(r['h1'] for r in st['real'])
        dh = [x['h1'] for x in st['doublet']]
        fit = sum(1 for x in dh if x > -90)
        L.append(f"{'H1 MY DEAR DORA (found by looking)':<32} real {f(h1):>8}  pair n={len(dh):>4} "
                 f"p={f(pval(h1, dh), 4):>7}  ({fit} of {len(dh)} copies had a DEAR DORA fit at all)")
    L += ["", "Line to beat: p < 0.0167 for C1, C2 or C3 (three tests); H1 is exploratory whatever its p.",
          "A real message would also show as a completed text that READS as English beyond the crib."]
    for t, name, med, dn, sn, cs, power in rows:
        st = state['st'][t]
        L += ["", "-" * 104, name, "-" * 104]
        if not st['real']:
            L.append("real run not done yet")
            continue
        if t == 'C1':
            L.append("real per seed: " + ', '.join(f"{r['stat']:.3f} ({r['crib']})" for r in st['real']))
            b = max(st['real'], key=lambda r: r['stat'])
            L.append(f"placements tried on the real cipher: {b['n']}")
            for x in b.get('top', []):
                L.append(f"  [{x['score']:.3f}] crib {x['crib']} in {x['order']}")
                for a in range(0, 87, 29):
                    L.append("      " + x['text'][a:a + 29])
                L.append("      words: " + sd.segment(x['text']))
            if b.get('h1_text'):
                L.append(f"  H1 best completion ({b['h1']:.3f}):")
                for a in range(0, 87, 29):
                    L.append("      " + b['h1_text'][a:a + 29])
                L.append("      words: " + sd.segment(b['h1_text']))
        else:
            r = st['real'][0]
            L.append(f"real: best {r['stat']:.3f} with crib {r['crib'] or '(no crib fits any key of this kind)'}; "
                     f"crib-consistent keys found: {r['hits']:,}")
            if r['text']:
                L.append(f"  key: {r['key']}")
                for a in range(0, 87, 29):
                    L.append("      " + r['text'][a:a + 29])
                L.append("      words: " + sd.segment(r['text']))
            dh = [x['hits'] for x in st['doublet']]
            if dh:
                L.append(f"crib-consistent keys in pair-preserving copies: mean {statistics.mean(dh):,.0f}, "
                         f"p (real has as many or more) = {pval(r['hits'], dh):.4f}")
        for kind, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L.append(f"  vs {kind:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{pct(nul, .95):.3f}, max {max(nul):.3f}")
        if cs:
            c = st['controls']
            L.append(f"controls: n={len(c)}, median stat {statistics.median(cs):.3f}, planted crib found "
                     f"{sum(1 for x in c if x['found'] == x['planted'])}/{len(c)}, median letters right "
                     f"{statistics.median(x['acc'] for x in c):.0%}")
            ex = max(c, key=lambda x: x['acc'])
            L.append(f"  example: planted {ex['plain'][:50]}")
            L.append(f"           found   {ex['text'][:50]}")
    L += ["", "-" * 104, "ELIGIBILITY (decided on 100 pair-preserving copies, never on the real cipher)", "-" * 104]
    for w, raw, group, mf, ok in state['elig_table']:
        L.append(f"  {'USED   ' if ok else 'skipped'} {w:<28} {group:<26} fits per copy {mf:6.1f}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='cribs_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--real-seeds', type=int, default=5)
    ap.add_argument('--r1', type=int, default=2, help='first-pass restarts per placement')
    ap.add_argument('--i1', type=int, default=6000, help='first-pass steps per restart')
    ap.add_argument('--keep', type=int, default=12, help='placements searched again deeply')
    ap.add_argument('--r2', type=int, default=10)
    ap.add_argument('--i2', type=int, default=20000)
    ap.add_argument('--copies', type=int, default=1000)
    ap.add_argument('--copies-c2', type=int, default=500)
    ap.add_argument('--copies-c3', type=int, default=1000)
    ap.add_argument('--controls', type=int, default=60)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    S.SEQ, S.ALTS = S.load_cipher()
    print(f"cipher read from Dorabella_transcription.csv ({len(S.ALTS)} uncertain symbols)")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))
    sd = S.StyleDict(S.base_dictionary(args.corpus), 'S0')

    def log(msg):
        line = f"{stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != S.SEQ:
            raise SystemExit("saved state was made with a different transcription; use --fresh")
    else:
        print("deciding which cribs are eligible (on scrambled copies only)...", flush=True)
        elig, table = eligibility(S.SEQ, S.ALTS)
        state = blank_state(args, elig, table)
    if args.report_only:
        with open(rp, 'w', encoding='utf-8') as fh:
            fh.write(build_report(state, sd))
        return
    state['params']['workers'] = args.workers
    S.keep_awake()
    pool = mp.Pool(args.workers, initializer=init_worker, initargs=(args.corpus,)) if args.workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, state, sd)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {args.workers}, "
            f"{len(state['eligible'])} eligible cribs")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stg = stages(state)
            if not stg:
                for t in ('C1', 'C2', 'C3'):
                    for kind in ('shuffle', 'doublet'):
                        state['targets'][t][kind] = int(state['targets'][t][kind] * 1.5)
                log("all targets met; raising copy targets by 50% and carrying on")
                continue
            _, _, test, kind = stg[0]
            t0 = time.time()
            do_unit(pool, state, args.corpus, test, kind)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['st'][test]
            msg = f"{test} {kind} in {H.fmt_secs(t1 - t0)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][test][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{stamp()} {msg}", flush=True)
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, sd))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished unit is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, sd))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


def selftest(args, pool, state, sd):
    print(f"\n{stamp()} SELF-TEST")
    print(f"{len(ALL_CRIBS)} cribs in the list; eligible after the rule: {len(state['eligible'])}")
    for w, raw, group, mf, ok in state['elig_table']:
        print(f"  {'USED   ' if ok else 'skipped'} {w:<28} fits per copy {mf:6.1f}")
    P = state['params']
    print("\nC1 control (planted English containing a crib, 2 misreads), then one pair-preserving copy:")
    s, a, plain, crib, order = plant_c1(state['eligible'], 3, 2, 99)
    t0 = time.time()
    r = c1_search(pool, args.corpus, s, a, state['eligible'], P, 7)
    dt = time.time() - t0
    acc = sum(x == y for x, y in zip(r['text'], plain)) / 87 if r['order'] == order else 0.0
    print(f"  planted crib {crib}; found crib {r['crib']}; score {r['stat']:.3f}; letters right {acc:.0%}; "
          f"{r['n']} placements; {H.fmt_secs(dt)}")
    print(f"  planted {plain[:60]}\n  found   {r['text'][:60]}")
    ns, na = S.null_tokens('doublet', S.SEQ, S.ALTS, 12345)
    t0 = time.time()
    rn = c1_search(pool, args.corpus, ns, na, state['eligible'], P, 8)
    dn = time.time() - t0
    print(f"  copy: best {rn['stat']:.3f} with {rn['crib']}; {rn['n']} placements; {H.fmt_secs(dn)}")
    print("\nC2 control (planted English with a crib, random grouped key), then one copy:")
    s2, plain2, crib2, order2 = plant_c2(state['eligible'], 3, 99)
    t0 = time.time()
    (b2, h2), = run(pool, group_job, [(args.corpus, s2, state['eligible'])])
    d2 = time.time() - t0
    print(f"  planted {crib2}; found {b2[1]}; score {b2[0]:.3f}; {h2} crib-consistent keys; {H.fmt_secs(d2)}")
    print(f"  planted {plain2[:60]}\n  found   {b2[2][:60]}")
    (b3, h3), = run(pool, group_job, [(args.corpus, ns, state['eligible'])])
    print(f"  copy: best {b3[0]:.3f} with {b3[1] or '-'}; {h3} crib-consistent keys")
    print("\nC3 control (planted English with a crib, random loop-group key), then one copy:")
    s3, plain3, crib3, order3 = plant_c3(state['eligible'], 3, 99)
    t0 = time.time()
    (b4, h4), = run(pool, loopgroup_job, [(args.corpus, s3, state['eligible'])])
    d3 = time.time() - t0
    print(f"  planted {crib3}; found {b4[1]}; score {b4[0]:.3f}; {h4} crib-consistent keys; {H.fmt_secs(d3)}")
    print(f"  planted {plain3[:60]}\n  found   {b4[2][:60]}")
    (b5, h5), = run(pool, loopgroup_job, [(args.corpus, ns, state['eligible'])])
    print(f"  copy: best {b5[0]:.3f} with {b5[1] or '-'}; {h5} crib-consistent keys")
    per_c1 = (dt + dn) / 2
    tot = (2 * args.copies + args.controls + args.real_seeds) * per_c1 + \
          (2 * args.copies_c2 + args.controls) * d2 / max(1, args.workers) + \
          (2 * args.copies_c3 + args.controls) * d3 / max(1, args.workers)
    print(f"\nEstimated time to reach the targets: about {H.fmt_secs(tot)}")


if __name__ == '__main__':
    mp.freeze_support()
    main()
