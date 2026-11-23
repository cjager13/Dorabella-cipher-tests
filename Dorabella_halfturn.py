#!/usr/bin/env python3
"""
Dorabella HALF-TURN test (roadmap items H4 and H1): what are the half-turned neighbours, and
does turning every second symbol explain them?

Needs Dorabella_notebook.py, Dorabella_three.py, Dorabella_language.py and everything they
need (Dorabella_hillclimb.py, Dorabella_tests.py, Dorabella_marathon.py,
Dorabella_transcription.csv, english.txt), plus numpy.

  python Dorabella_halfturn.py --selftest          # checks, one control per family, time estimate
  python Dorabella_halfturn.py --stop-when-done    # results in halfturn_out\\
  python Dorabella_halfturn.py --report-only

Ctrl+C any time; the same command resumes.

THE HABIT (found in the November 2026 audit)
Neighbouring symbols are the same shape turned half round (opposite directions) 27 times in
the Dorabella, against about 15 in English through Elgar's notebook key; 13 of them also share
a loop count, against about 6. Under the notebook key a half-turned symbol is always the letter
12 places along the 24-letter alphabet (A-N, E-R, H-U...). Pair-preserving copies keep the
half-turn count exactly, so no earlier test could see it.

PRE-REGISTERED (frozen 18 Nov 2026, before any run; three tests, so the line is p < 0.0167)
  H4  THE TURN PROFILE. For the 86 neighbour pairs, count the turn between their openings:
      none, 45, 90, 135 or 180 degrees (5 counts, written order, across line ends). Statistic:
      the Mahalanobis distance of the Dorabella's 5 counts from 200,000 English passages of
      87 letters through the notebook key; p = share of English passages at least as far out.
      Each count is also reported on its own, two-sided (description only), with the letter
      pairs that carry the half-turns. HONESTY NOTE: the counts were looked at on 18 Nov 2026
      while planning, so H4 states the habit exactly; it is not an independent discovery.
  H1A ALTERNATE SYMBOLS TURNED, notebook key, no search. Every second symbol is changed by one
      of the 48 symbol transforms (8 turns x mirror or not x loop count stepped 0/1/2), the
      others by another (48 x 48 = 2,304 settings), in each of the 5 reading orders, with the
      alternation running through the text or restarting at each line (10 configurations):
      23,040 readings. Statistic: the best quadgram score.
  H1B ALTERNATE SYMBOLS TURNED, any key. The same alternation with ANY one-symbol-one-letter
      key (the first symbol of each pair is left as it is, since the key absorbs it; 48 transforms
      for the second, the identity included). For each configuration the transforms are ranked
      without a key (the right one makes the symbols look like one simple substitution: on
      planted English it ranked first 34 times in 40, and never below 8th) and the top --keep
      are searched by annealing. Statistic: best quadgram score; real = median over 5 seeds.
  H1A and H1B: compared with pair-preserving copies (main) and plain shuffles put through the
  identical search; power from planted English written that way (0 or 2 misread symbols).
  A pass also needs a decrypt that reads as English.
  Transcription: the consensus reading (Hauer et al. 2025; symbols 23, 32, 34, 85, 86 read 14,
  11, 12, 14, 12), as in Dorabella_notebook.py. H4 also reports the range over all 1,024
  readings of the 10 uncertain symbols.
"""
import argparse
import itertools
import json
import os
import random
import re
import signal
import statistics
import time

import numpy as np

import Dorabella_notebook as NB          # applies the notebook key and the consensus reading
import Dorabella_language as LANG

Z = NB.Z
H, M = Z.H, Z.M
VERSION = 1
FAMILIES = ['H1A', 'H1B']
FAMILY_DESC = {'H1A': "notebook key, every second symbol transformed (2,304 settings x 10 configurations), no search",
               'H1B': "any key, every second symbol transformed (48 transforms x 10 configurations), screened then annealed"}
BONF = 0.05 / 3
DEFAULT_TARGETS = {'H1A': dict(shuffle=2000, doublet=5000, controls=100),
                   'H1B': dict(shuffle=300, doublet=600, controls=40)}
TURN_NAMES = ['none', '45', '90', '135', '180']
SEQ = None
CONFIGS = None
G2 = np.array(list(itertools.product(range(48), repeat=2)), dtype=np.int64)


# ============================================================== configurations
def configs():
    """(label, reading order R, phase 0/1 per reading position) x 10."""
    global CONFIGS
    if CONFIGS is None:
        line_of = [0] * 29 + [1] * 31 + [2] * 27
        out = []
        for o in Z.KNOWN_ORDERS:
            R = Z.known_order(o)
            out.append((f"{o}, alternation runs through", R, [i % 2 for i in range(87)]))
            cnt, ph = [0, 0, 0], []
            for p in R:
                ph.append(cnt[line_of[p]] % 2)
                cnt[line_of[p]] += 1
            out.append((f"{o}, alternation restarts each line", R, ph))
        CONFIGS = out
    return CONFIGS


def tlabel(g):
    return Z.TURN_LABEL[g]


# ============================================================== H4: the turn profile
def turn_counts(seq):
    c = [0] * 5
    for a, b in zip(seq, seq[1:]):
        t = (a % 10 - b % 10) % 8
        c[min(t, 8 - t)] += 1
    return c


def same_loop_halfturns(seq):
    return sum(1 for a, b in zip(seq, seq[1:]) if (a % 10 - b % 10) % 8 == 4 and a // 10 == b // 10)


def english_profiles(corpus, n, seed):
    """Turn counts and same-loop half-turns of n English passages (87 letters) through the notebook key."""
    letters = LANG.load_letters(corpus).replace('J', 'I').replace('V', 'U')
    code = {}
    for s, L in enumerate(Z.ELGAR):
        code[H.LETTERS[L]] = Z.code_of_slot(s)
    arr = np.array([code[ch] for ch in letters], dtype=np.int64)
    d = arr % 10
    loop = arr // 10
    t = (d[:-1] - d[1:]) % 8
    t = np.minimum(t, 8 - t)                       # turn between letter i and i+1
    same = ((t == 4) & (loop[:-1] == loop[1:])).astype(np.int64)
    rng = np.random.default_rng(seed)
    starts = rng.integers(0, len(arr) - 87, size=n)
    idx = starts[:, None] + np.arange(86)[None, :]
    T = t[idx]
    counts = np.stack([(T == k).sum(axis=1) for k in range(5)], axis=1)
    return counts, same[idx].sum(axis=1), letters


def h4(corpus, n, seed):
    counts, same, letters = english_profiles(corpus, n, seed)
    seq, alts = LANG.load_cipher()
    alts = dict(alts)
    alts.setdefault(31, 11)                      # symbol 32: 18 or 11
    readings = np.array([turn_counts(s) for s in LANG.all_readings(seq, alts)])
    cons = turn_counts(SEQ)
    dx, p, d_ref = LANG.mahalanobis_p(cons, counts)
    per = [LANG.mahalanobis_p(r, counts)[1] for r in readings[::16]]
    per_bin = []
    for k in range(5):
        pv, side = LANG.two_sided(cons[k], counts[:, k])
        per_bin.append(dict(turn=TURN_NAMES[k], dora=cons[k], lo=int(readings[:, k].min()), hi=int(readings[:, k].max()),
                            mean=round(float(counts[:, k].mean()), 2),
                            q025=float(np.percentile(counts[:, k], 2.5)), q975=float(np.percentile(counts[:, k], 97.5)),
                            p=round(float(pv), 5), side=side))
    # letter pairs carrying the half-turns, against how often English writes each pair
    txt = Z.text_of([Z.ELGAR[Z.slot_of(c)] for c in SEQ])
    bg = {}
    for i in range(len(letters) - 1):
        k = letters[i:i + 2]
        bg[k] = bg.get(k, 0) + 1
    tot = sum(bg.values())
    pairs = []
    for i in range(86):
        a, b = SEQ[i], SEQ[i + 1]
        if (a % 10 - b % 10) % 8 == 4:
            k = txt[i:i + 2]
            pairs.append(dict(pos=i + 1, letters=k, same_loop=a // 10 == b // 10,
                              per_1000=round(1000 * bg.get(k, 0) / tot, 2)))
    s_d = same_loop_halfturns(SEQ)
    return dict(n=n, dist=round(dx, 2), p=float(p), p_lo=float(min(per)), p_hi=float(max(per)),
                d_med=round(float(np.median(d_ref)), 2), d_99=round(float(np.percentile(d_ref, 99)), 2),
                bins=per_bin, pairs=pairs, same=s_d, same_mean=round(float(same.mean()), 2),
                same_p=float((1 + np.sum(same >= s_d)) / (1 + len(same))))


# ============================================================== H1A: notebook key, no search
def h1a_job(seq, top=1):
    """Best (score, config, g0, g1) per configuration; returns the top entries overall."""
    E = np.array(Z.ELGAR, dtype=np.int64)
    LT = E[np.array(Z.TURN, dtype=np.int64)]                       # LT[g, slot] = letter
    res = []
    for ci, (_, R, ph) in enumerate(configs()):
        c = np.array([Z.slot_of(seq[r]) for r in R])
        ph = np.array(ph)
        P = np.empty((len(G2), 87), dtype=np.int64)
        cols = [np.where(ph == k)[0] for k in range(2)]
        P[:, cols[0]] = LT[G2[:, 0]][:, c[cols[0]]]
        P[:, cols[1]] = LT[G2[:, 1]][:, c[cols[1]]]
        sc = Z.np_scores(P)
        for b in np.argsort(-sc)[:top]:
            res.append((float(sc[b]), ci, int(G2[b, 0]), int(G2[b, 1]), Z.text_of(P[b].tolist())))
    res.sort(key=lambda x: -x[0])
    return res[:max(top, 1) * 3] if top > 1 else res[:1]


# ============================================================== H1B: any key, screened
def screen(c, ph):
    """Rank the 48 transforms of every second symbol without a key."""
    n, out = len(c), []
    for g1 in range(48):
        e = [Z.TURN[g1][c[q]] if ph[q] else c[q] for q in range(n)]
        cnt = [0] * 24
        for x in e:
            cnt[x] += 1
        ic = sum(v * (v - 1) for v in cnt) / (n * (n - 1))
        s2, s3, rep = set(), set(), 0
        for q in range(n - 1):
            t2 = e[q] * 24 + e[q + 1]
            rep += t2 in s2
            s2.add(t2)
            if q < n - 2:
                t3 = t2 * 24 + e[q + 2]
                rep += 2 * (t3 in s3)
                s3.add(t3)
        out.append((ic * 100 + rep * 0.3, g1))
    out.sort(reverse=True)
    return [g for _, g in out]


def h1b_job(seq, keep, restarts, iters, seed):
    best = (-1e18, None, None, None)
    for ci, (_, R, ph) in enumerate(configs()):
        c = [Z.slot_of(seq[r]) for r in R]
        for g1 in screen(c, ph)[:keep]:
            e = [Z.TURN[g1][c[q]] if ph[q] else c[q] for q in range(87)]
            s, p, _ = Z.anneal('masc', e, [0] * 87, restarts, iters, Z.derive_seed(seed, ci, g1))
            if s > best[0]:
                best = (s, ci, g1, Z.text_of(p))
    return best


# ============================================================== controls
def plant(fam, idx, noise, seed):
    rng = random.Random(seed)
    texts = M.CONTROL_TEXTS
    raw = texts[idx % len(texts)] + ' ' + texts[(idx + 7) % len(texts)]
    plain = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))[:87]
    L = [H.LIDX[ch] for ch in plain]
    ci = rng.randrange(10)
    _, R, ph = configs()[ci]
    if fam == 'H1A':
        g = (rng.randrange(48), rng.randrange(48))
        slots = [Z.TURN_INV[g[ph[i]]][Z.ELGAR_INV[x][0]] for i, x in enumerate(L)]
    else:
        g = (0, rng.randrange(1, 48))
        key = list(range(24))
        rng.shuffle(key)
        kinv = Z.inverse(key)
        slots = [Z.TURN_INV[g[ph[i]]][kinv[x]] for i, x in enumerate(L)]
    seq = [0] * 87
    for i, r in enumerate(R):
        seq[r] = Z.code_of_slot(slots[i])
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, ci, g


# ============================================================== workers and batches
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    Z.init_tables(corpus)


def null_seq(kind, fam, n, seed0):
    rng = random.Random(Z.derive_seed(seed0, 'null', fam, kind, n))
    if kind == 'shuffle':
        s = list(SEQ)
        rng.shuffle(s)
        return s
    return M.doublet_shuffle(list(SEQ), rng)


def do_batch(pool, state, fam, kind, count, args):
    st, seed0 = state['fam'][fam], state['params']['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append(SEQ)
            meta.append(n)
        elif kind == 'controls':
            noise = 0 if n % 2 == 0 else 2
            s, plain, ci, g = plant(fam, n, noise, Z.derive_seed(seed0, 'plant', fam, n))
            items.append(s)
            meta.append((noise, plain, ci, g))
        else:
            items.append(null_seq(kind, fam, n, seed0))
            meta.append(None)
    if fam == 'H1A':
        res = Z.run_jobs(pool, h1a_job, [(s, 5 if kind == 'real' else 1) for s in items])
    else:
        P = state['params']
        res = Z.run_jobs(pool, h1b_job, [(s, P['keep'], P['restarts'], P['iters'],
                                          Z.derive_seed(seed0, fam, kind, len(st[kind]) + k))
                                         for k, s in enumerate(items)])
    for r, mt in zip(res, meta):
        if fam == 'H1A':
            top = r
            s, ci, g0, g1, text = top[0]
            g = (g0, g1)
        else:
            s, ci, g1, text = r
            g, top = (0, g1), None
        rec = dict(stat=round(s, 5), config=ci, g=list(g))
        if kind == 'real':
            rec['text'] = text
            if top:
                rec['top'] = [dict(stat=round(x[0], 5), config=x[1], g=[x[2], x[3]], text=x[4]) for x in top]
        elif kind == 'controls':
            noise, plain, tci, tg = mt
            rec.update(noise=noise, plain=plain, text=text,
                       acc=round(sum(a == b for a, b in zip(text, plain)) / 87, 3),
                       found=int(ci == tci and tuple(g) == tuple(tg)))
        st[kind].append(rec)


def stages(state):
    out = []
    for f in FAMILIES:
        st, tg = state['fam'][f], state['targets'][f]
        need = state['params']['real_seeds'] if f == 'H1B' else 1
        if len(st['real']) < need:
            out.append((0, len(st['real']) / need, f, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], f, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            fr = len(st[kind]) / tg[kind]
            if fr < 1:
                out.append((2, fr, f, kind))
    return sorted(out)


def batch_size(state, f, kind, workers):
    if kind == 'real':
        return 1
    n = workers * (8 if f == 'H1A' else 1)
    return max(1, min(n, state['targets'][f][kind] - len(state['fam'][f][kind])))


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def cfg_label(ci, g):
    return f"{configs()[ci][0]} | first of each pair: {tlabel(g[0])}, second: {tlabel(g[1])}"


def build_report(state):
    L = ["DORABELLA HALF-TURN REPORT  (roadmap H4 and H1)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   sessions {state['sessions']}",
         NB.KEY_NOTE, NB.READ_NOTE,
         f"three pre-registered tests (H4, H1A, H1B): the line is p < {BONF:.4f}", ""]
    h = state.get('h4')
    L += ["=" * 104, "H4  THE TURN PROFILE  (86 neighbour pairs; English = 87-letter passages through the notebook key)", "=" * 104]
    if h:
        L.append(f"distance from English: {h['dist']:.1f} (typical English {h['d_med']:.1f}, 99th pct {h['d_99']:.1f});  "
                 f"p = {h['p']:.5f}  over {h['n']:,} passages   (range over the 1,024 readings {h['p_lo']:.5f} to {h['p_hi']:.5f})")
        L.append(f"verdict: {'PAST the line: the turn profile is not that of enciphered English' if h['p'] < BONF else 'not past the line'}")
        L.append("")
        L.append(f"{'turn (degrees)':<16} {'Dorabella':>9} {'readings':>9} {'English mean':>13} {'English 95% range':>18} {'p (2-sided)':>12}")
        for b in h['bins']:
            L.append(f"{b['turn']:<16} {b['dora']:>9} {str(b['lo']) + '-' + str(b['hi']):>9} {b['mean']:>13.2f} "
                     f"{str(int(b['q025'])) + '-' + str(int(b['q975'])):>18} {b['p']:>12.5f}  (Dorabella {b['side']})")
        L.append(f"half-turns that also share a loop count: {h['same']} (English mean {h['same_mean']}, p = {h['same_p']:.5f})")
        L.append("")
        L.append("the half-turned neighbours, read with the notebook key (per 1000 = how often English writes that pair):")
        for kind in (True, False):
            ps = [x for x in h['pairs'] if x['same_loop'] == kind]
            L.append(f"  {'same loops' if kind else 'different loops'} ({len(ps)}): " +
                     ', '.join(f"{x['letters']}@{x['pos']} ({x['per_1000']})" for x in ps))
        L.append("  (under the notebook key a half-turned symbol is the letter 12 places along: A-N, B-O, C-P, D-Q, E-R, "
                 "F-S, G-T, H-U, I-W, K-X, L-Y, M-Z)")
        L.append("  HONESTY NOTE: these counts were looked at while planning (18 Nov 2026); H4 states the habit exactly")
        L.append("  and replicates Hauer et al. (2021); it is not an independent discovery.")
    L.append("")
    for f in FAMILIES:
        st = state['fam'][f]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        q95 = Z.pct(dn, .95) if dn else None
        cs = st['controls']
        power = sum(1 for c in cs if q95 is not None and c['stat'] > q95) / len(cs) if (cs and q95 is not None) else None
        p_d = Z.pval(med, dn) if dn else None
        L += ["=" * 104, f"{f}  {FAMILY_DESC[f]}", "=" * 104]
        if med is None:
            L.append("real run not done yet")
            continue
        verdict = ('waiting for copies' if p_d is None else
                   'BEATS the copies past the line: read the decrypt by eye' if p_d < BONF else
                   'chance level, LOW POWER: inconclusive' if (power is not None and power < 0.5) else 'chance level')
        L.append(f"real {med:.3f}" + (f" (seeds: {', '.join(f'{x:.3f}' for x in real)})" if len(real) > 1 else '') +
                 f"   p vs pair-preserving {fmt(p_d, 4)} (n={len(dn)})   p vs plain shuffle "
                 f"{fmt(Z.pval(med, sn) if sn else None, 4)} (n={len(sn)})   power {'-' if power is None else f'{power:.0%}'}")
        L.append(f"verdict: {verdict}")
        for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L.append(f"  vs {name:<16} mean {statistics.mean(nul):.3f}, 95th pct {Z.pct(nul, .95):.3f}, max {max(nul):.3f}")
        b = max(st['real'], key=lambda r: r['stat'])
        shown = b.get('top') or [dict(stat=b['stat'], config=b['config'], g=b['g'], text=b['text'])]
        L.append("best real decrypts:")
        for x in shown[:5]:
            L.append(f"  [{x['stat']:.3f}] {cfg_label(x['config'], x['g'])}")
            for a in range(0, 87, 29):
                L.append("      " + x['text'][a:a + 29])
            L.append("      words: " + Z.show_words(x['text']))
        if dn:
            from collections import Counter
            wins = Counter(tuple(x['g']) for x in st['doublet']).most_common(3)
            L.append("transforms most often winning in the copies: " +
                     '; '.join(f"{tlabel(g[0])} / {tlabel(g[1])} ({k}x)" for g, k in wins))
        for nz in (0, 2):
            c2 = [c for c in cs if c['noise'] == nz]
            if c2:
                L.append(f"controls, {nz} misread symbols: n={len(c2)}, setting found {sum(c['found'] for c in c2)}/{len(c2)}, "
                         f"median letters right {statistics.median(c['acc'] for c in c2):.0%}, above null 95th pct "
                         f"{sum(1 for c in c2 if q95 is not None and c['stat'] > q95)}/{len(c2)}")
        if cs:
            ex = max(cs, key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:56]}")
            L.append(f"                                   found   {ex['text'][:56]}")
        L.append("")
    L += ["A family matters only if the real cipher beats the PAIR-PRESERVING copies past the line, its controls",
          "show it can find a message, and the decrypt reads as English.", "",
          "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    txt = Z.text_of([Z.ELGAR[Z.slot_of(c)] for c in SEQ])
    print(f"notebook key on the consensus reading: {txt[:40]}...  (widely printed: BPECAHTCKYFRQDRIRRHPPRDXYXGFS)")
    assert txt.startswith('BPECAHTCKYFRQDRIRRHPPRDXYXGFS')
    t0 = time.time()
    h = h4(args.corpus, 20000, args.seed)
    print(f"H4 on 20,000 English passages ({H.fmt_secs(time.time() - t0)}): turn counts {[b['dora'] for b in h['bins']]}, "
          f"English means {[b['mean'] for b in h['bins']]}, p = {h['p']:.5f}")
    for f in FAMILIES:
        seq, plain, ci, g = plant(f, 3, 0, Z.derive_seed('selftest', f))
        t0 = time.time()
        if f == 'H1A':
            s, rci, g0, g1, text = Z.run_jobs(pool, h1a_job, [(seq, 1)])[0][0]
            rg = (g0, g1)
        else:
            s, rci, g1, text = Z.run_jobs(pool, h1b_job, [(seq, args.keep, args.restarts, args.iters, 7)])[0]
            rg = (0, g1)
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(text, plain)) / 87
        nul = null_seq('doublet', f, 0, 12345)
        if f == 'H1A':
            ns = Z.run_jobs(pool, h1a_job, [(nul, 1)])[0][0][0]
        else:
            ns = Z.run_jobs(pool, h1b_job, [(nul, args.keep, args.restarts, args.iters, 8)])[0][0]
        print(f"\n{f}: planted control {s:.3f} (a pair-preserving copy: {ns:.3f}); setting found: "
              f"{'yes' if (rci == ci and tuple(rg) == tuple(g)) else 'no'}; letters right {acc:.0%}; {H.fmt_secs(dt)} per search")
        print(f"  true:  {cfg_label(ci, g)}")
        print(f"  found: {cfg_label(rci, rg)}")
        print(f"  planted {plain[:58]}\n  decoded {text[:58]}")
        tg = DEFAULT_TARGETS[f]
        n = sum(tg.values()) + (5 if f == 'H1B' else 1)
        par = workers if f == 'H1B' else workers * 8
        print(f"  estimated time for {n} searches, {workers} at a time: about {H.fmt_secs(n * dt / max(1, min(workers, par)))}")


# ============================================================== main
def main():
    global SEQ
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='halfturn_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--samples', type=int, default=200000, help='English passages for H4')
    ap.add_argument('--keep', type=int, default=4, help='H1B: transforms searched per configuration after the screen')
    ap.add_argument('--restarts', type=int, default=4)
    ap.add_argument('--iters', type=int, default=10000)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true')
    args = ap.parse_args()
    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    SEQ, source = Z.load_seq()
    print(f"cipher read from: {source}")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))

    def log(msg):
        line = f"{Z.stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    Z.init_tables(args.corpus)
    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != SEQ:
            raise SystemExit("saved state was made with a different transcription; use --fresh or another --out")
    else:
        state = dict(version=VERSION, created=str(Z.now()), elapsed=0.0, sessions=0, seq=SEQ,
                     params=dict(seed=args.seed, keep=args.keep, restarts=args.restarts, iters=args.iters, real_seeds=5),
                     targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                     fam={f: dict(real=[], doublet=[], shuffle=[], controls=[]) for f in FAMILIES}, log=[])
    if args.report_only:
        Z.write_text(rp, build_report(state))
        print(f"report written to {rp}")
        return
    workers = max(1, args.workers)
    Z.keep_awake()
    pool = Z.mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}")
        if not state.get('h4'):
            t0 = time.time()
            state['h4'] = h4(args.corpus, args.samples, args.seed)
            log(f"H4 done in {H.fmt_secs(time.time() - t0)}: p = {state['h4']['p']:.5f}")
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stg = stages(state)
            if not stg:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for f in FAMILIES:
                    for kind in ('doublet', 'shuffle'):
                        state['targets'][f][kind] = int(state['targets'][f][kind] * 1.5)
                log("all targets met; raising the copy targets by 50% and carrying on")
                continue
            _, _, f, kind = stg[0]
            n = batch_size(state, f, kind, workers)
            t0 = time.time()
            do_batch(pool, state, f, kind, n, args)
            now_t = time.time()
            state['elapsed'] += now_t - last
            last = now_t
            st = state['fam'][f]
            msg = f"{f} {kind} x{n} in {H.fmt_secs(now_t - t0)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][f][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{Z.stamp()} {msg}", flush=True)
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    Z.mp.freeze_support()
    main()
