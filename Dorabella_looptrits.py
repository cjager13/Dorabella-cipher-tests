#!/usr/bin/env python3
"""
Dorabella loop triples: do the LOOP COUNTS alone, read three at a time, spell a message?

Needs Dorabella_hillclimb.py, Dorabella_marathon.py, Dorabella_transcription.csv and
english.txt in the same folder.

  python Dorabella_looptrits.py --selftest       # a minute or two: checks + power preview
  python Dorabella_looptrits.py                  # full run (2 workers, so it can share the
                                                 # laptop with Dorabella_three.py)
  python Dorabella_looptrits.py --workers 8      # if the laptop is free

Output: looptrits_out/report.txt and looptrits_out/log.txt

THE IDEA
Each symbol has 1, 2 or 3 loops. Three loop counts make a base-3 number with 27
values: 26 letters plus a space. The cipher has 87 = 3 x 29 symbols, so the loop
counts alone could carry a 29-character message, with the directions as cover.

PRE-REGISTERED TESTS (frozen before any run; 3 tests, so Bonferroni p < 0.0167)

  L1  Letter-likeness, no key needed. A 29-letter message written as triples
      repeats some triples a lot (E, T, A...). Statistic: index of coincidence of
      the 29 triples (forward reading, triples starting at symbol 1).
  L2  Counting keys. Triple = base-3 number (first or last digit most significant),
      read as A..Z plus space, with every rotation of that 27-character alphabet,
      over the 5 reading orders: 270 readings, no search. Statistic: best quadgram
      score (spaces dropped).
  L3  Any key. Any one-triple-one-letter key, found by annealing, over the 5 reading
      orders. Statistic: best quadgram score (real = median of 5 seeds). With only
      29 letters this is expected to have low power; the controls will say.

Each is compared with the IDENTICAL test on scrambled loop sequences: plain shuffles
of the 87 loop counts, and pair-preserving shuffles that keep every adjacent
loop-count pair. Planted English (29 characters, written as triples with a random
key or a random counting key) goes through the same tests to measure power.
p = (1+k)/(n+1).
"""
import argparse
import csv
import datetime
import math
import multiprocessing as mp
import os
import random
import signal
import statistics
import time
from collections import Counter

import Dorabella_hillclimb as H
import Dorabella_marathon as M
import Dorabella_tests as T

ALPHA27 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ_'
ORDERS = ['forward', 'fully reversed', 'each line reversed', 'boustrophedon A', 'boustrophedon B']
TESTS = ['L1', 'L2', 'L3']
BONF = 0.05 / len(TESTS)
QL = None


# ============================================================== data
def load_loops(path='Dorabella_transcription.csv'):
    with open(path, newline='', encoding='utf-8') as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r['position']))
    loops = [int(r['loops']) for r in rows]
    assert len(loops) == 87
    return loops


def order_loops(loops, order):
    return T.apply_order(list(loops), order)


def triples(loops):
    return [(loops[i] - 1) * 9 + (loops[i + 1] - 1) * 3 + (loops[i + 2] - 1)
            for i in range(0, len(loops) - 2, 3)]


def triples_lsb(loops):
    return [(loops[i + 2] - 1) * 9 + (loops[i + 1] - 1) * 3 + (loops[i] - 1)
            for i in range(0, len(loops) - 2, 3)]


def trits_of(v):
    return [v // 9 + 1, (v // 3) % 3 + 1, v % 3 + 1]


# ============================================================== scoring
def init_tables(corpus):
    global QL
    H.init_worker(corpus, None)
    QL = [H.QFLOOR] * (24 ** 4)
    for k, v in H.QUAD.items():
        QL[k] = v


def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    init_tables(corpus)


def qscore(p):
    if len(p) < 4:
        return -99.0
    Q = QL
    t = 0.0
    for s in range(len(p) - 3):
        t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
    return t / (len(p) - 3)


def letters_idx(text):
    return [H.LIDX[c] for c in H.fold_base(text.replace('_', ''))]


def ic(vals):
    n = len(vals)
    c = Counter(vals)
    return sum(v * (v - 1) for v in c.values()) / (n * (n - 1))


# ============================================================== the three tests
def test_l1(loops):
    return dict(stat=ic(triples(loops)), distinct=len(set(triples(loops))))


def l2_readings(loops):
    for o in ORDERS:
        lo = order_loops(loops, o)
        for sig, tr in (('first digit most significant', triples(lo)),
                        ('last digit most significant', triples_lsb(lo))):
            for rot in range(27):
                text = ''.join(ALPHA27[(v + rot) % 27] for v in tr)
                yield f"{o}, {sig}, alphabet rotated {rot}", text


def test_l2(loops):
    best = (-1e9, '', '')
    for lab, text in l2_readings(loops):
        s = qscore(letters_idx(text))
        if s > best[0]:
            best = (s, lab, text)
    return dict(stat=best[0], label=best[1], text=best[2],
                words=H.word_score(best[2].replace('_', '').translate(H.TOL)))


def anneal_small(sym, restarts, iters, rng, T0=10.0, T1=0.3):
    """Any key for a short symbol sequence. Injective while symbols <= 24 letters;
    if a scrambled copy has more than 24 distinct triples, the extras may share letters."""
    n = len(sym)
    uniq = sorted(set(sym))
    idx = {u: i for i, u in enumerate(uniq)}
    c = [idx[x] for x in sym]
    S = len(uniq)
    fac = (T1 / T0) ** (1.0 / iters)
    best_all = (-1e18, None)
    for _ in range(restarts):
        key = list(range(24))
        rng.shuffle(key)
        key += [rng.randrange(24) for _ in range(max(0, S - 24))]
        p = [key[x] for x in c]
        cur = qscore(p) * (n - 3)
        best, bp = cur, p[:]
        Tm = T0
        for _ in range(iters):
            i = rng.randrange(S)
            if i >= 24:
                old = key[i]
                key[i] = rng.randrange(24)
                j = None
            else:
                j = rng.randrange(23)
                if j >= i:
                    j += 1
                key[i], key[j] = key[j], key[i]
            p2 = [key[x] for x in c]
            s2 = qscore(p2) * (n - 3)
            d = s2 - cur
            if d >= 0 or rng.random() < math.exp(d / Tm):
                cur, p = s2, p2
                if cur > best:
                    best, bp = cur, p[:]
            else:
                if j is None:
                    key[i] = old
                else:
                    key[i], key[j] = key[j], key[i]
            Tm *= fac
        if best > best_all[0]:
            best_all = (best, bp)
    return best_all[0] / (n - 3), best_all[1]


def test_l3(loops, seed, restarts, iters):
    rng = random.Random(seed)
    best = (-1e9, None, '')
    for o in ORDERS:
        tr = triples(order_loops(loops, o))
        s, p = anneal_small(tr, restarts, iters, rng)
        if s > best[0]:
            best = (s, o, ''.join(H.LETTERS[x] for x in p))
    return dict(stat=best[0], label=best[1], text=best[2],
                words=H.word_score(best[2].translate(H.TOL)))


# ============================================================== nulls and controls
def null_loops(loops, kind, seed):
    rng = random.Random(seed)
    if kind == 'shuffle':
        s = list(loops)
        rng.shuffle(s)
        return s
    return M.doublet_shuffle(list(loops), rng)


def plant(test, k, seed):
    """29 characters of period English written as loop-count triples."""
    rng = random.Random(seed)
    raw = M.CONTROL_TEXTS[k % len(M.CONTROL_TEXTS)].upper()
    if test == 'L2':
        text = ''.join(ch if ch.isalpha() else '_' for ch in raw)
        text = ' '.join(text.split('_')).replace(' ', '_')[:29]
        o = rng.choice(ORDERS)
        lsb = rng.random() < 0.5
        rot = rng.randrange(27)
        vals = [(ALPHA27.index(ch) - rot) % 27 for ch in text]
        lo = [t for v in vals for t in (trits_of(v)[::-1] if lsb else trits_of(v))]
        inv = list(range(87))
        R = T.apply_order(inv, o)
        loops = [0] * 87
        for i, r in enumerate(R):
            loops[r] = lo[i]
        return loops, text
    text = ''.join(ch for ch in raw if ch.isalpha())[:29]
    letters = sorted(set(text))
    code = dict(zip(letters, rng.sample(range(27), len(letters))))
    loops = [t for ch in text for t in trits_of(code[ch])]
    return loops, text


# ============================================================== jobs
def job(test, loops, seed, restarts, iters):
    if test == 'L1':
        return test_l1(loops)
    if test == 'L2':
        return test_l2(loops)
    return test_l3(loops, seed, restarts, iters)


def derive_seed(*parts):
    return random.Random('|'.join(str(x) for x in parts)).getrandbits(48)


def pval(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


def pct(vals, q):
    v = sorted(vals)
    return v[min(len(v) - 1, int(q * len(v)))]


def run_jobs(pool, args_list):
    if pool is None:
        return [job(*a) for a in args_list]
    res = pool.starmap_async(job, args_list, chunksize=4)
    while True:
        try:
            return res.get(60)
        except mp.TimeoutError:
            continue


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='looptrits_out')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--nulls', type=int, default=2000, help='copies of each kind for L1 and L2')
    ap.add_argument('--nulls-l3', type=int, default=500, help='copies of each kind for L3')
    ap.add_argument('--controls', type=int, default=100)
    ap.add_argument('--real-seeds', type=int, default=5)
    ap.add_argument('--restarts', type=int, default=8)
    ap.add_argument('--iters', type=int, default=8000)
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if args.selftest:
        args.nulls, args.nulls_l3, args.controls, args.real_seeds = 100, 20, 10, 2

    os.makedirs(args.out, exist_ok=True)
    logf = open(os.path.join(args.out, 'log.txt'), 'a', encoding='utf-8')

    def log(msg=''):
        print(msg, flush=True)
        logf.write(msg + '\n')
        logf.flush()

    t0 = time.time()
    init_tables(args.corpus)
    loops = load_loops()
    log(f"\n[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] START looptrits"
        f"{' SELFTEST' if args.selftest else ''}  workers {args.workers}")
    log("loop counts: " + ''.join(map(str, loops)))
    tr = triples(loops)
    log(f"29 triples (base 3, first digit most significant): {' '.join(map(str, tr))}")
    log(f"distinct triples: {len(set(tr))}; most common: "
        + ', '.join(f"{''.join(map(str, trits_of(v)))} x{n}" for v, n in Counter(tr).most_common(5)))

    pool = mp.Pool(args.workers, initializer=init_worker, initargs=(args.corpus,)) if args.workers > 1 else None
    results = {}
    try:
        for test in TESTS:
            nn = args.nulls_l3 if test == 'L3' else args.nulls
            nseeds = args.real_seeds if test == 'L3' else 1
            jobs = [(test, loops, derive_seed(args.seed, 'real', test, s), args.restarts, args.iters)
                    for s in range(nseeds)]
            for kind in ('shuffle', 'doublet'):
                jobs += [(test, null_loops(loops, kind, derive_seed(args.seed, kind, test, i)),
                          derive_seed(args.seed, 'nullsearch', test, kind, i), args.restarts, args.iters)
                         for i in range(nn)]
            plants = [plant(test, i, derive_seed(args.seed, 'plant', test, i)) for i in range(args.controls)]
            jobs += [(test, pl[0], derive_seed(args.seed, 'ctrl', test, i), args.restarts, args.iters)
                     for i, pl in enumerate(plants)]
            ts = time.time()
            out = run_jobs(pool, jobs)
            real = out[:nseeds]
            sh = out[nseeds:nseeds + nn]
            db = out[nseeds + nn:nseeds + 2 * nn]
            ct = out[nseeds + 2 * nn:]
            results[test] = (real, sh, db, ct, plants)
            log(f"[{datetime.datetime.now():%H:%M:%S}] {test} done in {H.fmt_secs(time.time() - ts)}")
    except KeyboardInterrupt:
        log("Stopped by Ctrl+C; reporting what finished.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()

    log("\n" + "=" * 90)
    log("RESULTS  (p vs pair-preserving copies is the one that matters; "
        f"Bonferroni line p < {BONF:.4f})")
    log("=" * 90)
    for test, (real, sh, db, ct, plants) in results.items():
        rstat = statistics.median(r['stat'] for r in real)
        dn = [x['stat'] for x in db]
        sn = [x['stat'] for x in sh]
        cs = [x['stat'] for x in ct]
        q95 = pct(dn, 0.95)
        power = sum(1 for c in cs if c > q95) / len(cs)
        p_d, p_s = pval(rstat, dn), pval(rstat, sn)
        if p_d <= BONF and p_s <= BONF:
            verdict = 'BEATS both nulls: look closely' + (' (but low power)' if power < 0.5 else '')
        elif p_d <= 0.05:
            verdict = 'below 0.05, not past Bonferroni'
        elif power < 0.5:
            verdict = 'chance level, LOW POWER: inconclusive'
        else:
            verdict = 'chance level'
        name = {'L1': 'letter-likeness (IC of triples)', 'L2': 'counting keys (270 readings)',
                'L3': 'any key (annealing)'}[test]
        log(f"\n{test}  {name}")
        seeds_txt = ', '.join('%.3f' % r['stat'] for r in real)
        log(f"  real {rstat:.4f}" + (f"  (seeds: {seeds_txt})" if len(real) > 1 else ''))
        log(f"  pair-preserving: n={len(dn)}, mean {statistics.mean(dn):.4f}, 95th pct {q95:.4f}, "
            f"max {max(dn):.4f}, p = {p_d:.4f}")
        log(f"  plain shuffle:   n={len(sn)}, mean {statistics.mean(sn):.4f}, max {max(sn):.4f}, p = {p_s:.4f}")
        log(f"  planted English: median {statistics.median(cs):.4f}; power (above null 95th pct) {power:.0%}")
        if test != 'L1':
            b = max(real, key=lambda r: r['stat'])
            log(f"  best real reading: [{b['label']}] {b['text']}   words {b['words']}")
            wd = [x['words'] for x in db]
            log(f"  word score p vs pair-preserving: {pval(b['words'], wd):.4f}")
            ex = max(range(len(ct)), key=lambda i: ct[i]['stat'])
            log(f"  example control: planted {plants[ex][1]}  ->  found {ct[ex]['text']}")
        else:
            log(f"  real distinct triples {real[0]['distinct']}; pair-preserving mean distinct "
                f"{statistics.mean(x['distinct'] for x in db):.1f}")
            # diagnostic (not pre-registered): the same IC with triples starting at symbol 2 or 3,
            # compared with the same frames in the pair-preserving copies
            for fr in (1, 2):
                rf = ic(triples(loops[fr:]))
                nf = [ic(triples(null_loops(loops, 'doublet', derive_seed(args.seed, 'doublet', test, i))[fr:]))
                      for i in range(len(db))]
                log(f"  diagnostic, triples starting at symbol {fr + 1}: IC {rf:.4f}, p vs pair-preserving "
                    f"{pval(rf, nf):.4f} (a real message in frame 1 should look random here)")
        log(f"  verdict: {verdict}")
    log(f"\n[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] END looptrits, {H.fmt_secs(time.time() - t0)}")
    if results:
        with open(os.path.join(args.out, 'report.txt'), 'w', encoding='utf-8') as fh:
            with open(os.path.join(args.out, 'log.txt'), encoding='utf-8') as lf:
                txt = lf.read()
            fh.write(txt[txt.rfind('START looptrits') - 22:])


if __name__ == '__main__':
    mp.freeze_support()
    main()
