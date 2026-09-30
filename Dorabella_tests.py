#!/usr/bin/env python3
"""
Dorabella test suite. Needs Dorabella_hillclimb.py and english.txt in the same folder.
  python Dorabella_tests.py | Tee-Object Dorabella_tests_out.txt
"""
import argparse
import heapq
import itertools
import multiprocessing as mp
import os
import random
import statistics
from statistics import NormalDist

import Dorabella_hillclimb as H

# position (1-based) -> (as transcribed, alternate)
AMBIG = {23: (15, 14), 32: (11, 18), 34: (11, 12), 78: (15, 14), 86: (13, 12)}
ORDERS = ['forward', 'fully reversed', 'each line reversed',
          'boustrophedon A', 'boustrophedon B']
CANON = ('triplet', (0, 1, 2), 2, 1)   # A,B,C = 13,23,33; letters advance clockwise
TABS = None


def get_tabs():
    """Key family: keypad layout (loop-order x rotation x direction) + affine shifts of it."""
    global TABS
    if TABS is None:
        tabs = []
        for perm in itertools.permutations(range(3)):
            for s in range(8):
                for sense in (1, -1):
                    key = {}
                    for loop in (1, 2, 3):
                        for d in range(1, 9):
                            t = ((d - 1 - s) * sense) % 8
                            key[loop * 10 + d] = 3 * t + perm[loop - 1]
                    tabs.append((('triplet', perm, s, sense), key))
        for a in (1, 5, 7, 11, 13, 17, 19, 23):
            for b in range(24):
                if a == 1 and b == 0:
                    continue
                key = {}
                for loop in (1, 2, 3):
                    for d in range(1, 9):
                        x = 3 * ((d - 3) % 8) + (loop - 1)
                        key[loop * 10 + d] = (a * x + b) % 24
                tabs.append((('affine', a, b), key))
        TABS = tabs
    return TABS


def apply_order(seq, name):
    l = [seq[a:b] for a, b in H.LINE_SPANS]
    if name == 'forward':
        return list(seq)
    if name == 'fully reversed':
        return seq[::-1]
    if name == 'each line reversed':
        return l[0][::-1] + l[1][::-1] + l[2][::-1]
    if name == 'boustrophedon A':
        return l[0] + l[1][::-1] + l[2]
    return l[0][::-1] + l[1] + l[2][::-1]


def all_seqs(base):
    positions = sorted(AMBIG)
    for bits in itertools.product((0, 1), repeat=len(positions)):
        s = list(base)
        changed = []
        for pos, b in zip(positions, bits):
            if b:
                s[pos - 1] = AMBIG[pos][1]
                changed.append(f"{pos}={AMBIG[pos][1]}")
        vlabel = ','.join(changed) or 'as transcribed'
        for oname in ORDERS:
            yield f"{oname} | {vlabel}", apply_order(s, oname)


def scan(seqs, tabs, keep=0):
    best = (-1e18, None, None, None)
    top = []
    n = 0
    for label, seq in seqs:
        for cfg, key in tabs:
            p = [key[c] for c in seq]
            s = H.SCORE_FN(p)
            n += 1
            if s > best[0]:
                best = (s, label, cfg, p)
            if keep:
                if len(top) < keep:
                    heapq.heappush(top, (s, n, label, cfg, p))
                elif s > top[0][0]:
                    heapq.heapreplace(top, (s, n, label, cfg, p))
    return best, top, n


def null_scan(seed, n, expanded):
    """Null: shuffle the symbols, then run the IDENTICAL search; return best score of each."""
    tabs = get_tabs()
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        s = list(H.BASE_SEQ)
        rng.shuffle(s)
        seqs = list(all_seqs(s)) if expanded else [('forward', s)]
        out.append(scan(seqs, tabs)[0][0])
    return out


def letters_of(p):
    return ''.join(H.LETTERS[i] for i in p)


def show(text):
    H.print_lines(text)
    print("   segmented (folded): " + H.segment(text.translate(H.TOL)))


def part_c(args, pool, workers, timer):
    print("\n=== PART C: KEYPAD-LAYOUT KEYS ===")
    tabs = get_tabs()
    n_keys = len(tabs) + 1
    base = list(H.BASE_SEQ)
    canon_key = dict(next(k for c, k in tabs if c == CANON))
    scored = []
    for cfg, key in tabs:
        p = [key[c] for c in base]
        scored.append((H.SCORE_FN(p), cfg, p))
    scored.sort(key=lambda x: -x[0])
    timer.lap("level 1 scan", n_keys, "keys")

    rank = next(i for i, (s, cfg, p) in enumerate(scored, 1) if cfg == CANON)
    print(f"\nYour recollected key (A,B,C = 13,23,33): score {scored[rank - 1][0]:.4f}, "
          f"rank {rank} of {len(scored)}")
    show(letters_of(scored[rank - 1][2]))
    print(f"\nTop 5 of {len(scored)} keys on the as-transcribed cipher:")
    for r, (s, cfg, p) in enumerate(scored[:5], 1):
        text = letters_of(p)
        print(f"\n#{r} {cfg}  score={s:.4f}  words={H.word_score(text.translate(H.TOL))}")
        show(text)

    counts = H.split_counts(args.c1_shuffles, workers)
    jobs = [(args.seed + 31 * w, n, False) for w, n in enumerate(counts) if n > 0]
    null1 = [x for r in H.run_parallel(pool, null_scan, jobs) for x in r]
    timer.lap("level 1 null", len(null1) * len(tabs), "keys")
    print(f"\nLEVEL 1 ({len(tabs)} keys, as transcribed): best {scored[0][0]:.4f}; "
          f"null mean {statistics.mean(null1):.4f}, max {max(null1):.4f}; "
          f"p = {H.pvalue(scored[0][0], null1):.4f} (n={len(null1)})")

    seqs = list(all_seqs(base))
    best, top, n = scan(seqs, tabs, keep=args.top)
    timer.lap("level 2 scan", n, "key-decodes")
    print(f"\nLevel 2: {len(seqs)} sequences x {len(tabs)} keys. Top {args.top}:")
    for r, (s, _, label, cfg, p) in enumerate(sorted(top, reverse=True), 1):
        text = letters_of(p)
        print(f"\n#{r} {label}  {cfg}  score={s:.4f}  words={H.word_score(text.translate(H.TOL))}")
        show(text)
    counts = H.split_counts(args.c2_shuffles, workers)
    jobs = [(args.seed + 977 * w, k, True) for w, k in enumerate(counts) if k > 0]
    null2 = [x for r in H.run_parallel(pool, null_scan, jobs) for x in r]
    timer.lap("level 2 null", len(null2) * n, "key-decodes")
    print(f"\nLEVEL 2 (all orders, variants, keys): best {best[0]:.4f}; "
          f"null mean {statistics.mean(null2):.4f}, max {max(null2):.4f}; "
          f"p = {H.pvalue(best[0], null2):.4f} (n={len(null2)})")
    print("\nHow to read Part C: p < ~0.02 at either level would be the first real signal so far; "
          "then read the decode by eye. The nulls repeat the identical search on shuffled symbols, "
          "so multiple looks are already corrected for.")


def plant(rng):
    plain = H.fold_base(''.join(c for c in H.PLANT.upper() if c.isalpha()))[:87]
    letters = sorted(set(plain))
    codes = rng.sample(H.CODES_ALL, len(letters))
    code_of = dict(zip(letters, codes))
    return plain, [code_of[c] for c in plain]


def corrupt(seq, k, rng):
    s = list(seq)
    for pos in rng.sample(range(len(s)), k):
        loop, d = divmod(s[pos], 10)
        d = (d - 1 + rng.choice((-1, 1))) % 8 + 1
        s[pos] = loop * 10 + d
    return s


def part_a(args, pool, workers, timer):
    print("\n=== PART A: DOES THE SOLVER SURVIVE TRANSCRIPTION ERRORS? ===")
    print("Planted English, random key; k symbols nudged one direction step.")
    print(f"{'k':>3} {'mean recovery':>14} {'worst':>7} {'best':>7} {'ceiling':>8}")
    evals = 0
    for k in args.corrupt:
        accs = []
        for t in range(args.trials):
            rng = random.Random(args.seed + 100 * k + t)
            plain, seq = plant(rng)
            seq = corrupt(seq, k, rng)
            uniq, cidx = H.prep_cipher(seq)
            finals, ev = H.search_parallel(pool, workers, args.seed + t, cidx,
                                           len(uniq), args.restarts)
            evals += ev
            _, key = max(finals, key=lambda x: x[0])
            text = H.decrypt_text(cidx, key)
            accs.append(sum(a == b for a, b in zip(text, plain)) / len(plain))
        print(f"{k:>3} {statistics.mean(accs):>13.0%} {min(accs):>7.0%} "
              f"{max(accs):>7.0%} {(87 - k) / 87:>8.0%}")
    timer.lap("Part A", evals, "score-evals")
    print("\nHow to read: recovery near the ceiling means a few misread symbols can't explain a "
          "null result. Collapse at k=2-4 means our nulls say little about your transcription.")


def part_b(args, pool, workers, timer):
    print("\n=== PART B: HILL-CLIMB OVER AMBIGUOUS SYMBOLS x READING ORDERS ===")
    seqs = list(all_seqs(H.BASE_SEQ))
    n_cfg = len(seqs)
    results, evals = [], 0
    for i, (label, seq) in enumerate(seqs):
        uniq, cidx = H.prep_cipher(seq)
        finals, ev = H.search_parallel(pool, workers, args.seed + i, cidx,
                                       len(uniq), args.restarts)
        evals += ev
        s, key = max(finals, key=lambda x: x[0])
        results.append((s, label, H.decrypt_text(cidx, key)))
    timer.lap(f"{n_cfg} real hill-climbs", evals, "score-evals")
    results.sort(key=lambda r: -r[0])

    uniq0, cidx0 = H.prep_cipher(H.BASE_SEQ)
    sh = H.run_parallel(pool, H.shuffle_job,
                        [(args.seed + 5000 + i, cidx0, len(uniq0), args.restarts)
                         for i in range(args.shuffles)])
    timer.lap("Part B null", sum(x[2] for x in sh), "score-evals")
    null = [x[0] for x in sh]
    nd = NormalDist(statistics.mean(null), statistics.pstdev(null) or 1.0)
    best = results[0][0]
    z = (best - nd.mean) / nd.stdev
    p_raw = 1 - nd.cdf(best)
    print(f"\nBest config score {best:.4f}; single-search null mean {nd.mean:.4f}, sd {nd.stdev:.4f}")
    print(f"z = {z:.2f}; normal-approx p = {p_raw:.4g}; corrected for {n_cfg} configs: "
          f"{min(1.0, p_raw * n_cfg):.4g}")
    print(f"z needed for corrected p = 0.02: {NormalDist().inv_cdf(1 - 0.02 / n_cfg):.2f}")
    for r, (s, label, text) in enumerate(results[:args.top], 1):
        print(f"\n#{r} {label}  score={s:.4f}  words={H.word_score(text.translate(H.TOL))}")
        show(text)
    print("\nHow to read Part B: only a z above the 'needed' value counts as signal here.")


def main():
    timer = H.RunTimer("dorabella_tests")
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument('--corpus', default='english.txt')
        ap.add_argument('--restarts', type=int, default=80)
        ap.add_argument('--shuffles', type=int, default=100)
        ap.add_argument('--c1-shuffles', type=int, default=2000)
        ap.add_argument('--c2-shuffles', type=int, default=300)
        ap.add_argument('--trials', type=int, default=5)
        ap.add_argument('--corrupt', type=int, nargs='+', default=[0, 2, 4, 6, 8])
        ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
        ap.add_argument('--top', type=int, default=5)
        ap.add_argument('--seed', type=int, default=1897)
        ap.add_argument('--skip-a', action='store_true')
        ap.add_argument('--skip-b', action='store_true')
        ap.add_argument('--skip-c', action='store_true')
        args = ap.parse_args()

        if not os.path.exists(args.corpus):
            raise SystemExit(f"corpus not found: {args.corpus}")
        for pos, (a, _) in AMBIG.items():
            assert H.BASE_SEQ[pos - 1] == a, f"position {pos} is not {a}"

        H.build_scorer(args.corpus)
        H.build_words(None)
        workers = max(1, args.workers)
        timer.lap("setup")
        pool = None
        if workers > 1:
            pool = mp.Pool(workers, initializer=H.init_worker, initargs=(args.corpus, None))
            timer.lap("worker pool created")
        try:
            if not args.skip_c:
                part_c(args, pool, workers, timer)
            if not args.skip_a:
                part_a(args, pool, workers, timer)
            if not args.skip_b:
                part_b(args, pool, workers, timer)
        finally:
            if pool is not None:
                pool.close()
                pool.join()
    finally:
        timer.finish()


if __name__ == '__main__':
    mp.freeze_support()
    main()