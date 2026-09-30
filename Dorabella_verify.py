#!/usr/bin/env python3
"""
Dorabella verification: is Part B's z=4.98 real, or an artifact of a mismatched null?
Needs Dorabella_hillclimb.py, Dorabella_tests.py and english.txt in the same folder.

  python Dorabella_verify.py | Tee-Object Dorabella_verify_out.txt
  python Dorabella_verify.py --cal 10 --reps 8        (quicker, weaker check)

Fixes two flaws in Part B's null:
 1. Matched null: each of the 32 symbol-variant frequency profiles gets its own
    shuffled-null mean/sd (variants change symbol frequencies, which alone moves scores).
 2. Identical-search null: the SAME 160-config search (32 variants x 5 reading orders)
    is repeated on shuffled copies; the real best is compared with those bests.
    Ambiguous symbols travel with their alternates when shuffled.
"""
import argparse
import itertools
import multiprocessing as mp
import os
import random
import statistics

import Dorabella_hillclimb as H
import Dorabella_tests as T

POSITIONS = sorted(T.AMBIG)                       # 1-based ambiguous positions
BITS = list(itertools.product((0, 1), repeat=len(POSITIONS)))


def base_tokens():
    """(code, alternate, ambiguity index or -1) for each of the 87 symbols."""
    toks = []
    for i, c in enumerate(H.BASE_SEQ, 1):
        if i in T.AMBIG:
            toks.append((c, T.AMBIG[i][1], POSITIONS.index(i)))
        else:
            toks.append((c, c, -1))
    return toks


def shuffled_tokens(seed):
    rng = random.Random(seed)
    t = base_tokens()
    rng.shuffle(t)
    return t


def realize(tokens, bits):
    return [alt if (a >= 0 and bits[a]) else code for code, alt, a in tokens]


def vlabel(bits):
    ch = [f"{p}={T.AMBIG[p][1]}" for p, b in zip(POSITIONS, bits) if b]
    return ','.join(ch) or 'as transcribed'


def climb_job(seed, seq, restarts):
    uniq, cidx = H.prep_cipher(seq)
    finals, ev = H.run_restarts(seed, cidx, len(uniq), restarts)
    return max(s for s, _ in finals), ev


def real_job(seed, seq, restarts):
    uniq, cidx = H.prep_cipher(seq)
    finals, ev = H.run_restarts(seed, cidx, len(uniq), restarts)
    s, key = max(finals, key=lambda x: x[0])
    return s, H.decrypt_text(cidx, key), ev


def order_means(zs, configs):
    out = {}
    for o in T.ORDERS:
        out[o] = statistics.mean(z for z, (_, oo) in zip(zs, configs) if oo == o)
    return out


def main():
    timer = H.RunTimer("dorabella_verify")
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument('--corpus', default='english.txt')
        ap.add_argument('--restarts', type=int, default=80)
        ap.add_argument('--cal', type=int, default=16, help='shuffles per variant (matched null)')
        ap.add_argument('--reps', type=int, default=20, help='full-search replicates on shuffles')
        ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
        ap.add_argument('--seed', type=int, default=1897)
        ap.add_argument('--top', type=int, default=6)
        args = ap.parse_args()

        if not os.path.exists(args.corpus):
            raise SystemExit(f"corpus not found: {args.corpus}")
        for pos, (a, _) in T.AMBIG.items():
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
            configs = [(vi, o) for vi in range(len(BITS)) for o in T.ORDERS]

            # ---- real data: all 160 configs
            toks = base_tokens()
            jobs = [(args.seed + i, T.apply_order(realize(toks, BITS[vi]), o), args.restarts)
                    for i, (vi, o) in enumerate(configs)]
            real = H.run_parallel(pool, real_job, jobs)
            timer.lap("real 160 configs", sum(r[2] for r in real), "score-evals")

            # ---- matched null per variant
            cal_jobs = []
            for vi in range(len(BITS)):
                for k in range(args.cal):
                    t = shuffled_tokens(args.seed + 10_000 + vi * 1000 + k)
                    cal_jobs.append((args.seed + 50_000 + vi * 1000 + k,
                                     realize(t, BITS[vi]), args.restarts))
            res = H.run_parallel(pool, climb_job, cal_jobs)
            timer.lap("matched null calibration", sum(r[1] for r in res), "score-evals")
            mu, sd, idx = {}, {}, 0
            for vi in range(len(BITS)):
                sc = [res[idx + k][0] for k in range(args.cal)]
                idx += args.cal
                mu[vi] = statistics.mean(sc)
                sd[vi] = statistics.pstdev(sc) or 1.0

            zs = [(real[i][0] - mu[vi]) / sd[vi] for i, (vi, o) in enumerate(configs)]
            order = sorted(range(len(configs)), key=lambda i: -zs[i])
            print(f"\n=== REAL CONFIGS RANKED BY VARIANT-MATCHED z (top {args.top}) ===")
            for rank, i in enumerate(order[:args.top], 1):
                vi, o = configs[i]
                text = real[i][1]
                print(f"\n#{rank} {o} | {vlabel(BITS[vi])}   raw={real[i][0]:.4f}   "
                      f"matched z={zs[i]:.2f}   words={H.word_score(text.translate(H.TOL))}")
                T.show(text)
            real_best_z = zs[order[0]]
            real_best_raw = max(r[0] for r in real)
            om = order_means(zs, configs)
            print("\nMean matched z by reading order (real data, over the 32 variants):")
            for o in T.ORDERS:
                print(f"  {o:<20} {om[o]:+.2f}")
            print(f"\nReal best matched z = {real_best_z:.2f}   "
                  f"(Part B's unmatched z was 4.98)")
            print("If this is below about 3, Part B was mostly a null-mismatch artifact; "
                  "you can stop here (Ctrl+C) and skip the replicates.\n", flush=True)

            # ---- identical-search replicates on shuffled data
            rep_z, rep_raw, rep_om = [], [], {o: [] for o in T.ORDERS}
            evals = 0
            for r in range(args.reps):
                t = shuffled_tokens(args.seed + 900_000 + r)
                jobs = [(args.seed + 200_000 + r * 1000 + i,
                         T.apply_order(realize(t, BITS[vi]), o), args.restarts)
                        for i, (vi, o) in enumerate(configs)]
                rs = H.run_parallel(pool, climb_job, jobs)
                evals += sum(x[1] for x in rs)
                z = [(rs[i][0] - mu[vi]) / sd[vi] for i, (vi, o) in enumerate(configs)]
                rep_z.append(max(z))
                rep_raw.append(max(x[0] for x in rs))
                for o, m in order_means(z, configs).items():
                    rep_om[o].append(m)
                print(f"  replicate {r + 1}/{args.reps}: best matched z = {max(z):.2f}", flush=True)
            timer.lap("replicates", evals, "score-evals")

            beat = sum(1 for z in rep_z if z >= real_best_z)
            beat_raw = sum(1 for s in rep_raw if s >= real_best_raw)
            print("\n=== REAL vs IDENTICAL SEARCH ON SHUFFLED DATA ===")
            print(f"Real best matched z: {real_best_z:.2f}")
            print(f"Shuffled bests (n={args.reps}): mean {statistics.mean(rep_z):.2f}, "
                  f"max {max(rep_z):.2f}")
            print(f"  -> p = {(1 + beat) / (1 + args.reps):.4f}  "
                  f"({beat} of {args.reps} shuffled searches matched or beat the real best; "
                  f"smallest possible p with this many replicates is {1 / (1 + args.reps):.3f})")
            print(f"Raw check: real best raw {real_best_raw:.4f}; shuffled bests mean "
                  f"{statistics.mean(rep_raw):.4f}, max {max(rep_raw):.4f}; "
                  f"{beat_raw} matched or beat it")
            if args.reps > 1:
                print("Per-order mean z: real vs spread across shuffled replicates (sd):")
                for o in T.ORDERS:
                    print(f"  {o:<20} real {om[o]:+.2f}   shuffled sd {statistics.pstdev(rep_om[o]):.2f}")
            print("\nHow to read this:")
            print(" * Real best z far above every shuffled best, AND words score above shuffles: "
                  "worth a much larger replication (50+ replicates) and a held-out test.")
            print(" * Real best z inside the shuffled range: Part B was an artifact of the "
                  "mismatched null.")
            print(" * A reading order whose real mean z stands well outside its shuffled sd "
                  "is a lead on direction, not proof of a solution.")
            timer.lap("report")
        finally:
            if pool is not None:
                pool.close()
                pool.join()
    finally:
        timer.finish()


if __name__ == '__main__':
    mp.freeze_support()
    main()