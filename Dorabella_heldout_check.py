#!/usr/bin/env python3
"""
Recheck of one held-out result from the marathon.

The held-out test learns a key on two lines and scores the third line with it. The
learned key depends on the climb's random seed, so one run of the real cipher can be
lucky. This script runs ONE configuration (default: each line reversed, with
32=18, 34=12, 78=14, 86=12, the marathon's pre-registered config) as follows:

  1. the real cipher with many different seeds, to give its typical score;
  2. many plain shuffles and pair-preserving shuffles, one seed each;
  3. a p-value for the real cipher's median, best and worst score against each.

Needs Dorabella_hillclimb.py, Dorabella_tests.py, Dorabella_marathon.py and english.txt
in the same folder.

  python Dorabella_heldout_check.py                       # 24 seeds, 240 copies of each kind
  python Dorabella_heldout_check.py --seeds 12 --nulls 100

Output: screen + heldout_check_out/log.txt
"""
import argparse
import datetime
import os
import statistics
import time
from multiprocessing import Pool

import Dorabella_hillclimb as H
import Dorabella_marathon as M
import Dorabella_tests as T

LOG = None
ORDER = 'each line reversed'
BITS = (0, 1, 1, 1, 1)          # positions 23, 32, 34, 78, 86: 1 = alternate reading
TOKS = M.make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})


def log(msg=''):
    print(msg, flush=True)
    LOG.write(msg + '\n')
    LOG.flush()


def real_job(seed):
    lines = M.oriented_lines(M.realize(TOKS, BITS), ORDER)
    return M.heldout_cfg(seed, lines, 80)[0]


def null_job(args):
    kind, k = args
    t = M.null_tokens(kind, TOKS, 50000 + k)
    lines = M.oriented_lines(M.realize(t, BITS), ORDER)
    return M.heldout_cfg(9000 + k, lines, 80)[0]


def p_of(real, nulls):
    return (1 + sum(x >= real for x in nulls)) / (1 + len(nulls))


def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--seeds', type=int, default=24)
    ap.add_argument('--nulls', type=int, default=240)
    ap.add_argument('--workers', type=int, default=os.cpu_count())
    args = ap.parse_args()
    os.makedirs('heldout_check_out', exist_ok=True)
    LOG = open(os.path.join('heldout_check_out', 'log.txt'), 'a', encoding='utf-8')
    alts = ','.join(f"{p}={T.AMBIG[p][1]}" for p, b in zip(sorted(T.AMBIG), BITS) if b)
    log(f"\n[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] START Dorabella_heldout_check  "
        f"config [{ORDER} | {alts}]  seeds={args.seeds} nulls={args.nulls}")
    log("Marathon's single real run of this config: -13.8362")

    with Pool(args.workers, initializer=M.init_worker, initargs=(args.corpus,)) as pool:
        t0 = time.time()
        real = sorted(pool.map(real_job, range(100, 100 + args.seeds)))
        med = statistics.median(real)
        log(f"real cipher, {args.seeds} seeds: median {med:.3f}, worst {real[0]:.3f}, best {real[-1]:.3f}  "
            f"({M.fmt_secs(time.time() - t0)})")
        log("   " + ' '.join(f"{x:.2f}" for x in real))
        for kind in ('shuffle', 'doublet'):
            t1 = time.time()
            nulls = pool.map(null_job, [(kind, k) for k in range(args.nulls)])
            log(f"vs {kind} ({len(nulls)} copies, median {statistics.median(nulls):.3f}, "
                f"max {max(nulls):.3f}):  p at real median {p_of(med, nulls):.4f} | "
                f"at real best {p_of(real[-1], nulls):.4f} | at real worst {p_of(real[0], nulls):.4f}  "
                f"({M.fmt_secs(time.time() - t1)})")
    log("The p at the real MEDIAN is the fair one; the best-of-seeds p shows how a lucky run looks.")
    LOG.close()


if __name__ == '__main__':
    main()
