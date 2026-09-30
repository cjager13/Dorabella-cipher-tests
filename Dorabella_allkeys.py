#!/usr/bin/env python3
"""
Exhaustive test of every 'sets of three' key for the Dorabella.

Key model: the 24 letters form 8 sets of three. Each set sits on one of the 8
directions (8! = 40,320 orders) and is read 1->2->3 loops either "up" or "down"
(2^8 = 256 patterns). Total: 10,321,920 keys, and every one is scored.

Two letter groupings are tested:
  yours    : ABC DEF GHI KLM NOP QRS TUW XYZ   (I=J, U=V)
  standard : ABC DEF GHI JKL MNO PQR STU WXY   (the published key lives in here;
             the script checks this at startup)

Needs Dorabella_hillclimb.py, Dorabella_tests.py, english.txt, and numpy
(pip install numpy).

  python Dorabella_allkeys.py                        # forward reading, 200 shuffled copies per grouping
  python Dorabella_allkeys.py --all-configs          # all 5 reading orders x 32 alternate readings
  python Dorabella_allkeys.py --all-configs --nulls 50
  python Dorabella_allkeys.py --grouping yours       # only one grouping

Output goes to the screen and to allkeys_out/log.txt (appended), plus
allkeys_out/top_keys.csv with the best keys found.

Scoring:
  n-gram score : quadgram English-likeness (higher = more English-like;
                 real English of this length ~ -9.2, shuffled Dorabella ~ -14 to -16)
  word score   : dictionary words found, longer words worth more. It is not what the
                 search optimises, so it is an independent check. It is computed on the
                 top --topk keys by n-gram score.
Both are compared with the IDENTICAL search on shuffled copies of the Dorabella,
so the p-values already account for trying 10 million keys.
"""
import argparse
import csv
import datetime
import itertools
import multiprocessing as mp
import os
import random
import statistics
import sys
import time

import numpy as np

import Dorabella_hillclimb as H
import Dorabella_tests as T

GROUPINGS = {
    'yours': ['ABC', 'DEF', 'GHI', 'KLM', 'NOP', 'QRS', 'TUW', 'XYZ'],
    'standard': ['ABC', 'DEF', 'GHI', 'JKL', 'MNO', 'PQR', 'STU', 'WXY'],
}
STANDARD_DECRYPT = ("BLTACEIARWUNISNFNNELLHSYWYDUOINIEYARQATNNTEDMINUNEHOMSYRRYUO"
                    "TOEHOTSHGDOTNEHMOSALDOEADYA")
PERMS = np.array(list(itertools.permutations(range(8))), dtype=np.int32)  # PERMS[p, s] = direction of set s
ROWS = np.arange(len(PERMS))
QT = None
LOG = None


# ------------------------------------------------------------------ logging
def log(msg=''):
    line = msg
    print(line, flush=True)
    if LOG:
        LOG.write(line + '\n')
        LOG.flush()


def stamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


# ------------------------------------------------------------------ core
def init(corpus):
    global QT
    H.build_scorer(corpus)
    H.build_words(None)
    QT = np.full(24 ** 4, H.QFLOOR, dtype=np.float64)
    for k, v in H.QUAD.items():
        QT[k] = v


def letters_of(sets):
    return [[H.LIDX[H.fold_base(ch)] for ch in s] for s in sets]


def key_tables(sets_idx, flip):
    """KL[p, (dir-1)*3 + (loop-1)] = letter index, for every direction order p."""
    KL = np.empty((len(PERMS), 24), dtype=np.int32)
    for s in range(8):
        down = (flip >> s) & 1
        for l in range(3):
            KL[ROWS, PERMS[:, s] * 3 + l] = sets_idx[s][2 - l if down else l]
    return KL


def code_index(seq):
    return np.array([(c % 10 - 1) * 3 + (c // 10 - 1) for c in seq], dtype=np.int64)


def search(seq, sets_idx, topk):
    """Top-k (score, perm, flip) over all 10,321,920 keys for one symbol sequence."""
    ci = code_index(seq)
    n = len(seq)
    best = []
    for flip in range(256):
        P = key_tables(sets_idx, flip)[:, ci]
        q = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
        S = QT[q].sum(axis=1) / (n - 3)
        k = min(topk, len(S))
        idx = np.argpartition(S, -k)[-k:]
        best.extend((float(S[i]), int(i), flip) for i in idx)
        best = sorted(best, reverse=True)[:topk]
    return best


def decode(seq, sets, p, flip):
    KL = key_tables(letters_of(sets), flip)[p]
    return ''.join(H.LETTERS[KL[i]] for i in code_index(seq))


def words(text):
    return H.word_score(text.translate(H.TOL))


def describe(sets, p, flip):
    """e.g. ABC^d3 = set ABC on direction 3, read up (A=1 loop); v = down (A=3 loops)."""
    return ' '.join(f"{sets[s]}{'v' if (flip >> s) & 1 else '^'}d{PERMS[p, s] + 1}" for s in range(8))


# ------------------------------------------------------------------ worker jobs
def job(label, seq, gname, topk):
    sets = GROUPINGS[gname]
    best = search(seq, letters_of(sets), topk)
    return [(s, label, p, f, words(decode(seq, sets, p, f))) for s, p, f in best]


def shuffled_configs(seed, all_configs):
    s = list(H.BASE_SEQ)
    random.Random(seed).shuffle(s)
    return list(T.all_seqs(s)) if all_configs else [('forward | as transcribed', s)]


def run_search(pool, configs, gname, topk):
    res = pool.starmap(job, [(lab, seq, gname, topk) for lab, seq in configs], chunksize=1)
    rows = sorted((r for rr in res for r in rr), reverse=True)[:topk]
    return rows          # (score, label, perm, flip, words)


# ------------------------------------------------------------------ checks
def check_standard():
    """The published key should be in the 'standard' space. Find the closest key."""
    sets = GROUPINGS['standard']
    tgt = np.array([H.LIDX[c] for c in H.fold_base(STANDARD_DECRYPT)])
    ci = code_index(H.BASE_SEQ)
    best = (0, 0, 0)
    for f in range(256):
        m = (key_tables(letters_of(sets), f)[:, ci] == tgt).sum(axis=1)
        i = int(m.argmax())
        if m[i] > best[0]:
            best = (int(m[i]), i, f)
    return best


def pval(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


# ------------------------------------------------------------------ main
def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='allkeys_out')
    ap.add_argument('--nulls', type=int, default=200, help='shuffled copies per grouping')
    ap.add_argument('--topk', type=int, default=50, help='keys kept (and word-scored) per search')
    ap.add_argument('--show', type=int, default=5, help='decrypts printed per grouping')
    ap.add_argument('--all-configs', action='store_true',
                    help='all 5 reading orders x 32 alternate readings (160x slower)')
    ap.add_argument('--grouping', choices=list(GROUPINGS) + ['both'], default='both')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--seed', type=int, default=1920)
    args = ap.parse_args()

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    os.makedirs(args.out, exist_ok=True)
    LOG = open(os.path.join(args.out, 'log.txt'), 'a', encoding='utf-8')
    t_start = time.time()
    log(f"\n[{stamp()}] START allkeys  {' '.join(sys.argv[1:]) or '(defaults)'}")
    log(f"workers {args.workers}, nulls {args.nulls}, topk {args.topk}, "
        f"{'all 160 configs' if args.all_configs else 'forward reading only'}")

    init(args.corpus)
    matched, p0, f0 = check_standard()
    log(f"check: closest 'standard' key to the published decrypt agrees at {matched}/87 positions "
        f"(the rest are ambiguous symbols)\n       {describe(GROUPINGS['standard'], p0, f0)}")

    configs = list(T.all_seqs(H.BASE_SEQ)) if args.all_configs else \
        [('forward | as transcribed', list(H.BASE_SEQ))]
    seqs = dict(configs)
    gnames = list(GROUPINGS) if args.grouping == 'both' else [args.grouping]

    csv_path = os.path.join(args.out, 'top_keys.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
        csv.writer(fh).writerow(['grouping', 'rank', 'ngram_score', 'word_score', 'config', 'key', 'decrypt'])

    with mp.Pool(args.workers, initializer=init, initargs=(args.corpus,)) as pool:
        for g in gnames:
            sets = GROUPINGS[g]
            log(f"\n{'=' * 78}\nGROUPING {g}: {' '.join(sets)}\n{'=' * 78}")
            t0 = time.time()
            rows = run_search(pool, configs, g, args.topk)
            log(f"[{stamp()}] real search: {len(configs)} sequence(s) x 10,321,920 keys "
                f"in {H.fmt_secs(time.time() - t0)}")

            with open(csv_path, 'a', newline='', encoding='utf-8') as fh:
                w = csv.writer(fh)
                for r, (s, lab, p, f, wd) in enumerate(rows, 1):
                    w.writerow([g, r, f"{s:.4f}", wd, lab, describe(sets, p, f), decode(seqs[lab], sets, p, f)])

            log("\nTop by n-gram score:")
            for s, lab, p, f, wd in rows[:args.show]:
                text = decode(seqs[lab], sets, p, f)
                log(f"\n  n-gram {s:.3f}  words {wd}  [{lab}]\n  key {describe(sets, p, f)}")
                for a, b in H.LINE_SPANS:
                    log("     " + text[a:b])
                log("     " + H.segment(text.translate(H.TOL)))
            by_words = sorted(rows, key=lambda r: (-r[4], -r[0]))
            log(f"\nMost words among the top {len(rows)} keys:")
            for s, lab, p, f, wd in by_words[:3]:
                text = decode(seqs[lab], sets, p, f)
                log(f"  words {wd}  n-gram {s:.3f}  [{lab}]  {H.segment(text.translate(H.TOL))}")

            real_ng = rows[0][0]
            real_wd_best = rows[0][4]
            real_wd_max = max(r[4] for r in rows)
            log(f"\n[{stamp()}] running the identical search on {args.nulls} shuffled copies...")
            null_ng, null_wd_best, null_wd_max = [], [], []
            t1 = time.time()
            batch = max(1, args.workers)          # shuffles searched at once, spread over all cores
            done = 0
            while done < args.nulls:
                ids = list(range(done, min(args.nulls, done + batch)))
                jobs, owner = [], []
                for i in ids:
                    for lab, seq in shuffled_configs(args.seed + 7919 * i, args.all_configs):
                        jobs.append((lab, seq, g, args.topk))
                        owner.append(i)
                res = pool.starmap(job, jobs, chunksize=1)
                for i in ids:
                    nr = sorted((r for o, rr in zip(owner, res) if o == i for r in rr), reverse=True)[:args.topk]
                    null_ng.append(nr[0][0])
                    null_wd_best.append(nr[0][4])
                    null_wd_max.append(max(r[4] for r in nr))
                done = ids[-1] + 1
                el = time.time() - t1
                log(f"  [{stamp()}] {done}/{args.nulls} shuffles, best n-gram so far "
                    f"{max(null_ng):.3f}, p = {pval(real_ng, null_ng):.4f}  "
                    f"(~{H.fmt_secs(el / done * (args.nulls - done))} left)")

            log(f"\nRESULT for grouping {g}:")
            log(f"  best n-gram     real {real_ng:.3f}   shuffles mean {statistics.mean(null_ng):.3f}, "
                f"max {max(null_ng):.3f}   p = {pval(real_ng, null_ng):.4f}")
            log(f"  words at best   real {real_wd_best}   shuffles mean {statistics.mean(null_wd_best):.1f}, "
                f"max {max(null_wd_best)}   p = {pval(real_wd_best, null_wd_best):.4f}")
            log(f"  max words (top {args.topk})  real {real_wd_max}   shuffles mean "
                f"{statistics.mean(null_wd_max):.1f}, max {max(null_wd_max)}   "
                f"p = {pval(real_wd_max, null_wd_max):.4f}")
            log(f"  (p floor with {args.nulls} shuffles = {1 / (1 + args.nulls):.4f}; "
                "real English of this length scores about -9.2)")

    log("\nHOW TO READ THIS")
    log(" * A readable decrypt with n-gram near -9 to -10 would be a solution. Scores around")
    log("   -13 to -15 are what 10 million tries squeeze out of random symbol order.")
    log(" * Small p on n-gram alone: the symbol order fits these keys better than chance.")
    log("   Small p on words too: stronger. Large p on both: no sets-of-three key reads it as English.")
    log(f"\n[{stamp()}] END allkeys   total {H.fmt_secs(time.time() - t_start)}")
    log(f"log: {os.path.abspath(os.path.join(args.out, 'log.txt'))}   keys: {os.path.abspath(csv_path)}")
    LOG.close()


if __name__ == '__main__':
    mp.freeze_support()
    main()
