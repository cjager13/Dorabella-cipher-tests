#!/usr/bin/env python3
"""
Loop-group keys: the loop count picks the letter group and the direction picks the
letter within it (e.g. 1 loop = A-H, 2 loops = I/J-Q, 3 loops = R-Z).

Every such key is tried: 6 ways to give the groups to the loop counts x each group read
forwards or backwards (2^3) x each group rotated to any of 8 starting directions (8^3)
= 24,576 keys, in all 160 reading configurations (5 reading orders x 32 alternate
readings). The identical search is repeated on plain shuffles and on pair-preserving
shuffles of the cipher.

Needs numpy, Dorabella_hillclimb.py, Dorabella_tests.py, Dorabella_marathon.py and
english.txt in the same folder.

  python Dorabella_loopgroups.py                 # 200 copies of each kind
  python Dorabella_loopgroups.py --nulls 50      # quicker

Output: screen + loopgroups_out/log.txt
"""
import argparse
import datetime
import itertools
import os
import statistics
import time

import numpy as np

import Dorabella_hillclimb as H
import Dorabella_marathon as M
import Dorabella_tests as T

GROUPS = ['ABCDEFGH', 'IKLMNOPQ', 'RSTUWXYZ']   # I = I/J, U = U/V
LOG = None


def log(msg=''):
    print(msg, flush=True)
    LOG.write(msg + '\n')
    LOG.flush()


def build_keys():
    keys, rows = [], []
    for perm in itertools.permutations(range(3)):
        for rev in itertools.product((0, 1), repeat=3):
            for rots in itertools.product(range(8), repeat=3):
                row = [0] * 24
                for loop in range(3):
                    g = list(GROUPS[perm[loop]])
                    g = g[::-1] if rev[loop] else g
                    g = g[rots[loop]:] + g[:rots[loop]]
                    for d in range(8):
                        row[loop * 8 + d] = H.LIDX[g[d]]
                rows.append(row)
                keys.append((perm, rev, rots))
    return keys, np.array(rows, dtype=np.int64)


def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--nulls', type=int, default=200)
    ap.add_argument('--seed', type=int, default=1000)
    args = ap.parse_args()
    os.makedirs('loopgroups_out', exist_ok=True)
    LOG = open(os.path.join('loopgroups_out', 'log.txt'), 'a', encoding='utf-8')
    log(f"\n[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] START Dorabella_loopgroups  nulls={args.nulls}")

    H.build_scorer(args.corpus)
    H.build_words(None)
    QT = np.full(24 ** 4, H.QFLOOR)
    for k, v in H.QUAD.items():
        QT[k] = v
    keys, KL = build_keys()
    log(f"keys: {len(KL):,}   reading configurations: 160")

    def ci(seq):
        return np.array([(c // 10 - 1) * 8 + (c % 10 - 1) for c in seq])

    def best_over(seqs):
        best = (-1e9, None)
        for lab, seq in seqs:
            P = KL[:, ci(seq)]
            q = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
            S = QT[q].mean(1)
            i = int(S.argmax())
            if S[i] > best[0]:
                best = (float(S[i]), (lab, i, seq))
        return best

    def text_of(i, seq):
        return ''.join(H.LETTERS[x] for x in KL[i][ci(seq)])

    toks = M.make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})

    def configs(t):
        return [(f"{o} | readings {M.BITS[vi]}", T.apply_order(M.realize(t, M.BITS[vi]), o))
                for vi in range(len(M.BITS)) for o in T.ORDERS]

    t0 = time.time()
    real, (lab, i, seq) = best_over(configs(toks))
    text = text_of(i, seq)
    rw = H.word_score(text.translate(H.TOL))
    perm, rev, rots = keys[i]
    groups = ', '.join(f"{n + 1} loop(s) = {GROUPS[perm[n]]}{' reversed' if rev[n] else ''} rot {rots[n]}"
                       for n in range(3))
    log(f"real best {real:.3f} [{lab}]  words {rw}  ({time.time() - t0:.1f}s)")
    log(f"key: {groups}")
    for a, b in H.LINE_SPANS:
        log("   " + text[a:b])
    log("   " + H.segment(text.translate(H.TOL)))

    for kind in ('shuffle', 'doublet'):
        ns, nw = [], []
        t1 = time.time()
        for k in range(args.nulls):
            t = M.null_tokens(kind, toks, args.seed + k)
            s, (_, ii, sq) = best_over(configs(t))
            ns.append(s)
            nw.append(H.word_score(text_of(ii, sq).translate(H.TOL)))
        p = (1 + sum(x >= real for x in ns)) / (1 + len(ns))
        pw = (1 + sum(x >= rw for x in nw)) / (1 + len(nw))
        log(f"vs {kind} ({len(ns)} copies): mean {statistics.mean(ns):.3f}, max {max(ns):.3f}, "
            f"p = {p:.4f} | words mean {statistics.mean(nw):.1f}, p = {pw:.3f}  "
            f"({H.fmt_secs(time.time() - t1)})")
    log("Real English of this length scores about -9.2.")
    LOG.close()


if __name__ == '__main__':
    main()
