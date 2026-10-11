#!/usr/bin/env python3
"""
Dorabella "is it language?" test: were the symbols chosen BY A MESSAGE, or BY THE HAND?

Needs Dorabella_transcription.csv and english.txt in the same folder, plus numpy.

  python Dorabella_language.py                        # about 5-10 minutes
  python Dorabella_language.py --samples 500000       # tighter p-values
  python Dorabella_language.py --human Dorabella_human_samples.txt
                                                      # add hand-drawn "fake ciphers"

Output: language_out/report.txt (and log.txt)

THE IDEA
A key scatters a message's letters over the 8 directions and 3 loop counts with no
regard for how they LOOK next to each other. Someone drawing symbols to look pleasing
leaves visual habits: small turns, mirror pairs, alternating loop counts, no doubles.
So we measure the visual habits of the Dorabella and ask where it falls among
thousands of real English passages put through a key.

IMPORTANT: scrambled copies cannot be the comparison here (pair-preserving shuffles
keep exactly the pair habits being measured). The comparison is REAL ENGLISH THROUGH A KEY.

PRE-REGISTERED FEATURES (adjacent symbol pairs, written order; frozen before any run)
  F1 same-loop rate        share of neighbours with the same loop count
  F2 longest loop-change   longest run in which every neighbour has a different loop
                           count (Massey's "twelve in a row")
  F3 left-right mirrors    neighbours with the same loop count that are left-right
                           mirror images of each other (e.g. opening right, then left)
  F4 up-down mirrors       the same for up-down mirror images
  F5 mean turn             average turn between neighbours' openings (0 to 4 eighths)
  F6 45-degree turns       share of neighbours whose openings differ by one eighth
  F7 doubled symbols       share of neighbours that are the identical symbol
                           (English has LL, EE, SS, OO...)
  F8 index of coincidence  how unevenly the 24 symbols are used (key-independent)
  F9 repeated pairs        how often adjacent symbol pairs repeat (key-independent)

REFERENCE SETS (all 87 symbols long)
  R1 English (Dracula) through ELGAR'S OWN KEY        <- the main comparison
  R2 English through a random one-symbol-one-letter key (any key, any language-like text)
  R3 symbols drawn at random, all 24 equally likely    (no message, no habits)
  R4 the Dorabella's own symbols in random order        (same symbols, no order)
  H  optional: hand-drawn fake ciphers from volunteers (--human)

TESTS
  Primary (one test): the Dorabella's Mahalanobis distance from the R1 cloud over all
  9 features; p = share of R1 passages at least as far out.
  Secondary: each feature vs R1 and R2, two-sided; Bonferroni line 0.05/9 = 0.0056.
  The 9 uncertain symbols: every one of the 512 readings is measured; the report gives
  the median reading and the range, so no single reading decides the answer.
"""
import argparse
import csv
import datetime
import itertools
import os
import random
import re
import statistics

import numpy as np

FEATURES = ['F1 same-loop rate', 'F2 longest loop-change run', 'F3 left-right mirror pairs',
            'F4 up-down mirror pairs', 'F5 mean turn (eighths)', 'F6 45-degree turn rate',
            'F7 doubled-symbol rate', 'F8 index of coincidence', 'F9 repeated-pair excess']
NF = len(FEATURES)
BONF = 0.05 / NF

# Elgar's published key (my direction 1 = opening up). '^': 1 loop = first letter of the
# set; 'v': 1 loop = third letter. Reproduces the published decrypt (BLTACEIARWUN...).
ELGAR_SETS = [('ABC', 3, 'v'), ('DEF', 5, '^'), ('GHI', 1, 'v'), ('JKL', 7, '^'),
              ('MNO', 8, '^'), ('PQR', 6, 'v'), ('STU', 4, '^'), ('WXY', 2, 'v')]
ELGAR_CODE = {}
for letters, d, ud in ELGAR_SETS:
    for loop in (1, 2, 3):
        ch = letters[loop - 1] if ud == '^' else letters[3 - loop]
        ELGAR_CODE[ch] = loop * 10 + d
ELGAR_CODE['V'] = ELGAR_CODE['U']        # Elgar's key has no V or Z
ELGAR_CODE['Z'] = ELGAR_CODE['S']
ALL_CODES = [loop * 10 + d for loop in (1, 2, 3) for d in range(1, 9)]


# ============================================================== features
def features(seq):
    loops = [c // 10 for c in seq]
    dirs = [c % 10 - 1 for c in seq]          # 0..7, 0 = opening up
    n = len(seq)
    pairs = list(zip(range(n - 1), range(1, n)))
    same_loop = sum(loops[i] == loops[j] for i, j in pairs) / (n - 1)
    run = best = 1
    for i, j in pairs:
        run = run + 1 if loops[i] != loops[j] else 1
        best = max(best, run)
    lr = ud = 0
    turns = []
    for i, j in pairs:
        a, b = dirs[i], dirs[j]
        if loops[i] == loops[j] and a != b:
            if (-a) % 8 == b:
                lr += 1
            if (4 - a) % 8 == b:
                ud += 1
        t = abs(a - b) % 8
        turns.append(min(t, 8 - t))
    mean_turn = sum(turns) / len(turns)
    t45 = sum(t == 1 for t in turns) / len(turns)
    doubled = sum(seq[i] == seq[j] for i, j in pairs) / (n - 1)
    cnt = {}
    for c in seq:
        cnt[c] = cnt.get(c, 0) + 1
    ic = sum(v * (v - 1) for v in cnt.values()) / (n * (n - 1))
    bg = {}
    for i, j in pairs:
        k = (seq[i], seq[j])
        bg[k] = bg.get(k, 0) + 1
    rep = sum(v - 1 for v in bg.values())
    return [same_loop, best, lr, ud, mean_turn, t45, doubled, ic, rep]


# ============================================================== reference sets
def load_letters(corpus):
    with open(corpus, encoding='utf-8', errors='ignore') as fh:
        raw = fh.read().upper()
    return ''.join(re.findall(r'[A-Z]', raw))


def english_windows(letters, n, rng, length=87):
    for _ in range(n):
        s = rng.randrange(len(letters) - length)
        yield letters[s:s + length]


def make_r1(letters, n, rng):
    return [[ELGAR_CODE[ch] for ch in w] for w in english_windows(letters, n, rng)]


def make_r2(letters, n, rng):
    out = []
    for w in english_windows(letters, n, rng):
        w = w.replace('J', 'I').replace('V', 'U')
        key = dict(zip(sorted(set('ABCDEFGHIKLMNOPQRSTUWXYZ')), rng.sample(ALL_CODES, 24)))
        out.append([key[ch] for ch in w])
    return out


def make_r3(n, rng):
    return [[rng.choice(ALL_CODES) for _ in range(87)] for _ in range(n)]


def make_r4(seq, n, rng):
    out = []
    for _ in range(n):
        s = list(seq)
        rng.shuffle(s)
        out.append(s)
    return out


def load_human(path):
    """One sample per line: symbol codes (loops then direction, e.g. 23 37 24 ...),
    separated by spaces or commas. Lines starting with # are ignored."""
    out = []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            codes = [int(x) for x in re.split(r'[\s,;]+', line) if x]
            bad = [c for c in codes if c not in ALL_CODES]
            if bad:
                raise SystemExit(f"bad symbol codes {bad} in {path}")
            if len(codes) >= 20:
                out.append(codes)
    return out


# ============================================================== cipher
def load_cipher(path='Dorabella_transcription.csv'):
    with open(path, newline='', encoding='utf-8') as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r['position']))
    seq = [int(r['code']) for r in rows]
    alts = {int(r['position']) - 1: int(r['alternate_code'])
            for r in rows if r['sure'].strip().lower() == 'no' and r['alternate_code'].strip()}
    return seq, alts


def all_readings(seq, alts):
    pos = sorted(alts)
    for bits in itertools.product((0, 1), repeat=len(pos)):
        s = list(seq)
        for p, b in zip(pos, bits):
            if b:
                s[p] = alts[p]
        yield s


# ============================================================== stats
def two_sided(x, ref):
    ref = np.asarray(ref)
    lo = (1 + np.sum(ref <= x + 1e-12)) / (1 + len(ref))
    hi = (1 + np.sum(ref >= x - 1e-12)) / (1 + len(ref))
    return min(1.0, 2 * min(lo, hi)), ('low' if lo < hi else 'high')


def mahalanobis_p(x, ref):
    R = np.asarray(ref, dtype=float)
    mu = R.mean(axis=0)
    cov = np.cov(R, rowvar=False)
    inv = np.linalg.pinv(cov)
    d_ref = np.einsum('ij,jk,ik->i', R - mu, inv, R - mu)
    dx = float((np.asarray(x) - mu) @ inv @ (np.asarray(x) - mu))
    return dx, (1 + np.sum(d_ref >= dx)) / (1 + len(R)), d_ref


def fmt(v, f):
    return f"{v:.0f}" if f in (1, 2, 3, 8) else f"{v:.3f}"


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--samples', type=int, default=200000, help='passages per reference set')
    ap.add_argument('--human', help='file of hand-drawn fake ciphers, one per line')
    ap.add_argument('--out', default='language_out')
    ap.add_argument('--seed', type=int, default=1897)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    lines = []

    def log(msg=''):
        print(msg, flush=True)
        lines.append(msg)

    t0 = datetime.datetime.now()
    log(f"DORABELLA 'IS IT LANGUAGE?' REPORT\nwritten {t0:%Y-%m-%d %H:%M:%S}, "
        f"{args.samples:,} passages per reference set, seed {args.seed}")
    seq, alts = load_cipher()
    pub = ''.join({v: k for k, v in ELGAR_CODE.items() if k not in 'VZ'}[c] for c in seq)
    log(f"check: Elgar's key on the transcription gives {pub[:30]}... (published: BLTACEIARWUNISNF...)")

    rng = random.Random(args.seed)
    letters = load_letters(args.corpus)
    sets = {
        'R1 English, Elgar key': np.array([features(s) for s in make_r1(letters, args.samples, rng)]),
        'R2 English, random key': np.array([features(s) for s in make_r2(letters, args.samples, rng)]),
        'R3 random symbols': np.array([features(s) for s in make_r3(args.samples, rng)]),
        'R4 Dorabella shuffled': np.array([features(s) for s in make_r4(seq, args.samples, rng)]),
    }
    if args.human:
        hs = load_human(args.human)
        log(f"human samples: {len(hs)} from {args.human}")
        if hs:
            sets['H hand-drawn'] = np.array([features(s) for s in hs])

    readings = np.array([features(s) for s in all_readings(seq, alts)])
    dora = np.median(readings, axis=0)
    base = np.array(features(seq))
    log(f"Dorabella: {len(readings)} readings of the {len(alts)} uncertain symbols measured; "
        f"'median' below is the median over all of them")

    log("\n" + "=" * 110)
    log("PRIMARY TEST: is the Dorabella an unusual English passage under Elgar's key? (all 9 features at once)")
    log("=" * 110)
    R1 = sets['R1 English, Elgar key']
    dx, p, d_ref = mahalanobis_p(dora, R1)
    per_read = [mahalanobis_p(r, R1)[1] for r in readings[::16]]
    log(f"distance of the Dorabella from the English-through-Elgar's-key cloud: {dx:.1f} "
        f"(typical English passage: median {np.median(d_ref):.1f}, 99th pct {np.percentile(d_ref, 99):.1f})")
    log(f"p = {p:.5f}  (share of real English passages at least this far out; floor {1 / (1 + len(R1)):.6f})")
    log(f"over the uncertain-symbol readings, p ranges {min(per_read):.5f} to {max(per_read):.5f}")
    for name in ('R2 English, random key', 'R3 random symbols', 'R4 Dorabella shuffled', 'H hand-drawn'):
        if name in sets:
            d2, p2, _ = mahalanobis_p(dora, sets[name])
            log(f"  same distance test against {name:<24}: p = {p2:.5f}")

    log("\n" + "=" * 110)
    log("FEATURE BY FEATURE  (Dorabella = median over readings [range];  each set: mean [2.5%-97.5%])")
    log(f"two-sided p vs R1 and R2; Bonferroni line {BONF:.4f}")
    log("=" * 110)
    for f in range(NF):
        lo, hi = readings[:, f].min(), readings[:, f].max()
        log(f"\n{FEATURES[f]}:  Dorabella {fmt(dora[f], f)}  [{fmt(lo, f)}-{fmt(hi, f)}]   "
            f"(as transcribed {fmt(base[f], f)})")
        for name, arr in sets.items():
            col = arr[:, f]
            extra = ''
            if name.startswith('R1') or name.startswith('R2'):
                pv, side = two_sided(dora[f], col)
                flag = '  <-- past Bonferroni' if pv < BONF else ('  (p<0.05)' if pv < 0.05 else '')
                extra = f"   p = {pv:.5f} ({'Dorabella ' + side}){flag}"
            log(f"    {name:<24} {fmt(col.mean(), 4 if f in (1, 2, 3, 8) else f)} "
                f"[{fmt(np.percentile(col, 2.5), f)}-{fmt(np.percentile(col, 97.5), f)}]{extra}")

    log("\n" + "=" * 110)
    log("HOW TO READ IT")
    log("=" * 110)
    log("* Primary p tiny (below 0.001): the Dorabella's visual habits are not what English through")
    log("  Elgar's key produces. Then look at WHICH features are off, and whether the Dorabella")
    log("  sits nearer the hand-drawn samples (H) or random symbols (R3).")
    log("* Primary p ordinary: the Dorabella looks like a plausible enciphered English passage, so")
    log("  its visual habits give no reason to doubt it carries a message.")
    log("* R2 matters for features that depend on the key's layout (F1-F7): Elgar's key places")
    log("  letter groups by direction, so R1 and R2 can differ; F8-F9 do not depend on the key.")
    log(f"\nfinished {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    with open(os.path.join(args.out, 'report.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    with open(os.path.join(args.out, 'log.txt'), 'a', encoding='utf-8') as fh:
        fh.write(f"{t0:%Y-%m-%d %H:%M:%S} run, samples {args.samples}, human {args.human}\n")


if __name__ == '__main__':
    main()
