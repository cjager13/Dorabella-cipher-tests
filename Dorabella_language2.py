#!/usr/bin/env python3
"""
Dorabella "is it language?" test, corrected and extended (audit of 16 Nov 2026).

Needs Dorabella_language.py, Dorabella_transcription.csv, english.txt and numpy.

  python Dorabella_language2.py                     # about 5-10 minutes
  python Dorabella_language2.py --samples 500000    # tighter p-values

Output: language2_out/report.txt

WHAT WAS WRONG WITH Dorabella_language.py
  1. Its main comparison was "English through Elgar's key", but the key it used was the
     Cipher Mysteries pigpen-style key (BLTACE...), not the key in Elgar's notebook
     (BPECAH...; see Dorabella_notebook.py). This script uses the notebook key.
  2. Its mirror features (F3, F4) only counted left-right and up-down mirror images. The
     published oddity is a different one: neighbours that are the same shape turned half
     round (same loop count, opposite directions, like an E-shape next to a 3-shape).
     Hauer, Choi, Sundar, Hindle, Smallwood and Kondrak (HistoCrypt 2021) count 13 such pairs
     in the Dorabella against 3.6 on average in enciphered English (123 of 100,000 samples had
     13 or more); Massey noticed the same excess. The old test never measured it.

PRE-REGISTERED (frozen 16 Nov 2026)
  Features F1-F9 exactly as in Dorabella_language.py, plus
    F10 half-turn pairs, same loops   neighbours with the same loop count facing opposite ways
    F11 half-turn pairs, any loops    neighbours facing opposite ways, whatever their loops
  Reference sets, all 87 symbols long:
    R1 English (Dracula) through Elgar's NOTEBOOK key      <- the main comparison
    R2 English through a random one-symbol-one-letter key
    R3 symbols drawn at random
    R4 the Dorabella's own symbols in random order
  Tests (three, so the line is p < 0.05/3 = 0.0167):
    P1  the original primary test with the right key: Mahalanobis distance of the Dorabella
        from R1 over F1-F9; p = share of R1 passages at least as far out.
    P2  F10 against R1: share of English passages with at least as many (one-sided, the
        published direction). Its value against R2 is the replication of Hauer et al.
    P3  F11 against R1, one-sided the same way.
  Uncertain symbols: the 9 in the CSV plus position 32 (18 or 11), so 1,024 readings; the
  report gives the median reading, the range, the CSV's main reading and the consensus reading
  (Hauer et al. 2025: 23=14, 32=11, 34=12, 85=14, 86=12). p-values use the median reading.
  HONESTY NOTE: F10 and F11 were counted on the Dorabella during the audit, before this
  script was written. They are a replication of an independent published result, not a new
  discovery of this project, and should be read that way.
"""
import argparse
import datetime
import os
import random

import numpy as np

import Dorabella_language as L

NOTEBOOK_SETS = [('ABC', 3), ('DEF', 4), ('GHI', 5), ('KLM', 6), ('NOP', 7), ('QRS', 8), ('TUW', 1), ('XYZ', 2)]
NB_CODE = {}
for _letters, _d in NOTEBOOK_SETS:
    for _k, _ch in enumerate(_letters):
        NB_CODE[_ch] = (_k + 1) * 10 + _d
NB_CODE['J'] = NB_CODE['I']
NB_CODE['V'] = NB_CODE['U']
NB_LETTER = {v: k for k, v in NB_CODE.items() if k not in 'JV'}
CONSENSUS = {23: 14, 32: 11, 34: 12, 85: 14, 86: 12}
EXTRA_ALT = {32: 11}                         # 1-based position: alternate not in the CSV
FEATURES = L.FEATURES + ['F10 half-turn pairs, same loops', 'F11 half-turn pairs, any loops']
NF = len(FEATURES)
COUNTS = (1, 2, 3, 8, 9, 10)                 # features printed as whole numbers
BONF = 0.05 / 3


def features(seq):
    f = L.features(seq)
    same = anyl = 0
    for a, b in zip(seq, seq[1:]):
        if (a % 10 - b % 10) % 8 == 4:
            anyl += 1
            same += a // 10 == b // 10
    return f + [same, anyl]


def make_r1(letters, n, rng):
    return [[NB_CODE[ch] for ch in w] for w in L.english_windows(letters, n, rng)]


def fmt(v, f):
    return f"{v:.0f}" if f in COUNTS else f"{v:.3f}"


def upper_p(x, ref):
    ref = np.asarray(ref)
    return (1 + np.sum(ref >= x - 1e-12)) / (1 + len(ref))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--samples', type=int, default=200000)
    ap.add_argument('--out', default='language2_out')
    ap.add_argument('--seed', type=int, default=1897)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    lines = []

    def log(msg=''):
        print(msg, flush=True)
        lines.append(msg)

    t0 = datetime.datetime.now()
    log("DORABELLA 'IS IT LANGUAGE?' REPORT 2  (notebook key; half-turn pairs added)")
    log(f"written {t0:%Y-%m-%d %H:%M:%S}, {args.samples:,} passages per reference set, seed {args.seed}")
    seq, alts = L.load_cipher()
    for pos, code in EXTRA_ALT.items():
        alts.setdefault(pos - 1, code)
    cons = list(seq)
    for pos, code in CONSENSUS.items():
        cons[pos - 1] = code
    log(f"check: the notebook key on the consensus reading gives {''.join(NB_LETTER[c] for c in cons)[:30]}... "
        f"(widely printed: BPECAHTCKYFRQDRIRRHPPRDXYXGFS)")

    rng = random.Random(args.seed)
    letters = L.load_letters(args.corpus)
    sets = {
        'R1 English, notebook key': np.array([features(s) for s in make_r1(letters, args.samples, rng)]),
        'R2 English, random key': np.array([features(s) for s in L.make_r2(letters, args.samples, rng)]),
        'R3 random symbols': np.array([features(s) for s in L.make_r3(args.samples, rng)]),
        'R4 Dorabella shuffled': np.array([features(s) for s in L.make_r4(seq, args.samples, rng)]),
    }
    readings = np.array([features(s) for s in L.all_readings(seq, alts)])
    dora = np.median(readings, axis=0)
    base = np.array(features(seq))
    fc = np.array(features(cons))
    log(f"Dorabella: {len(readings)} readings of the {len(alts)} uncertain symbols measured")
    R1 = sets['R1 English, notebook key']
    R2 = sets['R2 English, random key']

    log("\n" + "=" * 110)
    log(f"HEADLINE  (three pre-registered tests; line to beat p < {BONF:.4f})")
    log("=" * 110)
    dx, p1, d_ref = L.mahalanobis_p(dora[:9], R1[:, :9])
    per = [L.mahalanobis_p(r[:9], R1[:, :9])[1] for r in readings[::32]]
    pc1 = L.mahalanobis_p(fc[:9], R1[:, :9])[1]
    log(f"P1 all of F1-F9 against English through the notebook key: distance {dx:.1f} "
        f"(typical English {np.median(d_ref):.1f}, 99th pct {np.percentile(d_ref, 99):.1f})")
    log(f"     p = {p1:.5f}   (range over readings {min(per):.5f} to {max(per):.5f}; consensus reading {pc1:.5f})")
    for tag, f in (('P2', 9), ('P3', 10)):
        log(f"{tag} {FEATURES[f]}: Dorabella {dora[f]:.0f} (readings {readings[:, f].min():.0f}-"
            f"{readings[:, f].max():.0f}; CSV main {base[f]:.0f}; consensus {fc[f]:.0f})")
        log(f"     English, notebook key: mean {R1[:, f].mean():.2f}, p (as many or more) = {upper_p(dora[f], R1[:, f]):.5f}"
            f"   | consensus reading p = {upper_p(fc[f], R1[:, f]):.5f}")
        log(f"     English, random key:   mean {R2[:, f].mean():.2f}, p = {upper_p(dora[f], R2[:, f]):.5f}"
            + ("   <- the Hauer et al. (2021) comparison: they report 13 against a mean of 3.64, p about 0.001"
               if f == 9 else ''))
        log(f"     its own symbols shuffled: mean {sets['R4 Dorabella shuffled'][:, f].mean():.2f}, "
            f"p = {upper_p(dora[f], sets['R4 Dorabella shuffled'][:, f]):.5f}")
    log("Reading it: a small P1 says the cipher's visual habits as a whole are not those of English")
    log("through Elgar's own key. A small P2 or P3 says its neighbours face opposite ways more often")
    log("than any enciphered English does: a habit of the hand, or a rule of the cipher that these")
    log("tests do not model. Neither is a decryption, and neither shows there is no message.")

    log("\n" + "=" * 110)
    log("FEATURE BY FEATURE  (Dorabella = median over readings [range];  each set: mean [2.5%-97.5%];")
    log("two-sided p against R1 and R2, for description only)")
    log("=" * 110)
    for f in range(NF):
        lo, hi = readings[:, f].min(), readings[:, f].max()
        log(f"\n{FEATURES[f]}:  Dorabella {fmt(dora[f], f)}  [{fmt(lo, f)}-{fmt(hi, f)}]   "
            f"(CSV main {fmt(base[f], f)}, consensus {fmt(fc[f], f)})")
        for name, arr in sets.items():
            col = arr[:, f]
            extra = ''
            if name.startswith('R1') or name.startswith('R2'):
                pv, side = L.two_sided(dora[f], col)
                extra = f"   p = {pv:.5f} (Dorabella {side})"
            log(f"    {name:<26} {col.mean():.3f} [{fmt(np.percentile(col, 2.5), f)}-"
                f"{fmt(np.percentile(col, 97.5), f)}]{extra}")

    log("\n" + "=" * 110)
    log("WHERE THE HALF-TURN PAIRS ARE (CSV main reading; positions of the first symbol of each pair)")
    log("=" * 110)
    for name, a, b in (('line 1', 0, 29), ('line 2', 29, 60), ('line 3', 60, 87)):
        pos = [i + 1 for i in range(a, b - 1) if (seq[i] % 10 - seq[i + 1] % 10) % 8 == 4]
        same = [i for i in pos if seq[i - 1] // 10 == seq[i] // 10]
        log(f"  {name}: {len(pos)} of {b - a - 1} neighbours face opposite ways (same loops at {same}); all at {pos}")
    log(f"\nfinished {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    with open(os.path.join(args.out, 'report.txt'), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
