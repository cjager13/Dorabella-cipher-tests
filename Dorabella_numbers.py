#!/usr/bin/env python3
"""
Numbers test: what if a few Dorabella symbols are DIGITS, not letters?

Idea (from the 1896 Courage cards, whose first card shows ten loop symbols, as if
Elgar were mapping the digits 0-9 onto his alphabet): if the message contains a date,
time or number, some symbols are digits. Every solver so far forced every symbol to be
a letter, which would spoil the English around each digit.

What it does:
  1. Candidate digit symbols = the rarest symbols (seen at most --max-count times;
     numbers are rare in a short note). Every set of 0..--max-digits of them is tried.
  2. Those symbols are treated as digits: shown as #, and skipped by the English scorer
     (no letter is forced on them).
  3. The remaining symbols are solved as a simple substitution by the same hill climb
     your other scripts use, in each reading order.
  4. The IDENTICAL search is repeated on shuffled copies (shuffling keeps every symbol's
     count, so the same candidates are tried), so the p-value is fair.

  python Dorabella_numbers.py --selftest        # planted English with numbers: can it find them?
  python Dorabella_numbers.py                   # the real Dorabella
  python Dorabella_numbers.py --digit-codes 31 35      # test one specific digit set you suspect
  python Dorabella_numbers.py --max-count 3 --max-digits 4   # wider (slower)

Needs Dorabella_hillclimb.py, Dorabella_tests.py and english.txt in the same folder.
Output: screen + numbers_out/log.txt + numbers_out/results.csv
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
from collections import Counter

import Dorabella_hillclimb as H
import Dorabella_tests as T

NL = H.NL
LOG = None
SELFTEST_TEXT = ("DEAR DORA WE COME TO YOU ON THE 14 BY THE 3 TRAIN AND HOPE THAT YOU "
                 "WILL MEET US AT THE STATION IF IT IS FINE WE SHALL WALK")


def log(msg=''):
    print(msg, flush=True)
    if LOG:
        LOG.write(msg + '\n')
        LOG.flush()


def stamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def build_display_words(corpus, min_count=3):
    """Add every corpus word seen min_count+ times to the word list (display only)."""
    import re
    H.build_words(None)
    with open(corpus, encoding='utf-8', errors='ignore') as fh:
        cnt = Counter(re.findall(r"[A-Z]+", fh.read().upper()))
    extra = {w.translate(H.TOL) for w, n in cnt.items() if n >= min_count and 2 <= len(w) <= H.MAX_WORD}
    extra |= {'WE', 'IS', 'IT', 'TO', 'ON', 'AT', 'BY', 'OF', 'IN', 'AN', 'AS', 'IF', 'OR', 'SO', 'UP', 'US', 'MY', 'ME', 'HE', 'BE', 'DO', 'GO', 'NO'}
    H.WORDS |= extra
    for w in extra:
        for k in range(1, len(w) + 1):
            H.PREFIXES.add(w[:k])
    H.MAXLEN = max(H.MAXLEN, max(len(w) for w in extra))


def init_worker(corpus):
    import signal
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    H.build_scorer(corpus)
    H.build_words(None)


# ------------------------------------------------------------------ solver
SKIP_CHARGE = -10.7


def solve(seed, seq, digits, restarts):
    """Hill climb with `digits` symbols excluded. Returns (score, decrypt with '#').
    Score = (sum of log-probs of the 4-letter windows that avoid digits + SKIP_CHARGE for
    every window a digit breaks) / all windows. The charge stops the search 'buying' a
    higher average by calling awkward letters digits."""
    digits = set(digits)
    letters = sorted(set(c for c in seq if c not in digits))
    idx = {c: i for i, c in enumerate(letters)}
    cidx = [idx.get(c, -1) for c in seq]
    n = len(seq)
    valid = [s for s in range(n - 3) if min(cidx[s:s + 4]) >= 0]
    if len(valid) < 20:
        return -99.0, '#' * n
    get = H.QUAD.get
    F = H.QFLOOR
    nv = len(valid)
    nsym = len(letters)

    def score(key):
        p = [key[c] if c >= 0 else 0 for c in cidx]
        t = 0.0
        for s in valid:
            t += get(((p[s] * NL + p[s + 1]) * NL + p[s + 2]) * NL + p[s + 3], F)
        return t / nv

    rng = random.Random(seed)
    best = (-1e9, None)
    for _ in range(restarts):
        key = list(range(NL))
        rng.shuffle(key)
        cur = score(key)
        improved = True
        while improved:
            improved = False
            for i in range(nsym):
                for j in range(i + 1, NL):
                    key[i], key[j] = key[j], key[i]
                    s = score(key)
                    if s > cur + 1e-12:
                        cur, improved = s, True
                    else:
                        key[i], key[j] = key[j], key[i]
        if cur > best[0]:
            best = (cur, key[:])
    key = best[1]
    text = ''.join(H.LETTERS[key[c]] if c >= 0 else '#' for c in cidx)
    charged = (best[0] * nv + SKIP_CHARGE * (n - 3 - nv)) / (n - 3)
    return charged, text


def job(seed, seq, digits, restarts, label, charge=-10.7):
    global SKIP_CHARGE
    SKIP_CHARGE = charge
    s, text = solve(seed, seq, digits, restarts)
    return s, text, tuple(digits), label


def candidates(seq, max_count, max_digits, fixed):
    if fixed:
        return [tuple(fixed)]
    cnt = Counter(seq)
    rare = sorted(c for c, n in cnt.items() if n <= max_count)
    out = [()]
    for k in range(1, max_digits + 1):
        out += list(itertools.combinations(rare, k))
    return out


def full_search(pool, seq, args, seed, orders):
    cands = candidates(seq, args.max_count, args.max_digits, args.digit_codes)
    jobs = []
    for oi, o in enumerate(orders):
        s = T.apply_order(list(seq), o)
        for ci, d in enumerate(cands):
            jobs.append((seed + 1000 * oi + ci, s, d, args.restarts, o, args.skip_charge))
    res = pool.starmap(job, jobs, chunksize=1) if pool else [job(*j) for j in jobs]
    return res, cands


def words(text):
    return H.word_score(text.replace('#', ' ').translate(H.TOL)) if '#' in text else \
        H.word_score(text.translate(H.TOL))


def seg(text):
    parts = text.split('#')
    return ' # '.join(H.segment(p.translate(H.TOL)) for p in parts)


def show(text):
    if len(text) == 87:
        for a, b in H.LINE_SPANS:
            log("     " + text[a:b])
    else:
        log("     " + text)
    log("     " + seg(text))


def pval(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


# ------------------------------------------------------------------ planted test
def planted(seed):
    rng = random.Random(seed)
    raw = SELFTEST_TEXT.replace(' ', '')
    plain = H.fold_base(raw)[:87]
    letters = sorted(set(c for c in plain if c.isalpha()))
    digits = sorted(set(c for c in plain if c.isdigit()))
    codes = rng.sample(H.CODES_ALL, len(letters) + len(digits))
    code_of = dict(zip(letters + digits, codes))
    return plain, [code_of[c] for c in plain], [code_of[d] for d in digits]


# ------------------------------------------------------------------ main
def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='numbers_out')
    ap.add_argument('--restarts', type=int, default=80)
    ap.add_argument('--nulls', type=int, default=100)
    ap.add_argument('--max-count', type=int, default=2, help='a digit candidate appears at most this often')
    ap.add_argument('--max-digits', type=int, default=3, help='most symbols treated as digits at once')
    ap.add_argument('--digit-codes', type=int, nargs='+', help='test only this digit set, e.g. 31 35')
    ap.add_argument('--orders', nargs='+', default=None,
                    help="reading orders (default all 5): forward 'fully reversed' ...")
    ap.add_argument('--skip-charge', type=float, default=-10.7,
                    help='log-prob charged per 4-letter window broken by a digit (lower = stricter)')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--seed', type=int, default=1896)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    LOG = open(os.path.join(args.out, 'log.txt'), 'a', encoding='utf-8')
    t0 = time.time()
    log(f"\n[{stamp()}] START Dorabella_numbers {' '.join(sys.argv[1:]) or '(defaults)'}")
    log(f"skip charge per broken window: {args.skip_charge}")
    H.build_scorer(args.corpus)
    build_display_words(args.corpus)
    orders = args.orders or T.ORDERS
    pool = mp.Pool(args.workers, initializer=init_worker, initargs=(args.corpus,)) if args.workers > 1 else None

    try:
        if args.selftest:
            plain, seq, true_digits = planted(args.seed)
            orders = ['forward']
            log(f"SELF-TEST planted: {plain}")
            log(f"digit symbols planted: {true_digits}   symbol counts: "
                f"{ {c: n for c, n in Counter(seq).items() if n <= args.max_count} }")
        else:
            seq = list(H.BASE_SEQ)
            cnt = Counter(seq)
            log("Dorabella symbol counts (rarest first): " +
                ', '.join(f"{c}x{n}" for c, n in sorted(cnt.items(), key=lambda x: (x[1], x[0]))))

        cands = candidates(seq, args.max_count, args.max_digits, args.digit_codes)
        log(f"digit sets tried: {len(cands)} (including 'no digits'), reading orders: {len(orders)}, "
            f"restarts: {args.restarts}")
        res, _ = full_search(pool, seq, args, args.seed, orders)
        log(f"[{stamp()}] real search done in {H.fmt_secs(time.time() - t0)}")
        res.sort(key=lambda r: -r[0])

        with open(os.path.join(args.out, 'results.csv'), 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['score', 'words', 'order', 'digit_codes', 'decrypt'])
            for s, text, d, o in res:
                w.writerow([f"{s:.4f}", words(text), o, ' '.join(map(str, d)) or '-', text])

        log("\nTop results (# = a symbol treated as a digit):")
        for s, text, d, o in res[:8]:
            log(f"\n  score {s:.3f}  words {words(text)}  [{o}]  digits = {list(d) or 'none'}")
            show(text)
        base = [r for r in res if not r[2]]
        if base:
            log(f"\nFor comparison, no digits: best {max(r[0] for r in base):.3f}")
        by_k = {}
        for s, text, d, o in res:
            by_k[len(d)] = max(by_k.get(len(d), -99), s)
        log("Best score by number of digit symbols: " +
            ', '.join(f"{k}: {v:.3f}" for k, v in sorted(by_k.items())))

        if args.selftest:
            s, text, d, o = res[0]
            acc = sum(a == b for a, b in zip(text, plain) if b.isalpha()) / sum(c.isalpha() for c in plain)
            log(f"\nSELF-TEST: best digit set {list(d)} (true {true_digits}); "
                f"letters recovered {acc:.0%}")

        real = res[0][0]
        real_w = words(res[0][1])
        log(f"\n[{stamp()}] identical search on {args.nulls} shuffled copies...")
        null_s, null_w = [], []
        t1 = time.time()
        rng = random.Random(args.seed + 99)
        for i in range(args.nulls):
            sh = list(seq)
            rng.shuffle(sh)
            r, _ = full_search(pool, sh, args, args.seed + 7919 * (i + 1), orders)
            b = max(r, key=lambda x: x[0])
            null_s.append(b[0])
            null_w.append(words(b[1]))
            if (i + 1) % 5 == 0 or i + 1 == args.nulls:
                el = time.time() - t1
                log(f"  [{stamp()}] {i + 1}/{args.nulls}: null best so far {max(null_s):.3f}, "
                    f"p = {pval(real, null_s):.4f}  (~{H.fmt_secs(el / (i + 1) * (args.nulls - i - 1))} left)")
        log(f"\nRESULT: best score real {real:.3f}; shuffled mean {statistics.mean(null_s):.3f}, "
            f"max {max(null_s):.3f}; p = {pval(real, null_s):.4f} (floor {1 / (1 + args.nulls):.4f})")
        log(f"        words at best: real {real_w}; shuffled mean {statistics.mean(null_w):.1f}, "
            f"max {max(null_w)}; p = {pval(real_w, null_w):.4f}")
        log("\nNote: on 87 symbols the exact digit symbols are often a near-tie (the top few results")
        log("differ by ~0.001); what matters is whether the LETTERS read as English.")
        log("\nHow to read: a hit is a decrypt that reads as English around the # marks, scores")
        log("well above the shuffled copies, AND whose best digit set clearly beats 'no digits'.")
        log("Removing symbols always helps a little (fewer awkward windows); the shuffled copies")
        log("get the same help, so only the p-value tells you whether it is more than that.")
    except KeyboardInterrupt:
        log("\nstopped by Ctrl+C")
        if pool:
            pool.terminate()
            pool = None
    finally:
        if pool:
            pool.close()
            pool.join()
        log(f"[{stamp()}] END   total {H.fmt_secs(time.time() - t0)}")
        LOG.close()


if __name__ == '__main__':
    mp.freeze_support()
    main()
