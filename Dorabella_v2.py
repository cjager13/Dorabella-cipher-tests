#!/usr/bin/env python3
"""
Dorabella v2: your rechecked transcription, with the 10 "on the fence" symbols
solved as part of the key.

The hill climb (the same kind your other scripts use) now has two kinds of move:
  * swap two letters in the key            (as before)
  * flip a fence symbol to its other reading  (new)
so it picks the best letters AND the best reading of each fence symbol together.
Each restart is also 'kicked' (shaken and re-climbed) 10 times: on planted English
this recovers 97-99% of letters, versus 69-98% for plain restarts.
All 2^10 = 1,024 reading combinations are covered without trying them one by one.
Every reading order is searched.

Fair test: the IDENTICAL search is repeated on scrambled copies of the cipher, and
fence symbols keep their two readings when scrambled, so the copies get exactly the
same freedom. Two kinds of copy:
  shuffle - symbols in random order
  doublet - random order that keeps every adjacent symbol pair (the tougher test)

  python dora_v2.py --selftest            # planted English with 10 fence symbols
  python dora_v2.py                       # the real Dorabella (v2 transcription)
  python dora_v2.py --nulls 200 --kind shuffle
  python dora_v2.py --lang de --selftest   # German (needs german.txt)
  python dora_v2.py --lang de              # German: the real Dorabella

Needs Dorabella_hillclimb.py, Dorabella_tests.py and english.txt in the same folder.
Output: screen + v2_out/log.txt + v2_out/nulls.csv
"""
import argparse
import csv
import datetime
import multiprocessing as mp
import os
import random
import re
import signal
import statistics
import sys
import time
from collections import Counter, defaultdict

import Dorabella_hillclimb as H
import Dorabella_tests as T

NL = H.NL

# ---------------------------------------------------------------- v2 transcription
# From your sheet (SortDirection tab), 2026-09-29. Your direction numbering.
LINE1 = [23,37,24,33,13,25,11,33,16,22,34,28,18,14,28,35,28,28,25,37,37,28,15,12,22,12,15,34,38]
LINE2 = [11,28,11,25,11,33,16,26,33,24,28,28,24,25,15,18,11,28,34,28,25,21,38,18,14,12,16,16,12,34,38]
LINE3 = [24,38,25,21,38,24,14,21,31,15,38,24,28,25,21,18,38,15,33,37,15,38,25,33,15,13,33]
SEQ = LINE1 + LINE2 + LINE3
# 1-based position: alternate reading (the main reading is the one in SEQ)
FENCE = {10: 23, 23: 14, 24: 11, 25: 23, 32: 18, 34: 12, 67: 15, 78: 14, 85: 14, 86: 12}
assert len(SEQ) == 87

SELFTEST_TEXT = {
    'en': ("MY DEAR DORA I WAS SO GLAD OF YOUR LETTER AND THE NEWS THAT YOU ARE ALL WELL "
           "AT HOME WE SHALL COME OVER ON THURSDAY"),
    'de': ("LIEBE DORA ICH DANKE DIR HERZLICH FUER DEINEN LIEBEN BRIEF WIR KOMMEN AM "
           "DONNERSTAG MIT DEM FRUEHEN ZUGE UND FREUEN UNS SEHR AUF DAS WIEDERSEHEN"),
}
UMLAUTS = {'Ä': 'AE', 'Ö': 'OE', 'Ü': 'UE', 'ß': 'SS', 'ẞ': 'SS', '\u0364': 'E',
           'À': 'A', 'Á': 'A', 'Â': 'A', 'È': 'E', 'É': 'E', 'Ê': 'E', 'Ë': 'E',
           'Î': 'I', 'Ï': 'I', 'Ô': 'O', 'Û': 'U', 'Ç': 'C', 'Ñ': 'N'}


def prepare_corpus(path, lang, out_dir):
    """English: used as-is (same as all your other scripts). German: strip Project
    Gutenberg's English licence text and spell out umlauts (ä->AE, ö->OE, ü->UE, ß->SS),
    then save a cleaned copy that the scorer reads."""
    if lang == 'en':
        return path
    with open(path, encoding='utf-8', errors='ignore') as fh:
        raw = fh.read()
    parts = re.split(r'\*\*\* ?START OF[^\n]*\n', raw)
    if len(parts) > 1:
        raw = '\n'.join(re.split(r'\*\*\* ?END OF', p)[0] for p in parts[1:])
    up = raw.upper()
    text = ''.join(UMLAUTS.get(ch, ch) for ch in up)
    clean = os.path.join(out_dir, 'corpus_' + lang + '_clean.txt')
    with open(clean, 'w', encoding='utf-8') as fh:
        fh.write(text)
    return clean

LOG = None


def log(msg=''):
    print(msg, flush=True)
    if LOG:
        LOG.write(msg + '\n')
        LOG.flush()


def stamp():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def build_display_words(corpus, min_count=3, lang='en'):
    """Word list for the reading aid: corpus words seen 3+ times (+ the embedded English
    list when the language is English)."""
    H.build_words(None)
    if lang != 'en':
        H.WORDS, H.PREFIXES, H.MAXLEN = set(), set(), 1
    with open(corpus, encoding='utf-8', errors='ignore') as fh:
        cnt = Counter(re.findall(r"[A-Z]+", fh.read().upper()))
    extra = {w.translate(H.TOL) for w, n in cnt.items() if n >= min_count and 2 <= len(w) <= H.MAX_WORD}
    H.WORDS |= extra
    for w in extra:
        for k in range(1, len(w) + 1):
            H.PREFIXES.add(w[:k])
    H.MAXLEN = max(H.MAXLEN, max(len(w) for w in extra))


def init_worker(corpus, lang='en'):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    H.build_scorer(corpus)
    build_display_words(corpus, lang=lang)


# ---------------------------------------------------------------- tokens & nulls
def base_tokens(seq=SEQ, fence=FENCE):
    """token = (main code, alternate code or None)"""
    return [(c, fence.get(i + 1)) for i, c in enumerate(seq)]


def doublet_shuffle(seq, rng):
    """Random order keeping every adjacent pair count (random Eulerian path)."""
    n = len(seq)
    out_edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        out_edges[a].append(b)
    last = seq[-1]
    in_tree, nxt = {last}, {}
    for v in list(out_edges):
        if v == last:
            continue
        u = v
        while u not in in_tree:
            k = rng.randrange(len(out_edges[u]))
            nxt[u] = k
            u = out_edges[u][k]
        u = v
        while u not in in_tree:
            in_tree.add(u)
            u = out_edges[u][nxt[u]]
    lists = {}
    for v, es in out_edges.items():
        if v == last:
            l = list(es)
            rng.shuffle(l)
        else:
            k = nxt[v]
            l = es[:k] + es[k + 1:]
            rng.shuffle(l)
            l.append(es[k])
        lists[v] = l
    ptr = defaultdict(int)
    out, u = [seq[0]], seq[0]
    for _ in range(n - 1):
        w = lists[u][ptr[u]]
        ptr[u] += 1
        out.append(w)
        u = w
    assert Counter(zip(out, out[1:])) == Counter(zip(seq, seq[1:]))
    return out


def null_tokens(kind, toks, seed):
    rng = random.Random(seed)
    if kind == 'shuffle':
        t = list(toks)
        rng.shuffle(t)
        return t
    new = doublet_shuffle([m for m, _ in toks], rng)
    fence_by_code = defaultdict(list)
    for tk in toks:
        if tk[1] is not None:
            fence_by_code[tk[0]].append(tk)
    out = [(c, None) for c in new]
    for c, lst in fence_by_code.items():
        slots = [i for i, x in enumerate(new) if x == c]
        for i, tk in zip(rng.sample(slots, len(lst)), lst):
            out[i] = tk
    return out


# ---------------------------------------------------------------- joint climb
def climb_joint(toks, restarts, seed, kicks=0):
    """Best (score, letters text, chosen codes) over key swaps + fence flips.
    Each restart climbs to a peak, then gets `kicks` shakes (3 random swaps/flips followed
    by a new climb, kept only if it ends higher): iterated local search."""
    codes = sorted({m for m, _ in toks} | {a for _, a in toks if a is not None})
    idx = {c: i for i, c in enumerate(codes)}
    nsym = len(codes)
    cm = [idx[m] for m, _ in toks]
    ca = [idx[a] if a is not None else -1 for _, a in toks]
    fences = [i for i, a in enumerate(ca) if a >= 0]
    score_fn = H.SCORE_FN
    rng = random.Random(seed)

    def settle(key, ch):
        cur_idx = [ca[i] if ch.get(i) else cm[i] for i in range(len(toks))]
        cur = score_fn([key[c] for c in cur_idx])
        improved = True
        while improved:
            improved = False
            for i in range(nsym):
                for j in range(i + 1, NL):
                    key[i], key[j] = key[j], key[i]
                    s = score_fn([key[c] for c in cur_idx])
                    if s > cur + 1e-12:
                        cur, improved = s, True
                    else:
                        key[i], key[j] = key[j], key[i]
            for f in fences:
                ch[f] = not ch[f]
                cur_idx[f] = ca[f] if ch[f] else cm[f]
                s = score_fn([key[c] for c in cur_idx])
                if s > cur + 1e-12:
                    cur, improved = s, True
                else:
                    ch[f] = not ch[f]
                    cur_idx[f] = ca[f] if ch[f] else cm[f]
        return cur, cur_idx

    best = (-1e9, None, None)
    for _ in range(restarts):
        key = list(range(NL))
        rng.shuffle(key)
        ch = {f: rng.random() < 0.5 for f in fences}
        cur, cur_idx = settle(key, ch)
        for _ in range(kicks):
            k2, c2 = key[:], dict(ch)
            for _ in range(3):
                if fences and rng.random() < 0.25:
                    f = rng.choice(fences)
                    c2[f] = not c2[f]
                else:
                    i, j = rng.randrange(nsym), rng.randrange(NL)
                    k2[i], k2[j] = k2[j], k2[i]
            s2, idx2 = settle(k2, c2)
            if s2 > cur + 1e-12:
                key, ch, cur, cur_idx = k2, c2, s2, idx2
        if cur > best[0]:
            text = ''.join(H.LETTERS[key[c]] for c in cur_idx)
            chosen = [codes[c] for c in cur_idx]
            best = (cur, text, chosen)
    return best


def order_job(toks, order, restarts, seed, kicks=10):
    t = T.apply_order(list(toks), order)
    s, text, chosen = climb_joint(t, restarts, seed, kicks)
    return order, s, text, chosen, list(t)


def full_search_serial(toks, restarts, seed, kicks=10):
    """All reading orders, in one process (used for null copies)."""
    best = None
    for k, o in enumerate(T.ORDERS):
        r = order_job(toks, o, restarts, seed + 101 * k, kicks)
        if best is None or r[1] > best[1]:
            best = r
    return best


def null_job(kind, toks, restarts, seed, kicks=10):
    t = null_tokens(kind, toks, seed)
    o, s, text, chosen, _ = full_search_serial(t, restarts, seed + 7, kicks)
    return kind, s, words(text), o


def words(text):
    return H.word_score(text.translate(H.TOL))


def show(text):
    for a, b in H.LINE_SPANS:
        log("     " + text[a:b])
    log("     " + H.segment(text.translate(H.TOL)))


def fence_report(ordered_toks, chosen):
    """Which reading was chosen at each fence position (reported by original position)."""
    # map back: ordered_toks[i] is a token; find its original position by identity of order
    return [(m, a, c) for (m, a), c in zip(ordered_toks, chosen) if a is not None]


def pval(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


# ---------------------------------------------------------------- planted test
def planted(seed, lang='en'):
    rng = random.Random(seed)
    plain = H.fold_base(''.join(c for c in SELFTEST_TEXT[lang].upper() if c.isalpha()))[:87]
    letters = sorted(set(plain))
    code_of = dict(zip(letters, rng.sample(H.CODES_ALL, len(letters))))
    true = [code_of[c] for c in plain]
    toks, truth = [], {}
    for i, c in enumerate(true):
        pos = i + 1
        if pos in FENCE:
            loop, d = divmod(c, 10)
            other = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
            if rng.random() < 0.5:
                toks.append((other, c))      # transcriber's main reading is the WRONG one
            else:
                toks.append((c, other))
            truth[pos] = c
        else:
            toks.append((c, None))
    return plain, toks, truth


# ---------------------------------------------------------------- main
def main():
    global LOG
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default=None, help='default english.txt (en) or german.txt (de)')
    ap.add_argument('--lang', choices=['en', 'de'], default='en',
                    help='de = German: umlauts spelled out, Gutenberg licence removed, German self-test')
    ap.add_argument('--out', default=None, help='default v2_out (en) or v2_de_out (de)')
    ap.add_argument('--restarts', type=int, default=40, help='climb restarts per reading order')
    ap.add_argument('--kicks', type=int, default=10, help='shakes per restart (iterated local search)')
    ap.add_argument('--nulls', type=int, default=200, help='scrambled copies per kind')
    ap.add_argument('--kind', choices=['shuffle', 'doublet', 'both'], default='both')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--seed', type=int, default=1897)
    args = ap.parse_args()

    args.corpus = args.corpus or ('english.txt' if args.lang == 'en' else 'german.txt')
    args.out = args.out or ('v2_out' if args.lang == 'en' else 'v2_de_out')
    os.makedirs(args.out, exist_ok=True)
    LOG = open(os.path.join(args.out, 'log.txt'), 'a', encoding='utf-8')
    t0 = time.time()
    log(f"\n[{stamp()}] START dora_v2 {' '.join(sys.argv[1:]) or '(defaults)'}")
    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    corpus = prepare_corpus(args.corpus, args.lang, args.out)
    H.build_scorer(corpus)
    build_display_words(corpus, lang=args.lang)
    ref_plain = H.fold_base(''.join(c for c in SELFTEST_TEXT[args.lang] if c.isalpha()))[:87]
    ref_score = H.SCORE_FN([H.LIDX[c] for c in ref_plain])
    log(f"language: {args.lang}   scorer: {H.MODE}   corpus: {args.corpus}")
    log(f"yardstick: a real {'German' if args.lang == 'de' else 'English'} note of this length scores about {ref_score:.2f}")
    pool = mp.Pool(args.workers, initializer=init_worker, initargs=(corpus, args.lang)) if args.workers > 1 else None

    def run(fn, jobs):
        if pool is None:
            return [fn(*j) for j in jobs]
        return pool.starmap_async(fn, jobs, chunksize=1).get(10 ** 6)

    try:
        if args.selftest:
            plain, toks, truth = planted(args.seed, args.lang)
            log(f"SELF-TEST planted: {plain}")
            wrong = sum(1 for i, (m, a) in enumerate(toks) if a is not None and m != truth[i + 1])
            log(f"10 fence symbols; the main reading is wrong at {wrong} of them")
        else:
            toks = base_tokens()
            log(f"v2 transcription: 87 symbols, fence positions {sorted(FENCE)}")

        # ---- real search: one job per reading order
        jobs = [(toks, o, args.restarts, args.seed + 101 * k, args.kicks) for k, o in enumerate(T.ORDERS)]
        res = run(order_job, jobs)
        res.sort(key=lambda r: -r[1])
        log(f"[{stamp()}] search done in {H.fmt_secs(time.time() - t0)}\n")
        log("Best per reading order (letters AND fence readings solved together):")
        for o, s, text, chosen, ot in res:
            log(f"\n  [{o}]  score {s:.3f}  words {words(text)}")
            show(text)
        o, real, text, chosen, ot = res[0]
        # fence choices, reported in original positions
        orig_pos = T.apply_order(list(range(1, 88)), o)
        picks = sorted((p, m, a, c) for p, (m, a), c in zip(orig_pos, ot, chosen) if a is not None)
        log(f"\nFence readings chosen by the best solution [{o}]:")
        log("  " + '  '.join(f"{p}:{c}{'(alt)' if c == a and a != m else ''}" for p, m, a, c in picks))

        if args.selftest:
            ref = ''.join(T.apply_order(list(plain), o))
            acc = sum(x == y for x, y in zip(text, ref)) / 87
            right = sum(1 for p, m, a, c in picks if c == truth[p])
            log(f"\nSELF-TEST: letters recovered {acc:.0%}, fence readings right {right}/10, "
                f"order picked: {o}")

        # ---- nulls
        kinds = ['shuffle', 'doublet'] if args.kind == 'both' else [args.kind]
        real_w = words(text)
        rows = []
        for kind in kinds:
            log(f"\n[{stamp()}] identical search on {args.nulls} {kind} copies...")
            t1 = time.time()
            done, ns, nw = 0, [], []
            batch = max(1, args.workers) * 4
            while done < args.nulls:
                k = min(batch, args.nulls - done)
                out = run(null_job, [(kind, toks, args.restarts, args.seed + 7919 * (done + i + 1) +
                                      (0 if kind == 'shuffle' else 500000), args.kicks) for i in range(k)])
                for _, s, w, oo in out:
                    ns.append(s)
                    nw.append(w)
                    rows.append((kind, s, w, oo))
                done += k
                el = time.time() - t1
                log(f"  [{stamp()}] {done}/{args.nulls}: best copy so far {max(ns):.3f}, "
                    f"p = {pval(real, ns):.4f}  (~{H.fmt_secs(el / done * (args.nulls - done))} left)")
            log(f"\nRESULT vs {kind}: real best {real:.3f}; copies mean {statistics.mean(ns):.3f}, "
                f"max {max(ns):.3f}; p = {pval(real, ns):.4f} (floor {1 / (1 + len(ns)):.4f})")
            log(f"   words at best: real {real_w}; copies mean {statistics.mean(nw):.1f}, "
                f"max {max(nw)}; p = {pval(real_w, nw):.4f}")
        with open(os.path.join(args.out, 'nulls.csv'), 'w', newline='', encoding='utf-8') as fh:
            w = csv.writer(fh)
            w.writerow(['kind', 'best_score', 'words_at_best', 'order'])
            w.writerows(rows)

        log(f"\nHow to read: real text of this length scores about {ref_score:.1f} (see yardstick above). A hit")
        log("reads as English AND beats both kinds of copy on score and words. Beating only the")
        log("'shuffle' copies means structure, not necessarily language (see the marathon).")
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
