#!/usr/bin/env python3
"""
Dorabella quotation test: does any phrase from texts Elgar and Dora shared (Elgar's
own song words up to 1897, and later the Bible) FIT the cipher, and does fixing it make
the REST of the cipher read as English?

Needs Dorabella_cribs.py and everything it needs (Dorabella_styles.py,
Dorabella_allkeys.py, Dorabella_hillclimb.py, Dorabella_marathon.py, Dorabella_tests.py,
Dorabella_transcription.csv, english.txt, numpy), plus a folder of plain-text sources.

  python Dorabella_quotes.py --selftest                 # phrase counts, power, timing
  python Dorabella_quotes.py                            # quotes\\*.txt -> quotes_out\\
  python Dorabella_quotes.py --source bible --out bible_out --min-len 16
                       # the whole Bible: longer phrases, or the search never finishes
  python Dorabella_quotes.py --report-only / --fresh / --hours N   (as before)

Ctrl+C any time; the same command resumes.

PHRASES: every run of consecutive words in a source file that is MIN-MAX letters long
(default 12-30) after removing everything but letters (J->I, V->U, as in the cipher
alphabet). Duplicates are merged; each phrase remembers where it came from.

PRE-REGISTERED TESTS (two tests: line to beat p < 0.025)
  Q1  Best English score (quadgrams) of the cipher completed around ANY eligible phrase
      in ANY of the 5 reading orders, uncertain symbols free (the crib test's C1
      search). Each scrambled copy gets the identical search with the identical list.
  Q2  How many distinct phrases of 14+ letters fit the cipher's repeat pattern at all.
      A message quoting these texts should let more long phrases fit than chance.
  ELIGIBILITY: decided on 100 pair-preserving copies only (never on the real cipher):
  a phrase is used if it fits on average no more than 20 places per copy.
"""
import argparse
import glob
import json
import math
import multiprocessing as mp
import os
import random
import re
import statistics
import time
from collections import defaultdict

import Dorabella_cribs as C
import Dorabella_styles as S
import Dorabella_hillclimb as H

VERSION = 1
ORDERS = S.ORDERS
ORDER_PERM = S.ORDER_PERM
LONG = 14


# ---------------------------------------------------------------- phrases
def load_phrases(folder, lo, hi):
    files = sorted(glob.glob(os.path.join(folder, '*.txt')))
    if not files:
        raise SystemExit(f"no .txt files found in {folder}")
    src = defaultdict(set)
    per_file = {}
    for fn in files:
        with open(fn, encoding='utf-8', errors='ignore') as fh:
            raw = fh.read()
        s = raw.find('*** START OF')                       # drop Project Gutenberg boilerplate
        if s >= 0:
            raw = raw[raw.find('\n', s) + 1:]
        e = raw.find('*** END OF')
        if e >= 0:
            raw = raw[:e]
        raw = raw.upper().replace("'", '').replace('\u2019', '')
        words = [S.fold(w) for w in re.findall(r'[A-Z]+', raw)]
        name = os.path.splitext(os.path.basename(fn))[0]
        n0 = len(src)
        for i in range(len(words)):
            s = ''
            for j in range(i, len(words)):
                s += words[j]
                if len(s) > hi:
                    break
                if len(s) >= lo and j > i:          # at least two words
                    src[s].add(name)
        per_file[name] = (len(words), len(src) - n0)
    return {p: sorted(v) for p, v in src.items()}, per_file


def canon(seq):
    m, out = {}, []
    for x in seq:
        if x not in m:
            m[x] = len(m)
        out.append(m[x])
    return tuple(out)


def window_patterns(seq, alts, lengths):
    """{length: {pattern: count}} over all windows of all 5 reading orders, every
    combination of uncertain readings inside the window."""
    pats = {L: defaultdict(int) for L in lengths}
    for o in ORDERS:
        R = ORDER_PERM[o]
        for L in lengths:
            for i in range(87 - L + 1):
                win = R[i:i + L]
                opts = [[seq[r]] + ([alts[r]] if r in alts else []) for r in win]
                amb = [k for k, op in enumerate(opts) if len(op) > 1]
                base = [op[0] for op in opts]
                seen = set()
                for bits in range(1 << len(amb)):
                    s = list(base)
                    for b, k in enumerate(amb):
                        if bits >> b & 1:
                            s[k] = opts[k][1]
                    seen.add(canon(s))
                for pt in seen:
                    pats[L][pt] += 1
    return pats


CANON = {}


def fitting(phrases, seq, alts):
    """Phrases whose repeat pattern matches some window (pattern lookup, cached)."""
    lengths = sorted({len(p) for p in phrases})
    pats = window_patterns(seq, alts, lengths)
    out = {}
    for p in phrases:
        c = CANON.get(p)
        if c is None:
            c = CANON[p] = canon(p)
        k = pats[len(p)].get(c)
        if k:
            out[p] = k
    return out


def eligibility(phrases, n=100, seed=4242):
    rng = random.Random(seed)
    tot = defaultdict(int)
    for _ in range(n):
        s, a = S.null_tokens('doublet', S.SEQ, S.ALTS, rng.getrandbits(48))
        for p, k in fitting(phrases, s, a).items():
            tot[p] += k
    return [p for p in phrases if tot[p] / n <= C.ELIG_MAX_FITS]


# ---------------------------------------------------------------- one search
def q_search(pool, corpus, seq, alts, elig, P, seed):
    fits = fitting(elig, seq, alts)
    n_long = sum(1 for p in fits if len(p) >= LONG)
    if not fits:
        return dict(stat=-99.0, crib='', text='', order='', n=0, n_long=0, fits={}, top=[])
    r = C.c1_search(pool, corpus, seq, alts, sorted(fits), P, seed)
    r['n_long'] = n_long
    r['fits'] = {p: k for p, k in fits.items() if len(p) >= LONG}
    return r


# ---------------------------------------------------------------- state
def now():
    import datetime
    return datetime.datetime.now()


def stamp():
    return f"[{now():%Y-%m-%d %H:%M:%S}]"


def blank_state(args, elig, per_file, nphr):
    return dict(version=VERSION, created=str(now()), elapsed=0.0, sessions=0, seq=S.SEQ,
                source=args.source, lo=args.min_len, hi=args.max_len, n_phrases=nphr,
                per_file=per_file, eligible=elig,
                params=dict(seed=args.seed, real_seeds=args.real_seeds, r1=args.r1, i1=args.i1,
                            keep=args.keep, r2=args.r2, i2=args.i2, workers=args.workers),
                targets=dict(shuffle=args.copies, doublet=args.copies, controls=args.controls),
                st=dict(real=[], shuffle=[], doublet=[], controls=[]), log=[])


def next_stage(state):
    st, tg = state['st'], state['targets']
    if len(st['real']) < state['params']['real_seeds']:
        return 'real'
    if len(st['controls']) < min(10, tg['controls']):
        return 'controls'
    fr = {k: len(st[k]) / tg[k] for k in ('doublet', 'shuffle', 'controls')}
    k = min(fr, key=fr.get)
    return k if fr[k] < 1 else None


def do_unit(pool, state, corpus, kind, phrase_src):
    P = state['params']
    st = state['st']
    elig = state['eligible']
    n = len(st[kind])
    seed = S.derive_seed(P['seed'], 'Q', kind, n)
    if kind == 'real':
        r = q_search(pool, corpus, S.SEQ, S.ALTS, elig, P, seed)
        r['fits_src'] = {p: phrase_src.get(p, []) for p in r.get('fits', {})}
        st['real'].append(r)
    elif kind in ('shuffle', 'doublet'):
        s, a = S.null_tokens(kind, S.SEQ, S.ALTS, S.derive_seed(P['seed'], 'Qnull', kind, n))
        r = q_search(pool, corpus, s, a, elig, P, seed)
        st[kind].append(dict(stat=round(r['stat'], 5), n_long=r['n_long'], crib=r['crib']))
    else:
        noise = 0 if n % 2 == 0 else 2
        s, a, plain, crib, order = C.plant_c1(elig, n, noise, S.derive_seed(P['seed'], 'Qplant', n))
        r = q_search(pool, corpus, s, a, elig, P, seed)
        acc = sum(x == y for x, y in zip(r['text'], plain)) / 87 if r['order'] == order else 0.0
        st['controls'].append(dict(stat=round(r['stat'], 5), noise=noise, planted=crib, found=r['crib'],
                                   acc=round(acc, 3), plain=plain, text=r['text']))


# ---------------------------------------------------------------- report
def build_report(state, sd, phrase_src):
    st = state['st']
    pv, pct = C.pval, C.pct
    L = [f"DORABELLA QUOTATION TEST REPORT  (source folder: {state['source']})",
         f"written {now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"phrases of {state['lo']}-{state['hi']} letters: {state['n_phrases']:,}; eligible after the rule: "
         f"{len(state['eligible']):,}", ""]
    L.append("source files (words, new phrases):  " + "; ".join(
        f"{k} ({v[0]}, {v[1]})" for k, v in state['per_file'].items()))
    real = [r['stat'] for r in st['real']]
    med = statistics.median(real) if real else None
    dn = [x['stat'] for x in st['doublet']]
    sn = [x['stat'] for x in st['shuffle']]
    dl = [x['n_long'] for x in st['doublet']]
    rl = st['real'][0]['n_long'] if st['real'] else None
    cs = [c['stat'] for c in st['controls']]
    q95 = pct(dn, .95)
    power = sum(1 for c in cs if q95 is not None and c > q95) / len(cs) if (cs and q95 is not None) else None
    L += ["", "=" * 100, "HEADLINE (p vs PAIR-PRESERVING copies; line to beat p < 0.025)", "=" * 100,
          f"Q1 best completion   real {C.f(med)}   pair n={len(dn)} p={C.f(pv(med, dn), 4)}   "
          f"shuffle n={len(sn)} p={C.f(pv(med, sn), 4)}   power {'-' if power is None else f'{power:.0%}'}",
          f"Q2 long phrases fit  real {rl}   pair mean {C.f(statistics.mean(dl) if dl else None, 2)}   "
          f"p={C.f(pv(rl, dl), 4)}"]
    if st['real']:
        b = max(st['real'], key=lambda r: r['stat'])
        L += ["", "-" * 100, "REAL CIPHER", "-" * 100,
              "Q1 per seed: " + ', '.join(f"{r['stat']:.3f} ({r['crib']})" for r in st['real']),
              f"placements tried: {b['n']}"]
        for x in b.get('top', []):
            srcs = ', '.join(phrase_src.get(x['crib'], []))
            L.append(f"  [{x['score']:.3f}] phrase {x['crib']} ({srcs}) in {x['order']}")
            for a in range(0, 87, 29):
                L.append("      " + x['text'][a:a + 29])
            L.append("      words: " + sd.segment(x['text']))
        L.append(f"\nlong phrases ({LONG}+ letters) that fit the real cipher's pattern: {rl}")
        for p, k in sorted(st['real'][0].get('fits', {}).items(), key=lambda kv: -len(kv[0]))[:40]:
            L.append(f"  {p:<32} {k} place(s)   from {', '.join(phrase_src.get(p, []))}")
    for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
        if nul:
            L.append(f"Q1 vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                     f"{pct(nul, .95):.3f}, max {max(nul):.3f}")
    if st['controls']:
        c = st['controls']
        L.append(f"controls: n={len(c)}, median {statistics.median(cs):.3f}, planted phrase found "
                 f"{sum(1 for x in c if x['found'] == x['planted'])}/{len(c)}, median letters right "
                 f"{statistics.median(x['acc'] for x in c):.0%}")
        ex = max(c, key=lambda x: x['acc'])
        L.append(f"  example: planted {ex['plain'][:56]}\n           found   {ex['text'][:56]}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--source', default='quotes', help='folder of .txt source texts')
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='quotes_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--min-len', type=int, default=12)
    ap.add_argument('--max-len', type=int, default=30)
    ap.add_argument('--real-seeds', type=int, default=5)
    ap.add_argument('--r1', type=int, default=1)
    ap.add_argument('--i1', type=int, default=5000)
    ap.add_argument('--keep', type=int, default=12)
    ap.add_argument('--r2', type=int, default=10)
    ap.add_argument('--i2', type=int, default=20000)
    ap.add_argument('--copies', type=int, default=1000)
    ap.add_argument('--controls', type=int, default=60)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true',
                    help='exit when the copy targets are met (for running several sources back to back)')
    args = ap.parse_args()

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    S.SEQ, S.ALTS = S.load_cipher()
    phrase_src, per_file = load_phrases(args.source, args.min_len, args.max_len)
    phrases = sorted(phrase_src)
    print(f"cipher read; {len(phrases):,} phrases of {args.min_len}-{args.max_len} letters "
          f"from {len(per_file)} files in {args.source}")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))
    sd = S.StyleDict(S.base_dictionary(args.corpus), 'S0')

    def log(msg):
        line = f"{stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != S.SEQ or state['n_phrases'] != len(phrases):
            raise SystemExit("saved state was made with a different transcription or different source "
                             "texts; use --fresh (or another --out)")
    else:
        t0 = time.time()
        print("deciding which phrases are eligible (on 100 scrambled copies only)...", flush=True)
        elig = eligibility(phrases)
        print(f"  {len(elig):,} eligible ({H.fmt_secs(time.time() - t0)})")
        state = blank_state(args, elig, per_file, len(phrases))
    if args.report_only:
        S.write_text(rp, build_report(state, sd, phrase_src))
        return
    state['params']['workers'] = args.workers
    S.keep_awake()
    pool = mp.Pool(args.workers, initializer=C.init_worker, initargs=(args.corpus,)) if args.workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            P = state['params']
            print(f"\n{stamp()} SELF-TEST")
            for k, v in per_file.items():
                print(f"  {k:<32} {v[0]:>5} words, {v[1]:>6} new phrases")
            t0 = time.time()
            fr = fitting(state['eligible'], S.SEQ, S.ALTS)
            print(f"\npattern lookup on one cipher: {H.fmt_secs(time.time() - t0)}")
            print("\ncontrol (planted English containing a random eligible phrase, 2 misreads):")
            s, a, plain, crib, order = C.plant_c1(state['eligible'], 3, 2, 99)
            t0 = time.time()
            r = q_search(pool, args.corpus, s, a, state['eligible'], P, 7)
            dt = time.time() - t0
            acc = sum(x == y for x, y in zip(r['text'], plain)) / 87 if r['order'] == order else 0.0
            print(f"  planted {crib}; found {r['crib']}; score {r['stat']:.3f}; letters right {acc:.0%}; "
                  f"{r['n']} placements; {H.fmt_secs(dt)}")
            print(f"  planted {plain[:60]}\n  found   {r['text'][:60]}")
            ns, na = S.null_tokens('doublet', S.SEQ, S.ALTS, 12345)
            t0 = time.time()
            rn = q_search(pool, args.corpus, ns, na, state['eligible'], P, 8)
            dn = time.time() - t0
            print(f"  scrambled copy: best {rn['stat']:.3f}; {rn['n']} placements; {rn['n_long']} long "
                  f"phrases fit; {H.fmt_secs(dn)}")
            n = 2 * args.copies + args.controls + args.real_seeds
            print(f"\nEstimated time to reach the targets: about {H.fmt_secs(n * (dt + dn) / 2)}")
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {args.workers}, {len(state['eligible']):,} phrases")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            kind = next_stage(state)
            if kind is None:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for k in ('shuffle', 'doublet'):
                    state['targets'][k] = int(state['targets'][k] * 1.5)
                log("all targets met; raising copy targets by 50% and carrying on")
                continue
            t0 = time.time()
            do_unit(pool, state, args.corpus, kind, phrase_src)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['st']
            msg = f"{kind} in {H.fmt_secs(t1 - t0)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{stamp()} {msg}", flush=True)
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, sd, phrase_src))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished unit is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, sd, phrase_src))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    mp.freeze_support()
    main()
