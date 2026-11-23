#!/usr/bin/env python3
"""
Dorabella DOT test (roadmap item D1): what does the dot in line 3 mean?

The cipher has one dot among its symbols: in line 3, between the 5th and 6th symbols (written
positions 65 and 66; symbol 65 is code 38). The two dots in the date below the cipher ("14" and "97")
are ordinary punctuation and are not part of this test.

Needs Dorabella_notebook.py, Dorabella_three.py and everything they need (Dorabella_hillclimb.py,
Dorabella_tests.py, Dorabella_marathon.py, Dorabella_transcription.csv, english.txt).

  python Dorabella_dot.py --selftest
  python Dorabella_dot.py --stop-when-done      # results in dot_out\\
  python Dorabella_dot.py --report-only
  python Dorabella_dot.py --only D1A

Ctrl+C any time; the same command resumes.

PRE-REGISTERED (frozen 21 Nov 2026, before any run; three families, so the line is p < 0.0167)
  D1A  THE DOT IS "&". Elgar wrote "&" for "and" about 9 times in 10 in his private letters, and his
       cipher alphabet has no symbol for it. The plaintext is the decrypt with AND inserted at the
       dot (in the reading direction), scored as ordinary English (Dracula quadgrams).
  D1B  THE DOT IS A BREAK (a full stop, or the end of a word): the text is scored as ordinary
       English, but no four-letter group may span the dot.
  D1C  THE DOT ENDS THE MESSAGE: only the 65 symbols before it are text (the 22 after are filler);
       the 65 are read in the 5 reading orders and scored as ordinary English.
  For each: ANY one-symbol-one-letter key by annealing, over the 5 reading orders (forwards, fully
  reversed, each line reversed, the two boustrophedons; the dot stays between written positions 65
  and 66 in all of them). Statistic: best average quadgram score; real = median of 5 seeds.
  Compared with 1,000 pair-preserving copies (main) and 500 plain shuffles of the whole cipher, with
  the dot kept at the same place and the identical search. Power from planted English written that
  way (D1A: a passage with AND at the dot, AND left out of the symbols; D1B: any passage; D1C: a
  65-letter passage then 22 random symbols), with 0 or 2 misread symbols. A pass also needs a decrypt
  that reads as English. Transcription: the consensus reading (Hauer et al. 2025), as in
  Dorabella_notebook.py; fixed, not searched.
"""
import argparse
import json
import math
import os
import random
import signal
import statistics
import time

import Dorabella_notebook as NB          # consensus reading

Z = NB.Z
H, M = Z.H, Z.M
VERSION = 1
FAMILIES = ['D1A', 'D1B', 'D1C']
DESC = {'D1A': 'the dot is "&": AND inserted at the dot, any key',
        'D1B': 'the dot is a break: no four-letter group spans it, any key',
        'D1C': 'the dot ends the message: only the 65 symbols before it, any key'}
BONF = 0.05 / 3
DOT = 65                                   # the dot lies between written index 64 and 65 (positions 65 and 66)
DEFAULT_TARGETS = {f: dict(shuffle=500, doublet=1000, controls=40) for f in FAMILIES}
SEQ = None


# ============================================================== the text a family scores
def layout(fam, R):
    """For one reading order R (written indexes in reading order): a list of items, each
    ('s', written index) for a symbol or ('f', letter) for a fixed letter, and the list of
    segments [(start, end)] that four-letter groups must stay inside."""
    items, breaks = [], []
    for j, r in enumerate(R):
        if fam == 'D1C' and r >= DOT:
            continue
        items.append(('s', r))
        if j + 1 < len(R) and {r, R[j + 1]} == {DOT - 1, DOT}:
            if fam == 'D1A':
                items += [('f', H.LIDX['A']), ('f', H.LIDX['N']), ('f', H.LIDX['D'])]
            elif fam == 'D1B':
                breaks.append(len(items))
    segs, a = [], 0
    for b in breaks + [len(items)]:
        segs.append((a, b))
        a = b
    return items, segs


def anneal_items(seq, items, segs, restarts, iters, seed):
    """Any key; items as from layout(). Returns (average score, letters)."""
    Q = Z.QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(items)
    valid = [s for a, b in segs for s in range(a, b - 3)]
    vset = set(valid)
    slot_at = [Z.slot_of(seq[x]) if k == 's' else -1 for k, x in items]
    pos = [[] for _ in range(24)]
    for q, s in enumerate(slot_at):
        if s >= 0:
            pos[s].append(q)
    starts = [tuple(sorted({t for q in pos[s] for t in range(q - 3, q + 1) if t in vset})) for s in range(24)]
    used = [s for s in range(24) if pos[s]]
    fac = (Z.SA_T1 / Z.SA_T0) ** (1.0 / max(1, iters))

    def qs(p, ss):
        t = 0.0
        for s in ss:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    best_all = (-1e18, None)
    for _ in range(restarts):
        key = list(range(24))
        rng.shuffle(key)
        p = [key[s] if s >= 0 else x for s, (_, x) in zip(slot_at, items)]
        cur = qs(p, valid)
        best, bp = cur, p[:]
        Tm = Z.SA_T0
        for _it in range(iters):
            i = used[rr(len(used))]
            j = rr(23)
            if j >= i:
                j += 1
            ss = tuple(sorted(set(starts[i]) | set(starts[j])))
            o = qs(p, ss)
            li, lj = key[i], key[j]
            for q in pos[i]:
                p[q] = lj
            for q in pos[j]:
                p[q] = li
            d = qs(p, ss) - o
            if d >= 0 or rand() < exp(d / Tm):
                key[i], key[j] = lj, li
                cur += d
                if cur > best:
                    best, bp = cur, p[:]
            else:
                for q in pos[i]:
                    p[q] = li
                for q in pos[j]:
                    p[q] = lj
            Tm *= fac
        if best > best_all[0]:
            best_all = (best, bp)
    return best_all[0] / len(valid), Z.text_of(best_all[1])


def job(fam, seq, restarts, iters, seed):
    best = (-1e18, '', '')
    for k, o in enumerate(Z.KNOWN_ORDERS):
        items, segs = layout(fam, Z.known_order(o))
        s, text = anneal_items(seq, items, segs, restarts, iters, Z.derive_seed(seed, k))
        if s > best[0]:
            best = (s, text, o)
    return best


def show(fam, text, order):
    """The decrypt with the dot marked."""
    items, _ = layout(fam, Z.known_order(order))
    out = []
    for (k, x), ch in zip(items, text):
        out.append(ch.lower() if k == 'f' else ch)
    s = ''.join(out)
    if fam == 'D1B':
        R = Z.known_order(order)
        j = next(i for i in range(86) if {R[i], R[i + 1]} == {DOT - 1, DOT})
        s = s[:j + 1] + ' . ' + s[j + 1:]
    return s


# ============================================================== controls
def english(idx):
    texts = M.CONTROL_TEXTS
    out, j = '', idx
    while len(out) < 400:
        out += H.fold_base(''.join(ch for ch in texts[j % len(texts)].upper() if ch.isalpha()))
        j += 7
    return out


def plant(fam, idx, noise, seed):
    rng = random.Random(seed)
    o = Z.KNOWN_ORDERS[rng.randrange(5)]
    R = Z.known_order(o)
    items, _ = layout(fam, R)
    nsym = sum(1 for k, _ in items if k == 's')
    big = english(idx)
    if fam == 'D1A':
        k = next(i for i, (kind, _) in enumerate(items) if kind == 'f')          # where AND sits
        hits = [i for i in range(k, len(big) - (len(items) - k)) if big[i:i + 3] == 'AND']
        at = rng.choice(hits) if hits else k
        window = big[at - k:at - k + len(items)]
        letters = [ch for (kind, _), ch in zip(items, window) if kind == 's']
        plain = window
    else:
        start = rng.randrange(0, len(big) - nsym)
        letters = list(big[start:start + nsym])
        plain = ''.join(letters)
    key = list(range(24))
    rng.shuffle(key)
    kinv = Z.inverse(key)
    seq = [Z.code_of_slot(rng.randrange(24)) for _ in range(87)]       # D1C: the 22 after the dot are filler
    sym = [x for kind, x in items if kind == 's']
    for r, ch in zip(sym, letters):
        seq[r] = Z.code_of_slot(kinv[H.LIDX[ch]])
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, o


# ============================================================== batches
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    Z.init_tables(corpus)


def null_seq(kind, fam, n, seed0):
    rng = random.Random(Z.derive_seed(seed0, 'null', fam, kind, n))
    if kind == 'shuffle':
        s = list(SEQ)
        rng.shuffle(s)
        return s
    return M.doublet_shuffle(list(SEQ), rng)


def do_batch(pool, state, fam, kind, count):
    st, P = state['fam'][fam], state['params']
    seed0 = P['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append(SEQ)
            meta.append(None)
        elif kind == 'controls':
            noise = 0 if n % 2 == 0 else 2
            s, plain, o = plant(fam, n, noise, Z.derive_seed(seed0, 'plant', fam, n))
            items.append(s)
            meta.append((noise, plain, o))
        else:
            items.append(null_seq(kind, fam, n, seed0))
            meta.append(None)
    res = Z.run_jobs(pool, job, [(fam, s, P['restarts'], P['iters'], Z.derive_seed(seed0, fam, kind, len(st[kind]) + k))
                                 for k, s in enumerate(items)])
    for (s, text, o), mt in zip(res, meta):
        rec = dict(stat=round(s, 5), order=o)
        if kind == 'real':
            rec['text'] = text
        elif kind == 'controls':
            noise, plain, po = mt
            rec.update(noise=noise, plain=plain, text=text, right_order=int(o == po),
                       acc=round(sum(a == b for a, b in zip(text, plain)) / len(plain), 3) if o == po else 0.0)
        st[kind].append(rec)


def stages(state, only):
    out = []
    for f in (only or FAMILIES):
        st, tg = state['fam'][f], state['targets'][f]
        if len(st['real']) < state['params']['real_seeds']:
            out.append((0, len(st['real']) / state['params']['real_seeds'], f, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], f, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            fr = len(st[kind]) / tg[kind]
            if fr < 1:
                out.append((2, fr, f, kind))
    return sorted(out)


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def build_report(state):
    L = ["DORABELLA DOT REPORT  (roadmap D1: the dot in line 3, between written positions 65 and 66)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   sessions {state['sessions']}",
         NB.READ_NOTE,
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing {state['params']['restarts']} "
         f"restarts x {state['params']['iters']} steps per reading order",
         f"three pre-registered families: the line is p < {BONF:.4f}",
         "decrypts: in D1A the inserted 'and' is shown in small letters; in D1B the dot is shown as ' . '", ""]
    for f in FAMILIES:
        st = state['fam'][f]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        q95 = Z.pct(dn, .95) if dn else None
        cs = st['controls']
        power = sum(1 for c in cs if q95 is not None and c['stat'] > q95) / len(cs) if (cs and q95 is not None) else None
        p_d = Z.pval(med, dn) if (dn and med is not None) else None
        L += ["=" * 104, f"{f}  {DESC[f]}", "=" * 104]
        if med is None:
            L += ["real run not done yet", ""]
            continue
        verdict = ('waiting for copies' if p_d is None else
                   'BEATS the copies past the line: read the decrypt by eye' if p_d < BONF else
                   'chance level, LOW POWER: inconclusive' if (power is not None and power < 0.5) else 'chance level')
        L.append(f"real {med:.3f} (seeds: {', '.join(f'{x:.3f}' for x in real)})   p vs pair-preserving {fmt(p_d, 4)} "
                 f"(n={len(dn)})   p vs plain shuffle {fmt(Z.pval(med, sn) if sn else None, 4)} (n={len(sn)})   "
                 f"power {'-' if power is None else f'{power:.0%}'}")
        L.append(f"verdict: {verdict}")
        for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L.append(f"  vs {name:<16} mean {statistics.mean(nul):.3f}, 95th pct {Z.pct(nul, .95):.3f}, max {max(nul):.3f}")
        L.append("best real decrypts (best first, repeats left out):")
        seen = set()
        for r in sorted(st['real'], key=lambda r: -r['stat']):
            if r['text'] in seen:
                continue
            seen.add(r['text'])
            t = show(f, r['text'], r['order'])
            L.append(f"  [{r['stat']:.3f}] {r['order']}")
            L.append("      " + t)
            L.append("      words: " + Z.show_words(r['text']))
            if len(seen) == 3:
                break
        for nz in (0, 2):
            c2 = [c for c in cs if c['noise'] == nz]
            if c2:
                L.append(f"controls, {nz} misread symbols: n={len(c2)}, right reading order {sum(c['right_order'] for c in c2)}/{len(c2)}, "
                         f"median letters right {statistics.median(c['acc'] for c in c2):.0%}, above null 95th pct "
                         f"{sum(1 for c in c2 if q95 is not None and c['stat'] > q95)}/{len(c2)}")
        if cs:
            ex = max(cs, key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:60]}")
            L.append(f"                                   found   {ex['text'][:60]}")
        L.append("")
    L += ["A family matters only if the real cipher beats the PAIR-PRESERVING copies past the line, its controls",
          "show it can find a message, and the decrypt reads as English.", "",
          "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    print(f"symbols either side of the dot: {SEQ[DOT - 1]} . {SEQ[DOT]}  (written positions {DOT} and {DOT + 1})")
    for f in FAMILIES:
        items, segs = layout(f, Z.known_order('forward'))
        print(f"{f}: {len(items)} letters scored, segments {segs}")
    for f in FAMILIES:
        seq, plain, o = plant(f, 3, 0, Z.derive_seed('selftest', f))
        t0 = time.time()
        s, text, ro = Z.run_jobs(pool, job, [(f, seq, args.restarts, args.iters, 7)])[0]
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(text, plain)) / len(plain) if ro == o else 0.0
        ns = Z.run_jobs(pool, job, [(f, null_seq('doublet', f, 0, 12345), args.restarts, args.iters, 8)])[0][0]
        print(f"\n{f}: planted control {s:.3f} (a pair-preserving copy: {ns:.3f}); reading order "
              f"{'right' if ro == o else 'wrong'}; letters right {acc:.0%}; {H.fmt_secs(dt)} per search")
        print(f"  planted {plain[:62]}\n  decoded {text[:62]}")
        n = sum(DEFAULT_TARGETS[f].values()) + 5
        print(f"  estimated time for {n} searches, {workers} at a time: about {H.fmt_secs(n * dt / workers)}")


# ============================================================== main
def main():
    global SEQ
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='dot_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--restarts', type=int, default=8)
    ap.add_argument('--iters', type=int, default=15000)
    ap.add_argument('--only', nargs='+', choices=FAMILIES)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true')
    args = ap.parse_args()
    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    SEQ, source = Z.load_seq()
    print(f"cipher read from: {source}")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))

    def log(msg):
        line = f"{Z.stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    Z.init_tables(args.corpus)
    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != SEQ:
            raise SystemExit("saved state was made with a different transcription; use --fresh or another --out")
    else:
        state = dict(version=VERSION, created=str(Z.now()), elapsed=0.0, sessions=0, seq=SEQ,
                     params=dict(seed=args.seed, restarts=args.restarts, iters=args.iters, real_seeds=5),
                     targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                     fam={f: dict(real=[], doublet=[], shuffle=[], controls=[]) for f in FAMILIES}, log=[])
    if args.report_only:
        Z.write_text(rp, build_report(state))
        print(f"report written to {rp}")
        return
    workers = max(1, args.workers)
    Z.keep_awake()
    pool = Z.mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stg = stages(state, args.only)
            if not stg:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for f in (args.only or FAMILIES):
                    for kind in ('doublet', 'shuffle'):
                        state['targets'][f][kind] = int(state['targets'][f][kind] * 1.5)
                log("all targets met; raising the copy targets by 50% and carrying on")
                continue
            _, _, f, kind = stg[0]
            n = 1 if kind == 'real' else max(1, min(workers, state['targets'][f][kind] - len(state['fam'][f][kind])))
            t0 = time.time()
            do_batch(pool, state, f, kind, n)
            now_t = time.time()
            state['elapsed'] += now_t - last
            last = now_t
            st = state['fam'][f]
            msg = f"{f} {kind} x{n} in {H.fmt_secs(now_t - t0)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][f][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{Z.stamp()} {msg}", flush=True)
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    Z.mp.freeze_support()
    main()
