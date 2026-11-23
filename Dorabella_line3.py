#!/usr/bin/env python3
"""
Dorabella LINE 3 WORD-GAP test (roadmap item G1): is line 3 really written with word gaps?

WHY
In Dorabella_spaces.py the word-gap model (one symbol = the gap between words) beat its copies,
p = 0.005, but the decrypt did not read. Looking afterwards: the search used symbol 38 (three loops,
opening up-left) as the gap. It ends lines 1 and 2 and recurs through line 3, cutting it into
"words" of 1, 2, 5, 5, 4 and 5 letters, but never appears inside lines 1 and 2. So the pass may
belong to line 3 alone. This script tests that idea on each part of the cipher separately.

Needs Dorabella_spaces.py, Dorabella_cribs.py and everything they need (Dorabella_languages.py,
Dorabella_styles.py, Dorabella_hillclimb.py, Dorabella_marathon.py, Dorabella_tests.py,
Dorabella_transcription.csv, english.txt).

  python Dorabella_line3.py --selftest          # checks, one control per test, time estimate
  python Dorabella_line3.py --stop-when-done    # results in line3_out\\
  python Dorabella_line3.py --report-only
  python Dorabella_line3.py --only G1B          # work on one test

Ctrl+C any time; the same command resumes.

PRE-REGISTERED (frozen 19 Nov 2026, before any run; three tests, so the line is p < 0.0167)
  The model is the word-gap English of Dorabella_spaces.py (Dracula with each gap between words
  written as one extra letter, Z; real Z's written S), quadgram model from the first 85% of the
  book, planted text from the last 15%. Each part of the cipher is treated as a cipher of its own:
  its own search, and its own scrambled copies made from that part alone.
  G1A  LINE 3 ALONE (27 symbols), ANY key: the search finds which symbol, if any, is the gap.
       Reading orders: forwards and backwards.
  G1B  LINE 3 ALONE, SYMBOL 38 IS THE GAP: the same, with symbol 38 fixed as the gap and the rest of
       the key searched. Line 3 holds five 38s (at line positions 2, 5, 11, 17 and 22); a
       pair-preserving copy of line 3 has the same five, and they are fixed as the gap in the same way.
       HONESTY NOTE: this reading was found by looking at the Dorabella_spaces.py result, so G1B
       cannot be independent evidence for it. A pass means only that line 3, read this way, is more
       English than scrambled versions of itself; it still has to read.
  G1C  LINES 1 AND 2 (60 symbols), ANY key: is there any gap symbol there? Reading orders: forwards,
       backwards, each line backwards, and the two boustrophedons.
  For each: statistic = best quadgram score over the reading orders; real = median of 9 seeds;
  1,000 pair-preserving copies (main) and 500 plain shuffles of that part, through the identical
  search; power from planted word-gap English of the same length (0 or 2 misread symbols; for G1B
  the gap written with symbol 38). The uncertain symbols are those of Dorabella_spaces.py (the CSV's
  9, plus position 32 read as 18 or 11), chosen inside the search. A pass also needs a decrypt that
  reads as English words. Short texts give the search a lot of freedom, so the report states power.
"""
import argparse
import json
import os
import random
import signal
import statistics
import time

import Dorabella_spaces as SP              # the word-gap language (WS) and position 32 reopened
import Dorabella_cribs as C               # anneal_fixed: any key with some symbols fixed

LG = SP.LG
S, H, M = LG.S, LG.H, LG.M
VERSION = 1
TESTS = ['G1A', 'G1B', 'G1C']
DESC = {'G1A': 'line 3 alone (27 symbols), any key; the search picks the gap symbol',
        'G1B': 'line 3 alone, symbol 38 fixed as the word gap, the rest of the key searched',
        'G1C': 'lines 1 and 2 (60 symbols), any key; the search picks the gap symbol'}
BONF = 0.05 / 3
GAP = SP.GAP
GAP_CODE = 38
DEFAULT_TARGETS = {t: dict(shuffle=500, doublet=1000, controls=40) for t in TESTS}
SEQ = ALTS = None
HELD = None


# ============================================================== parts and reading orders
def part(test):
    return (60, 87) if test in ('G1A', 'G1B') else (0, 60)


def orders(test):
    """(name, permutation of the part's positions)."""
    if test in ('G1A', 'G1B'):
        r = list(range(27))
        return [('forwards', r), ('backwards', r[::-1])]
    l1, l2 = list(range(29)), list(range(29, 60))
    return [('forwards', l1 + l2), ('backwards', (l1 + l2)[::-1]), ('each line backwards', l1[::-1] + l2[::-1]),
            ('boustrophedon A', l1 + l2[::-1]), ('boustrophedon B', l1[::-1] + l2)]


def part_of(test, seq, alts):
    a, b = part(test)
    return list(seq[a:b]), {q - a: v for q, v in alts.items() if a <= q < b}


# ============================================================== search
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    LG.use_language('WS')


def job(test, seq, alts, restarts, iters, seed):
    LG.use_language('WS')
    fixed = {S.slot_of(GAP_CODE): H.LIDX[GAP]} if test == 'G1B' else {}
    best = (-1e18, '', '')
    for k, (name, R) in enumerate(orders(test)):
        c = [S.slot_of(seq[r]) for r in R]
        A = [S.slot_of(alts[r]) if r in alts else -1 for r in R]
        s, text = C.anneal_fixed(c, A, fixed, set(), restarts, iters, S.derive_seed(seed, k))
        if s > best[0]:
            best = (s, text, name)
    return best


# ============================================================== controls
def plant(test, idx, noise, seed):
    rng = random.Random(seed)
    a, b = part(test)
    n = b - a
    start = rng.randrange(0, len(HELD) - n)
    plain = HELD[start:start + n]
    slots = list(range(24))
    rng.shuffle(slots)
    if test == 'G1B':                                   # the gap is written with symbol 38
        g = S.slot_of(GAP_CODE)
        z = H.LIDX[GAP]
        k = slots.index(g)
        slots[k], slots[z] = slots[z], g
    name, R = orders(test)[rng.randrange(len(orders(test)))]
    seq = [0] * n
    for i, ch in enumerate(plain):
        seq[R[i]] = S.code_of_slot(slots[H.LIDX[ch]])
    _, palts = part_of(test, SEQ, ALTS)
    for q in rng.sample([i for i in range(n) if i not in palts], noise):
        seq[q] = S.neighbour(seq[q], rng)
    alts = {}
    for q in sorted(palts):
        wrong = S.neighbour(seq[q], rng)
        if test == 'G1B' and (seq[q] == GAP_CODE or wrong == GAP_CODE):
            continue
        if rng.random() < 0.5:
            alts[q] = wrong
        else:
            alts[q], seq[q] = seq[q], wrong
    return seq, alts, plain, name


# ============================================================== batches
def do_batch(pool, state, test, kind, count):
    st, P = state['st'][test], state['params']
    seed0 = P['seed']
    rseq, ralts = part_of(test, SEQ, ALTS)
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append((rseq, ralts))
            meta.append(None)
        elif kind == 'controls':
            noise = 0 if n % 2 == 0 else 2
            s, a, plain, name = plant(test, n, noise, S.derive_seed(seed0, 'plant', test, n))
            items.append((s, a))
            meta.append((noise, plain, name))
        else:
            items.append(S.null_tokens(kind, rseq, ralts, S.derive_seed(seed0, 'null', test, kind, n)))
            meta.append(None)
    args = [(test, s, a, P['restarts'], P['iters'], S.derive_seed(seed0, test, kind, len(st[kind]) + k))
            for k, (s, a) in enumerate(items)]
    res = _run(pool, args)
    for (s, text, name), mt in zip(res, meta):
        rec = dict(stat=round(s, 5), order=name)
        if kind == 'real':
            rec['text'] = text
        elif kind == 'controls':
            noise, plain, pname = mt
            rec.update(noise=noise, plain=plain, text=text, right_order=int(name == pname),
                       acc=round(sum(x == y for x, y in zip(text, plain)) / len(plain), 3) if name == pname else 0.0)
        st[kind].append(rec)


def _run(pool, args):
    if pool is None:
        return [job(*a) for a in args]
    res = pool.starmap_async(job, args, chunksize=1)
    while True:
        try:
            return res.get(60)
        except LG.mp.TimeoutError:
            continue


def stages(state, only):
    out = []
    for t in (only or TESTS):
        st, tg = state['st'][t], state['targets'][t]
        if len(st['real']) < state['params']['real_seeds']:
            out.append((0, len(st['real']) / state['params']['real_seeds'], t, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], t, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            fr = len(st[kind]) / tg[kind]
            if fr < 1:
                out.append((2, fr, t, kind))
    return sorted(out)


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def pval(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))]


def gapped(text):
    return text.replace(GAP, ' ').strip()


def build_report(state):
    L = ["DORABELLA LINE 3 WORD-GAP REPORT  (roadmap G1)",
         f"written {S.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   sessions {state['sessions']}",
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing {state['params']['restarts']} "
         f"restarts x {state['params']['iters']} steps per reading order; uncertain symbols chosen in the search",
         f"model: English with each word gap written as one letter ({GAP}); decrypts are shown with the gaps as spaces",
         f"three pre-registered tests: the line is p < {BONF:.4f}", ""]
    for t in TESTS:
        st = state['st'][t]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        q95 = pct(dn, .95) if dn else None
        cs = st['controls']
        power = sum(1 for c in cs if q95 is not None and c['stat'] > q95) / len(cs) if (cs and q95 is not None) else None
        p_d = pval(med, dn) if (dn and med is not None) else None
        L += ["=" * 104, f"{t}  {DESC[t]}", "=" * 104]
        if med is None:
            L += ["real run not done yet", ""]
            continue
        verdict = ('waiting for copies' if p_d is None else
                   'BEATS the copies past the line: read the decrypt by eye' if p_d < BONF else
                   'chance level, LOW POWER: inconclusive' if (power is not None and power < 0.5) else 'chance level')
        L.append(f"real {med:.3f} (seeds: {', '.join(f'{x:.3f}' for x in real)})   p vs pair-preserving {fmt(p_d, 4)} "
                 f"(n={len(dn)})   p vs plain shuffle {fmt(pval(med, sn) if sn else None, 4)} (n={len(sn)})   "
                 f"power {'-' if power is None else f'{power:.0%}'}")
        L.append(f"verdict: {verdict}")
        for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L.append(f"  vs {name:<16} mean {statistics.mean(nul):.3f}, 95th pct {pct(nul, .95):.3f}, max {max(nul):.3f}")
        L.append("best real decrypts (one per seed, best first):")
        seen = set()
        for r in sorted(st['real'], key=lambda r: -r['stat']):
            if r['text'] in seen:
                continue
            seen.add(r['text'])
            L.append(f"  [{r['stat']:.3f}] {r['order']}:  {gapped(r['text'])}")
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
            L.append(f"  example control ({ex['acc']:.0%} right): planted {gapped(ex['plain'])}")
            L.append(f"                                   found   {gapped(ex['text'])}")
        L.append("")
    L += ["A test matters only if the real part beats its PAIR-PRESERVING copies past the line, its controls",
          "show it can find a message, and the decrypt reads as English words.", "",
          "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, workers):
    print(f"\n{S.stamp()} SELF-TEST")
    s3, _ = part_of('G1A', SEQ, ALTS)
    print(f"line 3: {s3}  (symbol 38 at line positions {[i + 1 for i, c in enumerate(s3) if c == GAP_CODE]})")
    assert s3.count(GAP_CODE) == 5
    P = dict(restarts=args.restarts, iters=args.iters)
    for t in TESTS:
        seq, alts, plain, name = plant(t, 3, 0, S.derive_seed('selftest', t))
        t0 = time.time()
        s, text, o = _run(pool, [(t, seq, alts, P['restarts'], P['iters'], 7)])[0]
        dt = time.time() - t0
        acc = sum(x == y for x, y in zip(text, plain)) / len(plain) if o == name else 0.0
        rseq, ralts = part_of(t, SEQ, ALTS)
        ns, na = S.null_tokens('doublet', rseq, ralts, 12345)
        nul = _run(pool, [(t, ns, na, P['restarts'], P['iters'], 8)])[0]
        print(f"\n{t}: planted control {s:.3f} (a pair-preserving copy: {nul[0]:.3f}); reading order "
              f"{'right' if o == name else 'wrong'}; letters right {acc:.0%}; {H.fmt_secs(dt)} per search")
        print(f"  planted {gapped(plain)}\n  decoded {gapped(text)}")
        n = sum(DEFAULT_TARGETS[t].values()) + 9
        print(f"  estimated time for {n} searches, {workers} at a time: about {H.fmt_secs(n * dt / workers)}")


# ============================================================== main
def main():
    global SEQ, ALTS, HELD
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='line3_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--restarts', type=int, default=10)
    ap.add_argument('--iters', type=int, default=15000)
    ap.add_argument('--only', nargs='+', choices=TESTS)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true')
    args = ap.parse_args()
    if not os.path.exists('english.txt'):
        raise SystemExit("english.txt not found")
    SEQ, ALTS = S.load_cipher()
    LG.SEQ, LG.ALTS = SEQ, ALTS
    S.SEQ, S.ALTS = SEQ, ALTS
    HELD = LG.load_language('WS')[2]
    print(f"cipher read from Dorabella_transcription.csv ({len(ALTS)} uncertain symbols, position 32 included)")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))

    def log(msg):
        line = f"{S.stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != SEQ:
            raise SystemExit("saved state was made with a different transcription; use --fresh or another --out")
    else:
        state = dict(version=VERSION, created=str(S.now()), elapsed=0.0, sessions=0, seq=SEQ,
                     params=dict(seed=args.seed, restarts=args.restarts, iters=args.iters, real_seeds=9),
                     targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                     st={t: dict(real=[], doublet=[], shuffle=[], controls=[]) for t in TESTS}, log=[])
    if args.report_only:
        S.write_text(rp, build_report(state))
        print(f"report written to {rp}")
        return
    workers = max(1, args.workers)
    S.keep_awake()
    pool = LG.mp.Pool(workers, initializer=init_worker) if workers > 1 else None
    if pool is None:
        LG.use_language('WS')
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
                for t in (args.only or TESTS):
                    for kind in ('doublet', 'shuffle'):
                        state['targets'][t][kind] = int(state['targets'][t][kind] * 1.5)
                log("all targets met; raising the copy targets by 50% and carrying on")
                continue
            _, _, t, kind = stg[0]
            n = 1 if kind == 'real' else max(1, min(workers, state['targets'][t][kind] - len(state['st'][t][kind])))
            t0 = time.time()
            do_batch(pool, state, t, kind, n)
            now_t = time.time()
            state['elapsed'] += now_t - last
            last = now_t
            st = state['st'][t]
            msg = f"{t} {kind} x{n} in {H.fmt_secs(now_t - t0)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][t][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{S.stamp()} {msg}", flush=True)
            S.save_json(sp, state)
            S.write_text(rp, build_report(state))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            S.save_json(sp, state)
            S.write_text(rp, build_report(state))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    LG.mp.freeze_support()
    main()
