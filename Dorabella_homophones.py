#!/usr/bin/env python3
"""
Dorabella HOMOPHONES test: could several symbols stand for the SAME letter?

Needs Dorabella_three.py and what it needs (Dorabella_hillclimb.py, Dorabella_tests.py,
Dorabella_marathon.py), Dorabella_transcription.csv, english.txt and numpy in this folder.

  python Dorabella_homophones.py --selftest          # checks + one planted control per family + time estimate
  python Dorabella_homophones.py --stop-when-done    # results in homophones_out\\
  python Dorabella_homophones.py --report-only       # rewrite the report from saved state

Stop any time with Ctrl+C; the same command resumes.

THE IDEA
  On his 1920 notebook page Elgar wrote out "DO YOU GO TO LONDON TOMORROW?", ticked every O,
  and wrote "9 O's" and "23": nine of the 23 letters are one letter. A cipher-maker who counts
  a common letter is looking at the weakness of a one-symbol-one-letter cipher (a common
  letter shows itself), and the classic cure is to give common letters SEVERAL symbols
  (homophones), leaving rare letters (J, Q, X, Z) with none.
  Every key search so far was one-to-one (each symbol a different letter), so it could not
  have found such a key. This test lets several symbols share one letter.

PRE-REGISTERED FAMILIES (frozen 9 Oct 2026, before any run; 3 families, Bonferroni p < 0.0167)
  H2  any key in which a letter may have up to 2 symbols
  H3  any key in which a letter may have up to 3 symbols
  HF  homophones only for the common letters: E T A O up to 3 symbols, I N S H R up to 2,
      every other letter at most 1 (closest to Elgar's O-count)

  Each family: key found by annealing (6 restarts x 60000 steps), in each of the 5 reading
  orders (forward, fully reversed, each line reversed, two boustrophedons); statistic = best
  quadgram score over the 5 orders; real = median of 5 search seeds; symbols read as in the
  transcription. Compared with pair-preserving copies (the main comparison) and plain
  shuffles that get EXACTLY the same freedom. More freedom raises the copies' scores too, so
  planted controls (English written with homophones in that family's pattern, random reading
  order, with 0 or 2 misread symbols) measure how often a real message would stand out.
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import random
import signal
import statistics
import time

import Dorabella_three as Z

H, M = Z.H, Z.M
VERSION = 1
FAMILIES = ['H2', 'H3', 'HF']
FAMILY_DESC = {'H2': 'any key, up to 2 symbols per letter',
               'H3': 'any key, up to 3 symbols per letter',
               'HF': 'homophones for common letters only (ETAO up to 3, INSHR up to 2, others 1)'}
BONF = 0.05 / len(FAMILIES)
ORDERS = Z.KNOWN_ORDERS
DEFAULT_TARGETS = {f: dict(shuffle=500, doublet=1000, controls=40) for f in FAMILIES}
SA_T0, SA_T1 = 10.0, 0.3


def caps_of(fam):
    if fam == 'H2':
        return [2] * 24
    if fam == 'H3':
        return [3] * 24
    c = [1] * 24
    for ch in 'ETAO':
        c[H.LIDX[ch]] = 3
    for ch in 'INSHR':
        c[H.LIDX[ch]] = 2
    return c


CAPS = {f: caps_of(f) for f in FAMILIES}
ORDER_R = None


def order_perms():
    global ORDER_R
    if ORDER_R is None:
        ORDER_R = [Z.known_order(o) for o in ORDERS]
    return ORDER_R


# ============================================================== annealer
def h_anneal(c, caps, restarts, iters, seed):
    """c: symbol slots (0-23) in reading order. Many-to-one key: each letter L may be given to
    at most caps[L] symbols. Returns (average score, letters, key {slot: letter})."""
    Q = Z.QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(c)
    used = sorted(set(c))
    nu = len(used)
    pos = {s: [q for q in range(n) if c[q] == s] for s in used}
    starts = {s: tuple(sorted({t for q in pos[s] for t in range(max(0, q - 3), min(q, n - 4) + 1)}))
              for s in used}
    pair = {}
    fac = (SA_T1 / SA_T0) ** (1.0 / max(1, iters))
    allst = range(n - 3)

    def full(p):
        t = 0.0
        for s in allst:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    def part(p, ss):
        t = 0.0
        for s in ss:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    best_all = (-1e18, None, None)
    for _ in range(restarts):
        letters = list(range(24))
        rng.shuffle(letters)
        key = {s: letters[i] for i, s in enumerate(used)}      # one-to-one start (always allowed)
        cnt = [0] * 24
        for s in used:
            cnt[key[s]] += 1
        p = [key[x] for x in c]
        cur = full(p)
        best, bp, bk = cur, p[:], dict(key)
        Tm = SA_T0
        for _it in range(iters):
            i = used[rr(nu)]
            li = key[i]
            if rand() < 0.5:                       # give symbol i a new letter
                L = rr(24)
                if L == li or cnt[L] >= caps[L]:
                    Tm *= fac
                    continue
                ss = starts[i]
                o = part(p, ss)
                for q in pos[i]:
                    p[q] = L
                d = part(p, ss) - o
                if d >= 0 or rand() < exp(d / Tm):
                    key[i] = L
                    cnt[li] -= 1
                    cnt[L] += 1
                    cur += d
                    if cur > best:
                        best, bp, bk = cur, p[:], dict(key)
                else:
                    for q in pos[i]:
                        p[q] = li
            else:                                  # swap the letters of two symbols
                j = used[rr(nu)]
                lj = key[j]
                if lj == li:
                    Tm *= fac
                    continue
                a, b = (i, j) if i < j else (j, i)
                ss = pair.get((a, b))
                if ss is None:
                    ss = tuple(sorted(set(starts[a]) | set(starts[b])))
                    pair[(a, b)] = ss
                o = part(p, ss)
                for q in pos[i]:
                    p[q] = lj
                for q in pos[j]:
                    p[q] = li
                d = part(p, ss) - o
                if d >= 0 or rand() < exp(d / Tm):
                    key[i], key[j] = lj, li
                    cur += d
                    if cur > best:
                        best, bp, bk = cur, p[:], dict(key)
                else:
                    for q in pos[i]:
                        p[q] = li
                    for q in pos[j]:
                        p[q] = lj
            Tm *= fac
        assert abs(best - full(bp)) < 1e-6, "annealer drifted"
        if best > best_all[0]:
            best_all = (best, bp, bk)
    s, bp, bk = best_all
    return s / (n - 3), bp, bk


def homophone_note(key):
    """Letters that got more than one symbol, e.g. 'E x3, O x2'."""
    from collections import Counter
    c = Counter(key.values())
    multi = sorted(((v, k) for k, v in c.items() if v > 1), reverse=True)
    return ', '.join(f"{H.LETTERS[k]} x{v}" for v, k in multi) or 'none (one-to-one)'


# ============================================================== jobs
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    Z.init_tables(corpus)


def search_job(seq, fam, oi, restarts, iters, seed):
    R = order_perms()[oi]
    slots = [Z.slot_of(x) for x in seq]
    s, p, key = h_anneal([slots[r] for r in R], CAPS[fam], restarts, iters, seed)
    return s, Z.text_of(p), ORDERS[oi], homophone_note(key)


def run_family(pool, fam, items, P):
    """items: list of (seq, seed). Returns per item the list of (score, text, order, note)."""
    jobs, owner = [], []
    for k, (seq, seed) in enumerate(items):
        for oi in range(len(ORDERS)):
            jobs.append((seq, fam, oi, P['restarts'], P['iters'], Z.derive_seed(seed, fam, oi)))
            owner.append(k)
    res = Z.run_jobs(pool, search_job, jobs)
    out = [[] for _ in items]
    for k, r in zip(owner, res):
        out[k].append(r)
    return out


def summarize(res):
    s, text, order, note = max(res, key=lambda r: r[0])
    return dict(stat=round(s, 5), order=order, text=text, note=note)


# ============================================================== controls
def plant(fam, idx, noise, seed):
    """Planted English written with homophones in this family's pattern: every letter in the
    message gets a symbol, then spare symbols go to the commonest letters (up to the family's
    limit) until all 24 symbols are used; each occurrence uses one of its letter's symbols at
    random. Random reading order. Returns (codes, plaintext, order name)."""
    rng = random.Random(seed)
    texts = M.CONTROL_TEXTS
    raw = texts[idx % len(texts)] + ' ' + texts[(idx + 7) % len(texts)]
    plain = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))[:87]
    L = [H.LIDX[ch] for ch in plain]
    caps = CAPS[fam]
    from collections import Counter
    freq = Counter(L)
    want = {x: 1 for x in freq}
    spare = 24 - len(want)
    while spare > 0:
        cands = [x for x in sorted(freq, key=lambda x: -freq[x]) if want[x] < caps[x]]
        if not cands:
            break
        # the letter with the most uses per symbol so far gets the next symbol
        x = max(cands, key=lambda y: freq[y] / want[y])
        want[x] += 1
        spare -= 1
    slots = list(range(24))
    rng.shuffle(slots)
    homs, k = {}, 0
    for x in sorted(want):
        homs[x] = slots[k:k + want[x]]
        k += want[x]
    oi = rng.randrange(len(ORDERS))
    R = order_perms()[oi]
    seq = [0] * 87
    for i, r in enumerate(R):
        seq[r] = Z.code_of_slot(rng.choice(homs[L[i]]))
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    note = ', '.join(f"{H.LETTERS[x]} x{v}" for x, v in sorted(want.items(), key=lambda t: -t[1]) if v > 1)
    return seq, plain, ORDERS[oi], note


# ============================================================== state and batches
def blank_state(args):
    return dict(version=VERSION, created=str(Z.now()), elapsed=0.0, sessions=0,
                params=dict(seed=args.seed, real_seeds=args.real_seeds, restarts=args.restarts,
                            iters=args.iters),
                targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                fam={f: dict(real=[], shuffle=[], doublet=[], controls=[]) for f in FAMILIES},
                unit_secs={}, log=[])


def stages(state, only):
    out = []
    for f in (only or FAMILIES):
        st, tg = state['fam'][f], state['targets'][f]
        if len(st['real']) < state['params']['real_seeds']:
            out.append((0, -1.0, f, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], f, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            frac = len(st[kind]) / tg[kind]
            if frac < 1:
                out.append((2, frac, f, kind))
    return sorted(out)


def batch_size(state, f, kind, workers):
    if kind == 'real':
        return max(1, state['params']['real_seeds'] - len(state['fam'][f]['real']))
    n = max(1, math.ceil(2 * workers / len(ORDERS)))
    return max(1, min(n, state['targets'][f][kind] - len(state['fam'][f][kind])))


def do_batch(pool, state, SEQ, f, kind, count):
    st, P = state['fam'][f], state['params']
    seed0 = P['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append((SEQ, Z.derive_seed(seed0, 'real', f, n)))
        elif kind in ('shuffle', 'doublet'):
            rng = random.Random(Z.derive_seed(seed0, 'null', f, kind, n))
            if kind == 'shuffle':
                s = list(SEQ)
                rng.shuffle(s)
            else:
                s = M.doublet_shuffle(list(SEQ), rng)
            items.append((s, Z.derive_seed(seed0, 'nullsearch', f, kind, n)))
        else:
            noise = 0 if n % 2 == 0 else 2
            seq, plain, order, note = plant(f, n, noise, Z.derive_seed(seed0, 'plant', f, n))
            items.append((seq, Z.derive_seed(seed0, 'ctrlsearch', f, n)))
            meta.append((noise, plain, order, note))
    t0 = time.time()
    res = run_family(pool, f, items, P)
    dt = time.time() - t0
    for k, r in enumerate(res):
        sm = summarize(r)
        if kind == 'real':
            sm['all'] = [dict(score=round(x[0], 5), text=x[1], order=x[2], note=x[3])
                         for x in sorted(r, key=lambda x: -x[0])]
            st['real'].append(sm)
        elif kind == 'controls':
            noise, plain, order, note = meta[k]
            acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
            st['controls'].append(dict(stat=sm['stat'], noise=noise, found=int(sm['order'] == order),
                                       acc=round(acc, 3), text=sm['text'], plain=plain, note=note))
        else:
            st[kind].append(dict(stat=sm['stat'], order=sm['order']))
    state['unit_secs'][f] = dt / count
    return dt


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def verdict(p_d, power):
    if p_d is None:
        return 'waiting for copies'
    if p_d <= BONF:
        return 'BEATS the copies past Bonferroni: read the decrypt by eye'
    if p_d <= 0.05:
        return 'below 0.05, not past Bonferroni (expected by chance in 3 families)'
    if power is not None and power < 0.5:
        return 'chance level, LOW POWER: inconclusive'
    return 'chance level'


def segment(text):
    import Dorabella_date as D
    return D.segment(text)


def build_report(state):
    L = ["DORABELLA HOMOPHONES REPORT  (several symbols may stand for one letter)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing "
         f"{state['params']['restarts']} restarts x {state['params']['iters']} steps per reading order",
         f"cipher: {state.get('source', '')}",
         f"Bonferroni line for one of 3 families: p < {BONF:.4f}", "",
         "=" * 104,
         "HEADLINE (statistic = best quadgram score over the 5 reading orders; real = median of seeds;",
         "          power = planted homophonic controls above the 95th percentile of the pair-preserving copies)",
         "=" * 104,
         f"{'fam':<4} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} {'p_shuf':>7} {'power':>6}  verdict"]
    detail = []
    for f in FAMILIES:
        st = state['fam'][f]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        nd = [x['stat'] for x in st['doublet']]
        ns = [x['stat'] for x in st['shuffle']]
        p_d = Z.pval(med, nd) if med is not None and nd else None
        p_s = Z.pval(med, ns) if med is not None and ns else None
        q95 = Z.pct(nd, 0.95)
        ctrl = st['controls']
        power = (sum(1 for c in ctrl if c['stat'] > q95) / len(ctrl)) if (ctrl and q95 is not None) else None
        L.append(f"{f:<4} {fmt(med):>8} {len(nd):>6} {fmt(p_d, 4):>7} {len(ns):>6} {fmt(p_s, 4):>7} "
                 f"{('-' if power is None else f'{power:.0%}'):>6}  {verdict(p_d, power)}")
        detail.append((f, st, med, nd, ns, q95))
    L += ["", "A family matters only if the real cipher beats the pair-preserving copies past the",
          "Bonferroni line, its controls show it can find a message, and the decrypt reads as English."]
    for f, st, med, nd, ns, q95 in detail:
        L += ["", "-" * 104, f"{f}: {FAMILY_DESC[f]}", "-" * 104]
        if not st['real']:
            L.append("real run not done yet")
            continue
        L.append("real statistic per seed: " + ', '.join(f"{r['stat']:.3f}" for r in st['real'])
                 + f"   median {med:.3f}")
        for name, nul in (('pair-preserving', nd), ('plain shuffle', ns)):
            if nul:
                L.append(f"  vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{Z.pct(nul, .95):.3f}, max {max(nul):.3f}, p = {Z.pval(med, nul):.4f}")
        best = max(st['real'], key=lambda r: r['stat'])
        L.append("best real decrypts, one per reading order (seed with the highest statistic):")
        for t in best['all']:
            L.append(f"  [{t['score']:.3f}] {t['order']}   letters with several symbols: {t['note']}")
            for a in range(0, 87, 29):
                L.append("      " + t['text'][a:a + 29])
            L.append("      words: " + segment(t['text']))
        if st['controls']:
            for nz in sorted({c['noise'] for c in st['controls']}):
                cs = [c for c in st['controls'] if c['noise'] == nz]
                above = sum(1 for c in cs if q95 is not None and c['stat'] > q95)
                L.append(f"controls, {nz} misread symbols: n={len(cs)}, median stat "
                         f"{statistics.median(c['stat'] for c in cs):.3f}, reading order found "
                         f"{sum(c['found'] for c in cs)}/{len(cs)}, median letters right "
                         f"{statistics.median(c['acc'] for c in cs):.0%}, above null 95th pct "
                         f"{above}/{len(cs) if q95 is not None else '-'}")
            ex = max(st['controls'], key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right; homophones {ex['note']}):")
            L.append(f"      planted {ex['plain'][:58]}")
            L.append(f"      found   {ex['text'][:58]}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, P, SEQ, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    used = sorted({Z.slot_of(x) for x in SEQ})
    print(f"the cipher uses {len(used)} different symbols of 24")
    for f in FAMILIES:
        print(f"  {f}: {FAMILY_DESC[f]}  (letter limits sum to {sum(CAPS[f])})")
    # annealer sanity: a one-to-one key found by H2 must score at least as well as the masc search
    print("\nOne planted homophonic-English control per family, through the full family search:")
    est = {}
    for f in (args.only or FAMILIES):
        seq, plain, order, note = plant(f, 3, 0, Z.derive_seed('selftest', f))
        t0 = time.time()
        sm = summarize(run_family(pool, f, [(seq, 11)], P)[0])
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
        s = M.doublet_shuffle(list(SEQ), random.Random(Z.derive_seed('selftest-null', f)))
        nm = summarize(run_family(pool, f, [(s, 5)], P)[0])
        est[f] = dt
        print(f"  {f}: control {sm['stat']:.3f} (a pair-preserving copy of the cipher: {nm['stat']:.3f}), "
              f"reading order found: {'yes' if sm['order'] == order else 'no'}, letters right {acc:.0%}, "
              f"{H.fmt_secs(dt)} per search")
        print(f"      homophones planted: {note}")
        print(f"      planted: {plain[:58]}")
        print(f"      decoded: {sm['text'][:58]}")
    tot = 0.0
    print("\nEstimated time to reach the targets:")
    for f in (args.only or FAMILIES):
        n = sum(DEFAULT_TARGETS[f].values()) + args.real_seeds
        tot += n * est[f]
        print(f"  {f}: {n} searches x {H.fmt_secs(est[f])} = {H.fmt_secs(n * est[f])}")
    print(f"  total about {H.fmt_secs(tot)}")


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='homophones_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1920)
    ap.add_argument('--real-seeds', type=int, default=5)
    ap.add_argument('--restarts', type=int, default=6)
    ap.add_argument('--iters', type=int, default=60000)
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
    if os.path.exists(sp) and not args.fresh:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state.get('seq') not in (None, SEQ):
            raise SystemExit("saved state was made with a DIFFERENT transcription; use --fresh or another --out")
    else:
        state = blank_state(args)
    state['source'], state['seq'] = source, SEQ
    P = state['params']
    if args.report_only:
        Z.write_text(rp, build_report(state))
        print(f"report written to {rp}")
        return
    workers = max(1, args.workers)
    Z.keep_awake()
    pool = mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, P, SEQ, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}, cipher from {source}")
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
                    for kind in ('shuffle', 'doublet'):
                        state['targets'][f][kind] = int(state['targets'][f][kind] * 1.5)
                log("all targets met; raising targets by 50% and carrying on")
                continue
            _, _, f, kind = stg[0]
            n = batch_size(state, f, kind, workers)
            dt = do_batch(pool, state, SEQ, f, kind, n)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['fam'][f]
            msg = f"{f} {kind} x{n} in {H.fmt_secs(dt)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][f][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{Z.stamp()} {msg}", flush=True)
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            Z.save_json(sp, state)
            Z.write_text(rp, build_report(state))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    mp.freeze_support()
    main()
