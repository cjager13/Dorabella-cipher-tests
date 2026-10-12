#!/usr/bin/env python3
"""
Dorabella TICK-SHIFT test: Elgar's own key, with each group of three letters shifted up or down.

Needs Dorabella_three.py and what it needs, Dorabella_transcription.csv, english.txt and numpy.

  python Dorabella_tickshift.py --selftest          # checks + one planted control + time estimate
  python Dorabella_tickshift.py --stop-when-done    # results in tickshift_out\\
  python Dorabella_tickshift.py --report-only

THE IDEA (the investigator's, from Elgar's 1920 "circle and cross" drawings)
  Each group of three letters sits on an arc or an arm of a circle with a cross, and its
  three symbols are attached by ticks. The symbols themselves never change: each is drawn
  exactly as in Elgar's key. But the ticks of a group point either OUT (each letter is
  moved forward) or IN (each letter is moved back), and the move is either 1 letter or the
  symbol's loop count (1, 2 or 3). All three ticks of a group point the same way.
  Example: KLM ticks out, move 1 -> the K, L, M symbols mean L, M, N;
           KLM ticks in, move by loops -> K-1, L-2, M-3 = I, I, I (letters wrap round A..Z, which has no J or V).

PRE-REGISTERED FAMILY (frozen 11 Oct 2026, before any run; one family, so p < 0.05)
  TS  every tick-shift key, no search:
        letter groups, 2 ways   Elgar's own key (ABC DEF GHI JKL MNO PQR STU WXY, J read as I,
                                each group on Elgar's direction and running up or down as in
                                his key), or the plain alphabet groups (ABC DEF GHI KLM NOP QRS
                                TUW XYZ) on the same directions, up or down, in the same order
        ticks in or out for each of the 8 groups independently: 2^8 = 256
        move 1 or move by the loop count (the same rule for every group): 2
      = 1,024 keys, each read in the 5 reading orders (forward, fully reversed, each line
      reversed, two boustrophedons). A key may send two symbols to the same letter; that is
      allowed (it is how such a key would really behave). Statistic = best quadgram score;
      symbols read as in the transcription.
  Compared with pair-preserving copies (main) and plain shuffles put through all the same keys.
  Planted controls: English written with a random tick-shift key and reading order (0 or 2
  misreads); letters a key cannot write are avoided by choosing another key for that control.
"""
import argparse
import itertools
import json
import math
import multiprocessing as mp
import os
import random
import signal
import statistics
import time

import numpy as np

import Dorabella_three as Z

H, M = Z.H, Z.M
VERSION = 1
FAMILIES = ['TS']
FAMILY_DESC = {'TS': 'every tick-shift key on Elgar\'s key (1,024 keys x 5 reading orders), no search'}
BONF = 0.05
ORDERS = Z.KNOWN_ORDERS
DEFAULT_TARGETS = {'TS': dict(shuffle=2000, doublet=5000, controls=100)}
GROUPS = [('Elgar groups', ['ABC', 'DEF', 'GHI', 'JKL', 'MNO', 'PQR', 'STU', 'WXY']),
          ('alphabet groups', ['ABC', 'DEF', 'GHI', 'KLM', 'NOP', 'QRS', 'TUW', 'XYZ'])]
PLACES = [(d, ud) for _, d, ud in Z.ELGAR_SETS]       # Elgar's direction and up/down for each group
STEPS = ['move 1', 'move by loop count']
NKEYS = len(GROUPS) * 256 * len(STEPS)
CHUNK = 16384
KEYS = None
ORDER_R = None


def decode_index(ix):
    ix, step = divmod(ix, len(STEPS))
    g, bits = divmod(ix, 256)
    return g, bits, step


def key_of(ix):
    g, bits, step = decode_index(ix)
    E = [0] * 24
    for k, (letters, (d, ud)) in enumerate(zip(GROUPS[g][1], PLACES)):
        sign = 1 if bits >> k & 1 else -1
        for loop in (1, 2, 3):
            ch = letters[loop - 1] if ud == '^' else letters[3 - loop]
            base = H.LIDX[H.fold_base(ch)]
            E[(loop - 1) * 8 + (d - 1)] = (base + sign * (1 if step == 0 else loop)) % 24
    return E


def label(ix):
    g, bits, step = decode_index(ix)
    ticks = ' '.join(f"{letters}{'+' if bits >> k & 1 else '-'}" for k, letters in enumerate(GROUPS[g][1]))
    return f"{GROUPS[g][0]} | {STEPS[step]} | ticks {ticks}"


def all_keys():
    global KEYS
    if KEYS is None:
        KEYS = np.array([key_of(ix) for ix in range(NKEYS)], dtype=np.int32)
    return KEYS


def order_perms():
    global ORDER_R
    if ORDER_R is None:
        ORDER_R = [Z.known_order(o) for o in ORDERS]
    return ORDER_R


# ============================================================== scoring
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    Z.init_tables(corpus)
    all_keys()


def score_job(seq, oi, top=5):
    """Every key on one reading order. Returns the top (score, key index) pairs."""
    K = all_keys()
    QT = Z.QT
    R = order_perms()[oi]
    slots = np.array([Z.slot_of(seq[r]) for r in R], dtype=np.intp)
    best_s, best_i = [], []
    for a in range(0, NKEYS, CHUNK):
        P = K[a:a + CHUNK][:, slots]
        idx = ((P[:, :-3] * 24 + P[:, 1:-2]) * 24 + P[:, 2:-1]) * 24 + P[:, 3:]
        sc = QT[idx].sum(axis=1)
        k = min(top, len(sc))
        part = np.argpartition(-sc, k - 1)[:k]
        best_s.extend(sc[part].tolist())
        best_i.extend((part + a).tolist())
    order = sorted(range(len(best_s)), key=lambda j: -best_s[j])[:top]
    n = len(slots) - 3
    return [(best_s[j] / n, best_i[j], oi) for j in order]


def run_family(pool, items):
    """items: list of seqs. Returns per item the merged top list [(score, key, order index)]."""
    jobs, owner = [], []
    for k, seq in enumerate(items):
        for oi in range(len(ORDERS)):
            jobs.append((seq, oi))
            owner.append(k)
    res = Z.run_jobs(pool, score_job, jobs)
    out = [[] for _ in items]
    for k, r in zip(owner, res):
        out[k].extend(r)
    return [sorted(o, key=lambda x: -x[0]) for o in out]


def decrypt(seq, ix, oi):
    E = key_of(ix)
    R = order_perms()[oi]
    return Z.text_of([E[Z.slot_of(seq[r])] for r in R])


def summarize(seq, top):
    s, ix, oi = top[0]
    return dict(stat=round(s, 5), key=ix, order=ORDERS[oi], text=decrypt(seq, ix, oi))


# ============================================================== controls
def plant(idx, noise, seed):
    rng = random.Random(seed)
    texts = M.CONTROL_TEXTS
    raw = texts[idx % len(texts)] + ' ' + texts[(idx + 7) % len(texts)]
    plain = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))[:87]
    L = [H.LIDX[ch] for ch in plain]
    for _try in range(1000):
        ix = rng.randrange(NKEYS)
        E = key_of(ix)
        inv = [[s for s in range(24) if E[s] == x] for x in range(24)]
        if all(inv[x] for x in set(L)):
            break
    oi = rng.randrange(len(ORDERS))
    R = order_perms()[oi]
    seq = [0] * 87
    for i, r in enumerate(R):
        seq[r] = Z.code_of_slot(rng.choice(inv[L[i]]) if inv[L[i]] else rng.randrange(24))
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, ix, ORDERS[oi]


# ============================================================== state and batches
def blank_state(args):
    return dict(version=VERSION, created=str(Z.now()), elapsed=0.0, sessions=0,
                params=dict(seed=args.seed), targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                fam={'TS': dict(real=[], shuffle=[], doublet=[], controls=[])}, log=[])


def stages(state):
    st, tg = state['fam']['TS'], state['targets']['TS']
    if not st['real']:
        return [(0, -1.0, 'TS', 'real')]
    out = []
    if len(st['controls']) < min(10, tg['controls']):
        out.append((1, len(st['controls']) / tg['controls'], 'TS', 'controls'))
    for kind in ('doublet', 'shuffle', 'controls'):
        frac = len(st[kind]) / tg[kind]
        if frac < 1:
            out.append((2, frac, 'TS', kind))
    return sorted(out)


def do_batch(pool, state, SEQ, kind, count):
    st, seed0 = state['fam']['TS'], state['params']['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append(SEQ)
        elif kind in ('shuffle', 'doublet'):
            rng = random.Random(Z.derive_seed(seed0, 'null', 'TS', kind, n))
            if kind == 'shuffle':
                s = list(SEQ)
                rng.shuffle(s)
            else:
                s = M.doublet_shuffle(list(SEQ), rng)
            items.append(s)
        else:
            noise = 0 if n % 2 == 0 else 2
            seq, plain, ix, order = plant(n, noise, Z.derive_seed(seed0, 'plant', 'TS', n))
            items.append(seq)
            meta.append((noise, plain, ix, order))
    t0 = time.time()
    res = run_family(pool, items)
    dt = time.time() - t0
    for k, top in enumerate(res):
        sm = summarize(items[k], top)
        if kind == 'real':
            sm['top'] = [dict(score=round(s, 5), key=ix, order=ORDERS[oi], text=decrypt(SEQ, ix, oi))
                         for s, ix, oi in top[:10]]
            # where does Elgar's own key (forward reading) stand?
            st['real'].append(sm)
        elif kind == 'controls':
            noise, plain, ix, order = meta[k]
            acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
            st['controls'].append(dict(stat=sm['stat'], noise=noise,
                                       found=int(sm['key'] == ix and sm['order'] == order),
                                       acc=round(acc, 3), text=sm['text'], plain=plain))
        else:
            st[kind].append(dict(stat=sm['stat'], key=sm['key'], order=sm['order']))
    return dt


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def verdict(p_d, power):
    if p_d is None:
        return 'waiting for copies'
    if p_d <= BONF:
        return 'BEATS the copies (p < 0.05): read the decrypt by eye'
    if power is not None and power < 0.5:
        return 'chance level, LOW POWER: inconclusive'
    return 'chance level'


def segment(text):
    import Dorabella_date as D
    return D.segment(text)


def build_report(state):
    st = state['fam']['TS']
    L = ["DORABELLA TICK-SHIFT REPORT  (Elgar's key, each group moved up or down by its ticks)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"keys: {NKEYS:,} x {len(ORDERS)} reading orders; no search; cipher: {state.get('source', '')}",
         "one pre-registered family, so the line is p < 0.05", "",
         "Labels: 'KLM+' = that group's ticks point out (letters moved forward), 'KLM-' = in (moved back).", ""]
    real = [r['stat'] for r in st['real']]
    med = real[0] if real else None
    nd = [x['stat'] for x in st['doublet']]
    ns = [x['stat'] for x in st['shuffle']]
    p_d = Z.pval(med, nd) if med is not None and nd else None
    p_s = Z.pval(med, ns) if med is not None and ns else None
    q95 = Z.pct(nd, 0.95)
    ctrl = st['controls']
    power = (sum(1 for c in ctrl if c['stat'] > q95) / len(ctrl)) if (ctrl and q95 is not None) else None
    L += ["=" * 104,
          f"{'fam':<4} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} {'p_shuf':>7} {'power':>6}  verdict",
          f"{'TS':<4} {fmt(med):>8} {len(nd):>6} {fmt(p_d, 4):>7} {len(ns):>6} {fmt(p_s, 4):>7} "
          f"{('-' if power is None else f'{power:.0%}'):>6}  {verdict(p_d, power)}",
          "=" * 104]
    if st['real']:
        for name, nul in (('pair-preserving', nd), ('plain shuffle', ns)):
            if nul:
                L.append(f"  vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{Z.pct(nul, .95):.3f}, max {max(nul):.3f}, p = {Z.pval(med, nul):.4f}")
        elgar = state.get('elgar_forward')
        if elgar is not None:
            L.append(f"  for reference: Elgar's own key unshifted, read forwards, scores {elgar:.3f}")
        L.append("")
        L.append("best real decrypts:")
        for t in st['real'][0]['top']:
            L.append(f"  [{t['score']:.3f}] {t['order']} | {label(t['key'])}")
            for a in range(0, 87, 29):
                L.append("      " + t['text'][a:a + 29])
            L.append("      words: " + segment(t['text']))
    if ctrl:
        L.append("")
        for nz in sorted({c['noise'] for c in ctrl}):
            cs = [c for c in ctrl if c['noise'] == nz]
            above = sum(1 for c in cs if q95 is not None and c['stat'] > q95)
            L.append(f"controls, {nz} misread symbols: n={len(cs)}, median stat "
                     f"{statistics.median(c['stat'] for c in cs):.3f}, true key and order found "
                     f"{sum(c['found'] for c in cs)}/{len(cs)}, median letters right "
                     f"{statistics.median(c['acc'] for c in cs):.0%}, above null 95th pct "
                     f"{above}/{len(cs) if q95 is not None else '-'}")
        ex = max(ctrl, key=lambda c: c['acc'])
        L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:56]}")
        L.append(f"                                   found   {ex['text'][:56]}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(pool, SEQ, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    K = all_keys()
    print(f"{NKEYS:,} tick-shift keys built; {len({tuple(r) for r in K.tolist()}):,} distinct")
    # check: shifting Elgar's key back by the same rule gives Elgar's key
    ex = 255 * 2          # Elgar groups, all ticks out, move 1
    back = [(x - 1) % 24 for x in key_of(ex)]
    assert back == Z.ELGAR, "all-out move-1 key should be Elgar's key moved forward by one"
    print(f"  check: '{label(ex)}' is Elgar's key moved forward one letter: OK")
    one2one = sum(1 for r in K.tolist() if len(set(r)) == 24)
    print(f"  keys giving each symbol a different letter: {one2one:,}; keys with shared letters: {NKEYS - one2one:,}")
    seq, plain, ix, order = plant(3, 0, Z.derive_seed('selftest'))
    t0 = time.time()
    top = run_family(pool, [seq])[0]
    dt = time.time() - t0
    sm = summarize(seq, top)
    acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
    s = M.doublet_shuffle(list(SEQ), random.Random(Z.derive_seed('selftest-null')))
    nm = summarize(s, run_family(pool, [s])[0])
    print(f"\nplanted control: {sm['stat']:.3f} (a pair-preserving copy of the cipher: {nm['stat']:.3f}), "
          f"true key and order found: {'yes' if sm['key'] == ix and sm['order'] == order else 'no'}, "
          f"letters right {acc:.0%}, {H.fmt_secs(dt)} per copy")
    print(f"  true:    {order} | {label(ix)}")
    print(f"  found:   {sm['order']} | {label(sm['key'])}")
    print(f"  planted: {plain[:58]}")
    print(f"  decoded: {sm['text'][:58]}")
    n = sum(DEFAULT_TARGETS['TS'].values()) + 1
    print(f"\nEstimated time to reach the targets: {n} copies x {H.fmt_secs(dt)} = {H.fmt_secs(n * dt)}")


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='tickshift_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1920)
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
    state['elgar_forward'] = round(Z.score_letters([Z.ELGAR[Z.slot_of(c)] for c in SEQ]), 5)
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
            selftest(pool, SEQ, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}, cipher from {source}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stg = stages(state)
            if not stg:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for kind in ('shuffle', 'doublet'):
                    state['targets']['TS'][kind] = int(state['targets']['TS'][kind] * 1.5)
                log("all targets met; raising targets by 50% and carrying on")
                continue
            _, _, f, kind = stg[0]
            tg = state['targets']['TS']
            n = 1 if kind == 'real' else max(1, min(workers, tg[kind] - len(state['fam']['TS'][kind])))
            dt = do_batch(pool, state, SEQ, kind, n)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['fam']['TS']
            msg = f"TS {kind} x{n} in {H.fmt_secs(dt)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{tg[kind]})")
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
