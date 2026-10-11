#!/usr/bin/env python3
"""
Dorabella CIRCLE-AND-CROSS KEYS test: Elgar's 24 symbols laid out on a circle with a cross.

Needs Dorabella_three.py and what it needs, Dorabella_transcription.csv, english.txt and numpy.

  python Dorabella_circlekeys.py --selftest          # checks + one planted control + time estimate
  python Dorabella_circlekeys.py --stop-when-done    # results in circlekeys_out\\
  python Dorabella_circlekeys.py --report-only

THE IDEA (from Elgar's 1920 notebook page, as redrawn by the investigator)
  A circle with a cross gives "eight sets of three": three ticks on each of the 4 arms
  (facing up, down, left or right) and three strokes on the rim in each of the 4 quarters
  (facing the 4 in-between directions), each set with 1, 2 and 3 loops. That is exactly
  Elgar's 24 symbols. His own published key fills ONE cross first (ABC, DEF, GHI, JKL on
  directions 3, 5, 1, 7, all 90 degrees apart) and then the OTHER cross (MNO, PQR, STU, WXY
  on 8, 6, 4, 2) - an arms-then-rim layout, not a walk round eight spokes.
  Which side of an arm the ticks sit on (left/right, up/down) and whether the rim strokes point
  in or out decide which way each set faces; and each set of three letters can run "up"
  (1 loop = first letter) or "down" (1 loop = last letter). If Elgar drew a different picture
  for a different cipher, each picture is a different key. This test tries EVERY such picture.

PRE-REGISTERED FAMILY (frozen 9 Oct 2026, before any run; one family, so p < 0.05)
  CK  every circle-and-cross key, no search:
        letter groups, 2 ways   Elgar's own (ABC DEF GHI JKL MNO PQR STU WXY, J read as I)
                                or the plain alphabet (ABC DEF GHI KLM NOP QRS TUW XYZ)
        which cross the first four groups go on, 2 ways (the other four go on the other cross)
        the order of the first four groups round their cross, 4! = 24 ways
        the order of the last four groups round theirs, 4! = 24 ways
                                (every left/right, up/down and in/out choice of the drawing
                                is one of these 24 x 24 orders)
        each of the 8 groups "up" or "down", 2^8 = 256 ways
      = 589,824 keys (Elgar's own published key is one of them), each read in the 5 reading
      orders (forward, fully reversed, each line reversed, two boustrophedons): about 2.9 million
      readings per copy. Statistic = best quadgram score; symbols read as in the transcription.
  Compared with pair-preserving copies (main) and plain shuffles put through all the same keys.
  Planted controls: English written with a random circle key and reading order (0 or 2 misreads).
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
FAMILIES = ['CK']
FAMILY_DESC = {'CK': 'every circle-and-cross key (589,824 keys x 5 reading orders), no search'}
BONF = 0.05
ORDERS = Z.KNOWN_ORDERS
DEFAULT_TARGETS = {'CK': dict(shuffle=1000, doublet=2000, controls=40)}
GROUPS = [('Elgar groups', ['ABC', 'DEF', 'GHI', 'JKL', 'MNO', 'PQR', 'STU', 'WXY']),
          ('alphabet groups', ['ABC', 'DEF', 'GHI', 'KLM', 'NOP', 'QRS', 'TUW', 'XYZ'])]
CROSSES = [('first four on the odd cross (1,3,5,7)', [1, 3, 5, 7], [2, 4, 6, 8]),
           ('first four on the even cross (2,4,6,8)', [2, 4, 6, 8], [1, 3, 5, 7])]
PERMS = list(itertools.permutations(range(4)))          # 24
NKEYS = len(GROUPS) * len(CROSSES) * 24 * 24 * 256
CHUNK = 16384
KEYS = None
ORDER_R = None


def decode_index(ix):
    ix, bits = divmod(ix, 256)
    ix, pb = divmod(ix, 24)
    ix, pa = divmod(ix, 24)
    g, c = divmod(ix, 2)
    return g, c, pa, pb, bits


def key_parts(ix):
    """[(letters, direction, '^' or 'v')] for the 8 groups of key number ix."""
    g, c, pa, pb, bits = decode_index(ix)
    groups = GROUPS[g][1]
    _, ca, cb = CROSSES[c]
    out = []
    for k in range(4):
        out.append((groups[k], ca[PERMS[pa][k]], '^' if bits >> k & 1 else 'v'))
    for k in range(4):
        out.append((groups[4 + k], cb[PERMS[pb][k]], '^' if bits >> (4 + k) & 1 else 'v'))
    return out


def key_of(ix):
    E = [0] * 24
    for letters, d, ud in key_parts(ix):
        for loop in (1, 2, 3):
            ch = letters[loop - 1] if ud == '^' else letters[3 - loop]
            E[(loop - 1) * 8 + (d - 1)] = H.LIDX[H.fold_base(ch)]
    return E


def label(ix):
    g, c, pa, pb, bits = decode_index(ix)
    parts = key_parts(ix)
    lay = ' '.join(f"{l}@{d}{u}" for l, d, u in parts)
    tag = '  <- ELGAR\'S OWN KEY' if key_of(ix) == Z.ELGAR else ''
    return f"{GROUPS[g][0]} | {lay}{tag}"


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
    ix = rng.randrange(NKEYS)
    E = key_of(ix)
    inv = [[s for s in range(24) if E[s] == x] for x in range(24)]
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
                fam={'CK': dict(real=[], shuffle=[], doublet=[], controls=[])}, log=[])


def stages(state):
    st, tg = state['fam']['CK'], state['targets']['CK']
    if not st['real']:
        return [(0, -1.0, 'CK', 'real')]
    out = []
    if len(st['controls']) < min(10, tg['controls']):
        out.append((1, len(st['controls']) / tg['controls'], 'CK', 'controls'))
    for kind in ('doublet', 'shuffle', 'controls'):
        frac = len(st[kind]) / tg[kind]
        if frac < 1:
            out.append((2, frac, 'CK', kind))
    return sorted(out)


def do_batch(pool, state, SEQ, kind, count):
    st, seed0 = state['fam']['CK'], state['params']['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append(SEQ)
        elif kind in ('shuffle', 'doublet'):
            rng = random.Random(Z.derive_seed(seed0, 'null', 'CK', kind, n))
            if kind == 'shuffle':
                s = list(SEQ)
                rng.shuffle(s)
            else:
                s = M.doublet_shuffle(list(SEQ), rng)
            items.append(s)
        else:
            noise = 0 if n % 2 == 0 else 2
            seq, plain, ix, order = plant(n, noise, Z.derive_seed(seed0, 'plant', 'CK', n))
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
    st = state['fam']['CK']
    L = ["DORABELLA CIRCLE-AND-CROSS KEYS REPORT  (every 'eight sets of three' picture as the key)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"keys: {NKEYS:,} x {len(ORDERS)} reading orders; no search; cipher: {state.get('source', '')}",
         "one pre-registered family, so the line is p < 0.05", "",
         "Directions are numbered 1-8 in 45-degree steps as in all the other scripts; 'ABC@3v' means",
         "the letters ABC sit on direction 3 running 'down' (1 loop = C, 3 loops = A); '^' = 'up'.", ""]
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
          f"{'CK':<4} {fmt(med):>8} {len(nd):>6} {fmt(p_d, 4):>7} {len(ns):>6} {fmt(p_s, 4):>7} "
          f"{('-' if power is None else f'{power:.0%}'):>6}  {verdict(p_d, power)}",
          "=" * 104]
    if st['real']:
        for name, nul in (('pair-preserving', nd), ('plain shuffle', ns)):
            if nul:
                L.append(f"  vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{Z.pct(nul, .95):.3f}, max {max(nul):.3f}, p = {Z.pval(med, nul):.4f}")
        elgar = state.get('elgar_forward')
        if elgar is not None:
            L.append(f"  Elgar's own key, read forwards, scores {elgar:.3f} (it is one of the {NKEYS:,} keys)")
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
    elg = [ix for ix in range(NKEYS) if key_of(ix) == Z.ELGAR]
    print(f"{NKEYS:,} circle-and-cross keys built; Elgar's own key is key number(s) {elg}")
    assert elg, "Elgar's key should be one of the circle keys"
    print(f"  {label(elg[0])}")
    distinct = len({tuple(r) for r in K[:: 97].tolist()})
    print(f"  spot check: {distinct} distinct keys in a sample of {len(K[::97])}")
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
    n = sum(DEFAULT_TARGETS['CK'].values()) + 1
    print(f"\nEstimated time to reach the targets: {n} copies x {H.fmt_secs(dt)} = {H.fmt_secs(n * dt)}")


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='circlekeys_out')
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
                    state['targets']['CK'][kind] = int(state['targets']['CK'][kind] * 1.5)
                log("all targets met; raising targets by 50% and carrying on")
                continue
            _, _, f, kind = stg[0]
            tg = state['targets']['CK']
            n = 1 if kind == 'real' else max(1, min(workers, tg[kind] - len(state['fam']['CK'][kind])))
            dt = do_batch(pool, state, SEQ, kind, n)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['fam']['CK']
            msg = f"CK {kind} x{n} in {H.fmt_secs(dt)}"
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
