#!/usr/bin/env python3
"""
Dorabella TWINS test (roadmap items H2 and H3): could the half-turned neighbours come from the
way letters were written, rather than from the hand?

Needs Dorabella_notebook.py, Dorabella_three.py and everything they need (Dorabella_hillclimb.py,
Dorabella_tests.py, Dorabella_marathon.py, Dorabella_transcription.csv, english.txt), plus numpy.

  python Dorabella_twins.py --selftest          # checks, one control per family, time estimate
  python Dorabella_twins.py --stop-when-done    # results in twins_out\\
  python Dorabella_twins.py --report-only
  python Dorabella_twins.py --only H3           # work on one family

Ctrl+C any time; the same command resumes.

TWINS
A symbol's TWIN is the same shape turned half round: the same loop count, facing the opposite
way. The Dorabella has 13 neighbouring pairs that are twins (English through any key: about 4;
through Elgar's notebook key: about 6).

PRE-REGISTERED (frozen 19 Nov 2026, before any run; two families, so the line is p < 0.025)
  H2  A TWIN MARKS A DOUBLED LETTER (Hauer et al. 2021). Wherever a symbol is followed by its twin,
      the twin stands for the same letter again (LL, EE, SS written so that no symbol is drawn
      twice in a row). The rule is applied in each of the 5 reading orders (a twin takes the
      letter of the symbol before it), then ANY one-symbol-one-letter key is searched by
      annealing. Statistic: best quadgram score over the 5 orders; real = median of 5 seeds.
      Expected to fail: 13 twin pairs against about 3 doubled letters in 87 letters of English.
  H3  TWINS MAY SHARE A LETTER (shape homophones). A letter may be written either as its symbol or
      as the twin, chosen for looks; the letters that lose a symbol this way are rare ones that
      the message does without. The search is ANY key in which up to 6 of the 12 twin pairs may
      share one letter (which pairs, and whether, is part of the search), by annealing over the 5
      reading orders. Statistic: best quadgram score; real = median of 5 seeds.
  Both: compared with pair-preserving copies (main; they keep exactly the same 13 twin pairs) and
  plain shuffles, through the identical search. Power from planted English written that way:
  H2 with doubled letters written as a symbol and its twin; H3 with 1 to 6 twin pairs sharing a
  letter (the rarest letters left out of the message) and the twin chosen to face the previous
  symbol whenever it can (so the plants have the habit too); 0 or 2 misread symbols.
  A pass also needs a decrypt that reads as English.
  Transcription: the consensus reading (Hauer et al. 2025; symbols 23, 32, 34, 85, 86 read 14, 11,
  12, 14, 12), as in Dorabella_notebook.py; fixed, because which neighbours are twins depends on it.
"""
import argparse
import json
import math
import os
import random
import signal
import statistics
import time

import Dorabella_notebook as NB          # consensus reading (and the notebook key, not used here)

Z = NB.Z
H, M = Z.H, Z.M
VERSION = 1
FAMILIES = ['H2', 'H3']
FAMILY_DESC = {'H2': 'a twin after its symbol is the same letter again (doubled letters), any key',
               'H3': 'up to 6 twin pairs may share a letter (shape homophones), any key'}
BONF = 0.05 / 2
DEFAULT_TARGETS = {'H2': dict(shuffle=500, doublet=1000, controls=40),
                   'H3': dict(shuffle=500, doublet=1000, controls=40)}
RARITY = 'ETAOINSHRDLUCMWFGYPBKXQZ'      # the 24 letters, commonest first (J -> I, V -> U)
SEQ = None


# ============================================================== twins
def twin(s):
    return (s // 8) * 8 + (s % 8 + 4) % 8


def pair_of(s):
    return (s // 8) * 4 + s % 4          # 12 twin pairs


def pair_members(p):
    b, x = divmod(p, 4)
    return b * 8 + x, b * 8 + x + 4


def is_twin(a, b):
    return b == twin(a)


def orders():
    return [(o, Z.known_order(o)) for o in Z.KNOWN_ORDERS]


# ============================================================== H2
def collapse(c):
    """A twin after its symbol becomes that symbol (the same letter again). c: slots in reading order."""
    e = list(c)
    for q in range(1, len(c)):
        if is_twin(c[q - 1], c[q]):
            e[q] = e[q - 1]
    return e


def h2_job(seq, restarts, iters, seed):
    best = (-1e18, '', '')
    for k, (o, R) in enumerate(orders()):
        e = collapse([Z.slot_of(seq[r]) for r in R])
        s, p, _ = Z.anneal('masc', e, [0] * 87, restarts, iters, Z.derive_seed(seed, k))
        if s > best[0]:
            best = (s, Z.text_of(p), o)
    return best


# ============================================================== H3: annealing with twin merges
def twin_anneal(c, restarts, iters, seed, maxm):
    """Any key in which up to maxm twin pairs share one letter. Returns (avg score, letters, merged pairs)."""
    Q = Z.QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(c)
    pos, starts, used, _ = Z._index(c, n, 24)
    present = set(used)
    eligible = [p for p in range(12) if all(m in present for m in pair_members(p))]
    fac = (Z.SA_T1 / Z.SA_T0) ** (1.0 / max(1, iters))

    def qsum(p, ss):
        t = 0.0
        for s in ss:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    allst = range(n - 3)
    best_all = (-1e18, None, None)
    for _ in range(restarts):
        let = list(range(24))
        rng.shuffle(let)
        mer = [False] * 12
        free = []
        p = [let[x] for x in c]
        cur = qsum(p, allst)
        best, bp, bm = cur, p[:], []
        Tm = Z.SA_T0
        for _it in range(iters):
            r = rand()
            if r < 0.75 or (r < 0.85 and not free) or (r >= 0.85 and not eligible):
                s1, s2 = used[rr(len(used))], rr(24)
                g1 = (s1, twin(s1)) if mer[pair_of(s1)] else (s1,)
                g2 = (s2, twin(s2)) if mer[pair_of(s2)] else (s2,)
                if s2 in g1:
                    Tm *= fac
                    continue
                l1, l2 = let[s1], let[s2]
                changes = [(x, l2) for x in g1] + [(x, l1) for x in g2]
                undo_free = None
            elif r < 0.85:
                s1 = used[rr(len(used))]
                g1 = (s1, twin(s1)) if mer[pair_of(s1)] else (s1,)
                k = rr(len(free))
                changes = [(x, free[k]) for x in g1]
                undo_free = ('swapfree', k, let[s1])
            else:
                pp = eligible[rr(len(eligible))]
                a, b = pair_members(pp)
                if not mer[pp]:
                    if sum(mer) >= maxm:
                        Tm *= fac
                        continue
                    changes = [(b, let[a])]
                    undo_free = ('merge', pp, let[b])
                else:
                    if not free:
                        Tm *= fac
                        continue
                    k = rr(len(free))
                    changes = [(b, free[k])]
                    undo_free = ('unmerge', pp, k)
            slots = {x for x, _ in changes}
            ss = sorted({s for x in slots for s in starts[x]})
            old = qsum(p, ss)
            prev = [(x, let[x]) for x in slots]
            for x, l in changes:
                let[x] = l
                for q in pos[x]:
                    p[q] = l
            d = qsum(p, ss) - old
            if d >= 0 or rand() < exp(d / Tm):
                cur += d
                if undo_free:
                    kind = undo_free[0]
                    if kind == 'swapfree':
                        free[undo_free[1]] = undo_free[2]
                    elif kind == 'merge':
                        mer[undo_free[1]] = True
                        free.append(undo_free[2])
                    else:
                        mer[undo_free[1]] = False
                        free.pop(undo_free[2])
                if cur > best:
                    best, bp, bm = cur, p[:], [i for i in range(12) if mer[i]]
            else:
                for x, l in prev:
                    let[x] = l
                    for q in pos[x]:
                        p[q] = l
            Tm *= fac
        if best > best_all[0]:
            best_all = (best, bp, bm)
    return best_all[0] / (n - 3), Z.text_of(best_all[1]), best_all[2]


def h3_job(seq, restarts, iters, seed, maxm):
    best = (-1e18, '', '', [])
    for k, (o, R) in enumerate(orders()):
        c = [Z.slot_of(seq[r]) for r in R]
        s, text, merged = twin_anneal(c, restarts, iters, Z.derive_seed(seed, k), maxm)
        if s > best[0]:
            best = (s, text, o, merged)
    return best


def pair_label(p):
    a, b = pair_members(p)
    return f"{Z.code_of_slot(a)}/{Z.code_of_slot(b)}"


# ============================================================== controls
def english_letters(idx, need, keep=None):
    texts = M.CONTROL_TEXTS
    out, j = '', idx
    while len(out) < need:
        raw = texts[j % len(texts)]
        s = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))
        if keep is not None:
            s = ''.join(ch for ch in s if ch in keep)
        out += s
        j += 7
    return out[:need]


def plant(fam, idx, noise, seed):
    rng = random.Random(seed)
    o, R = orders()[rng.randrange(5)]
    if fam == 'H2':
        plain = english_letters(idx, 87)
        key = list(range(24))
        rng.shuffle(key)
        kinv = Z.inverse(key)
        slots = []
        for i, ch in enumerate(plain):
            L = H.LIDX[ch]
            if i and plain[i - 1] == ch:
                slots.append(twin(slots[-1]))          # doubled letter: the symbol, then its twin
            else:
                slots.append(kinv[L])
        info = dict(order=o)
    else:
        k = rng.randint(1, 6)
        pairs = rng.sample(range(12), k)
        merged_b = {pair_members(p)[1] for p in pairs}
        groups = [s for s in range(24) if s not in merged_b]          # 24 - k groups, named by first member
        letters = list(RARITY[:24 - k])
        rng.shuffle(letters)
        gl = {g: H.LIDX[l] for g, l in zip(groups, letters)}
        plain = english_letters(idx, 87, keep=set(RARITY[:24 - k]))
        inv = {v: g for g, v in gl.items()}
        slots = []
        for ch in plain:
            g = inv[H.LIDX[ch]]
            if pair_of(g) in pairs:                                    # twin chosen for looks
                a, b = pair_members(pair_of(g))
                if slots and twin(slots[-1]) in (a, b):
                    s = twin(slots[-1])
                else:
                    s = rng.choice((a, b))
            else:
                s = g
            slots.append(s)
        info = dict(order=o, merged=sorted(pairs))
    seq = [0] * 87
    for i, r in enumerate(R):
        seq[r] = Z.code_of_slot(slots[i])
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, info


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


def job_args(fam, seq, P, seed):
    if fam == 'H2':
        return (seq, P['restarts'], P['iters'], seed)
    return (seq, P['restarts'], P['iters'], seed, P['maxm'])


def run(pool, fam, argsets):
    return Z.run_jobs(pool, h2_job if fam == 'H2' else h3_job, argsets)


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
            s, plain, info = plant(fam, n, noise, Z.derive_seed(seed0, 'plant', fam, n))
            items.append(s)
            meta.append((noise, plain, info))
        else:
            items.append(null_seq(kind, fam, n, seed0))
            meta.append(None)
    res = run(pool, fam, [job_args(fam, s, P, Z.derive_seed(seed0, fam, kind, len(st[kind]) + k))
                          for k, s in enumerate(items)])
    for r, mt in zip(res, meta):
        s, text, o = r[0], r[1], r[2]
        rec = dict(stat=round(s, 5), order=o)
        if fam == 'H3':
            rec['merged'] = r[3]
        if kind == 'real':
            rec['text'] = text
        elif kind == 'controls':
            noise, plain, info = mt
            rec.update(noise=noise, plain=plain, text=text, right_order=int(o == info['order']),
                       acc=round(sum(a == b for a, b in zip(text, plain)) / 87, 3) if o == info['order'] else 0.0)
            if fam == 'H3':
                rec['true_merged'] = info['merged']
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


def batch_size(state, f, kind, workers):
    if kind == 'real':
        return 1
    return max(1, min(workers, state['targets'][f][kind] - len(state['fam'][f][kind])))


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def build_report(state):
    tw = [(i + 1, SEQ[i], SEQ[i + 1]) for i in range(86) if is_twin(Z.slot_of(SEQ[i]), Z.slot_of(SEQ[i + 1]))]
    L = ["DORABELLA TWINS REPORT  (roadmap H2 and H3: twin symbols as doubled letters or shared letters)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   sessions {state['sessions']}",
         NB.READ_NOTE,
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing {state['params']['restarts']} "
         f"restarts x {state['params']['iters']} steps per reading order, at most {state['params']['maxm']} shared pairs in H3",
         f"two pre-registered families: the line is p < {BONF:.4f}",
         f"twin neighbours in the cipher ({len(tw)}), position: codes  " + ', '.join(f"{p}: {a} {b}" for p, a, b in tw), ""]
    for f in FAMILIES:
        st = state['fam'][f]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        q95 = Z.pct(dn, .95) if dn else None
        cs = st['controls']
        power = sum(1 for c in cs if q95 is not None and c['stat'] > q95) / len(cs) if (cs and q95 is not None) else None
        p_d = Z.pval(med, dn) if dn else None
        L += ["=" * 104, f"{f}  {FAMILY_DESC[f]}", "=" * 104]
        if med is None:
            L.append("real run not done yet")
            L.append("")
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
        L.append("best real decrypts (one per seed, best first):")
        for r in sorted(st['real'], key=lambda r: -r['stat'])[:3]:
            extra = (f" | shared pairs: {', '.join(pair_label(p) for p in r['merged']) or 'none'}" if f == 'H3' else '')
            L.append(f"  [{r['stat']:.3f}] {r['order']}{extra}")
            for a in range(0, 87, 29):
                L.append("      " + r['text'][a:a + 29])
            L.append("      words: " + Z.show_words(r['text']))
        if f == 'H3' and dn:
            L.append(f"shared pairs used by the copies' best keys: mean {statistics.mean(len(x['merged']) for x in st['doublet']):.1f}; "
                     f"by the real cipher's: {statistics.mean(len(r['merged']) for r in st['real']):.1f}")
        for nz in (0, 2):
            c2 = [c for c in cs if c['noise'] == nz]
            if c2:
                L.append(f"controls, {nz} misread symbols: n={len(c2)}, right reading order {sum(c['right_order'] for c in c2)}/{len(c2)}, "
                         f"median letters right {statistics.median(c['acc'] for c in c2):.0%}, above null 95th pct "
                         f"{sum(1 for c in c2 if q95 is not None and c['stat'] > q95)}/{len(c2)}")
        if cs:
            ex = max(cs, key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:56]}")
            L.append(f"                                   found   {ex['text'][:56]}")
        L.append("")
    L += ["A family matters only if the real cipher beats the PAIR-PRESERVING copies past the line, its controls",
          "show it can find a message, and the decrypt reads as English.", "",
          "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    c = [Z.slot_of(x) for x in SEQ]
    tw = sum(1 for a, b in zip(c, c[1:]) if is_twin(a, b))
    print(f"twin neighbours in the cipher (consensus reading): {tw}")
    assert tw == 13
    for s in range(24):
        assert twin(twin(s)) == s and pair_of(s) == pair_of(twin(s))
    P = dict(restarts=args.restarts, iters=args.iters, maxm=args.maxm)
    for f in FAMILIES:
        seq, plain, info = plant(f, 3, 0, Z.derive_seed('selftest', f))
        t0 = time.time()
        r = run(pool, f, [job_args(f, seq, P, 7)])[0]
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(r[1], plain)) / 87 if r[2] == info['order'] else 0.0
        nul = run(pool, f, [job_args(f, null_seq('doublet', f, 0, 12345), P, 8)])[0]
        print(f"\n{f}: planted control {r[0]:.3f} (a pair-preserving copy: {nul[0]:.3f}); reading order "
              f"{'right' if r[2] == info['order'] else 'wrong'}; letters right {acc:.0%}; {H.fmt_secs(dt)} per search")
        if f == 'H3':
            print(f"  shared pairs planted {[pair_label(p) for p in info['merged']]}, found {[pair_label(p) for p in r[3]]}")
        print(f"  planted {plain[:58]}\n  decoded {r[1][:58]}")
        n = sum(DEFAULT_TARGETS[f].values()) + 5
        print(f"  estimated time for {n} searches, {workers} at a time: about {H.fmt_secs(n * dt / workers)}")


# ============================================================== main
def main():
    global SEQ
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='twins_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--restarts', type=int, default=8)
    ap.add_argument('--iters', type=int, default=15000)
    ap.add_argument('--maxm', type=int, default=6, help='H3: most twin pairs that may share a letter')
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
                     params=dict(seed=args.seed, restarts=args.restarts, iters=args.iters, maxm=args.maxm, real_seeds=5),
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
            n = batch_size(state, f, kind, workers)
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
