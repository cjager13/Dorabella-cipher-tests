#!/usr/bin/env python3
"""
Dorabella DATE test: is the date Elgar wrote on the cipher ("July 14 97") its key?

Needs Dorabella_three.py (and what it needs: Dorabella_hillclimb.py, Dorabella_tests.py,
Dorabella_marathon.py), Dorabella_transcription.csv and english.txt in the same folder.

  python Dorabella_date.py --selftest        # checks, one planted message per family, time estimate
  python Dorabella_date.py                   # start, or resume where it left off
  python Dorabella_date.py --stop-when-done  # exit when every target is met
  python Dorabella_date.py --report-only / --fresh / --hours N / --only DA DT

Stop any time with Ctrl+C; the same command resumes. Progress is saved after every batch to
date_out/state.json; date_out/report.txt is rewritten each time.

PRE-REGISTERED KEYS (frozen before any run; every month/day/year arrangement of the date)
  digit keys   71497  7141897  071497  07141897      month day year
               14797  1471897  140797  14071897      day month year (British order)
               970714 18970714                       year month day
               714    147                            month and day only
               195    97195    1897195               day of the year (195), YYDDD date code
  number keys  7,14,97   14,7,97                     the date's numbers used whole

HOW A KEY IS USED (each one is tried every way below; keys run along the cipher in the
order it was written, either straight through or restarting at each line)
  letter shift   each letter moved along the alphabet by the key value (Gronsfeld), + or -
  rotate         each symbol turned by the key value x 45 degrees, + or -
  loop step      each symbol's loop count stepped by the key value (mod 3), + or -
  columnar       the digits give the column order of a transposition (digit keys only),
                 read as decryption and as encryption

FAMILIES (3, so the Bonferroni line for one family is p < 0.05/3 = 0.0167)
  DE  Elgar's own published key, with every date-key setting above (no search needed)
  DA  ANY one-symbol-one-letter key (annealing), with every shift / rotate / loop-step setting
  DT  ANY key (annealing), with every date-key columnar transposition

THREE COMPARISONS for each family (statistic = best score over the family's settings;
real = median over search seeds)
  pair-preserving copies   scrambled cipher keeping every adjacent symbol pair (the main null)
  plain shuffles           scrambled cipher
  random-number keys       the REAL cipher, with every date key replaced by random digits of the
                           same length: do the date's numbers beat any other numbers?
Planted English enciphered with a random date setting shows each family can find a message.
Uncertain symbols are read as in the transcription (as in Dorabella_three.py).
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
FAMILIES = ['DE', 'DA', 'DT']
FAMILY_DESC = {'DE': "Elgar's key + date-key shifts, turns and transpositions (no search)",
               'DA': 'any key + date-key letter shift / rotate / loop step',
               'DT': 'any key + date-key columnar transposition'}
BONF = 0.05 / len(FAMILIES)
DIGIT_KEYS = [('71497', 'M D YY'), ('7141897', 'M D YYYY'), ('071497', 'MM DD YY'),
              ('07141897', 'MM DD YYYY'), ('14797', 'D M YY'), ('1471897', 'D M YYYY'),
              ('140797', 'DD MM YY'), ('14071897', 'DD MM YYYY'), ('970714', 'YY MM DD'),
              ('18970714', 'YYYY MM DD'), ('714', 'M D'), ('147', 'D M'),
              ('195', 'day of year'), ('97195', 'YYDDD'), ('1897195', 'YYYYDDD')]
NUMBER_KEYS = [((7, 14, 97), 'M,D,YY'), ((14, 7, 97), 'D,M,YY')]
LINE_OF = [0] * 29 + [1] * 31 + [2] * 27
DEFAULT_TARGETS = {'DE': dict(shuffle=5000, doublet=5000, randkey=5000, controls=40),
                   'DA': dict(shuffle=300, doublet=600, randkey=300, controls=40),
                   'DT': dict(shuffle=500, doublet=1000, randkey=500, controls=40)}


# ============================================================== members
def phases(n_vals, restart):
    if not restart:
        return [q % n_vals for q in range(87)]
    cnt, out = [0, 0, 0], []
    for q in range(87):
        out.append(cnt[LINE_OF[q]] % n_vals)
        cnt[LINE_OF[q]] += 1
    return out


def turn_index(rot, loop):
    return (rot % 8) * 3 + (loop % 3)          # Z.TURN order: mirror, rotation, loop step


def build_members(digit_keys, number_keys):
    """digit_keys: [(string, label)], number_keys: [(tuple, label)].
    Returns {'sub': [...], 'col': [...]}; each member a dict (label, type, vals, ph | R)."""
    keys = [([int(ch) for ch in s], f"{s} ({lab})") for s, lab in digit_keys]
    keys += [(list(t), f"{','.join(map(str, t))} ({lab})") for t, lab in number_keys]
    sub, seen = [], set()
    for vals, klab in keys:
        for restart in (False, True):
            ph = phases(len(vals), restart)
            plab = 'restarts each line' if restart else 'runs through'
            for mode in ('letter shift', 'rotate', 'loop step'):
                for sign in (1, -1):
                    if mode == 'letter shift':
                        v = [(sign * x) % 24 for x in vals]
                        typ = 'shift'
                    elif mode == 'rotate':
                        v = [turn_index(sign * x, 0) for x in vals]
                        typ = 'turn'
                    else:
                        v = [turn_index(0, sign * x) for x in vals]
                        typ = 'turn'
                    sig = (typ, tuple(v[ph[q]] for q in range(87)))
                    if sig in seen or all(x == 0 for x in v):
                        continue
                    seen.add(sig)
                    sub.append(dict(label=f"{klab} | {mode} {'+' if sign > 0 else '-'} | {plab}",
                                    type=typ, vals=v, ph=ph))
    col, seen = [], set()
    for vals, klab in keys[:len(digit_keys)]:
        w = len(vals)
        order = sorted(range(w), key=lambda i: (vals[i], i))
        tmp = []
        Z.add_both(tmp, f"{klab} | columnar width {w}", Z.columnar_enc(87, w, order))
        for lab, R in tmp:
            if tuple(R) in seen:
                continue
            seen.add(tuple(R))
            col.append(dict(label=lab, type='col', R=R))
    return {'sub': sub, 'col': col}


DATE_MEMBERS = build_members(DIGIT_KEYS, NUMBER_KEYS)


def random_members(seed):
    """The same structure with every date number replaced by random numbers of the same size."""
    rng = random.Random(seed)
    dk = [(''.join(str(rng.randrange(10)) for _ in s), lab) for s, lab in DIGIT_KEYS]
    nk = [(tuple(rng.randrange(1, 100) for _ in t), lab) for t, lab in NUMBER_KEYS]
    return build_members(dk, nk)


def family_members(fam, members):
    if fam == 'DE':
        return members['sub'] + members['col']
    return members['sub'] if fam == 'DA' else members['col']


# ============================================================== searching
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    Z.init_tables(corpus)


def de_job(seq, members):
    """Elgar's key, no search: one (score, text, label) per member."""
    slots = [Z.slot_of(x) for x in seq]
    out = []
    for m in members:
        if m['type'] == 'shift':
            p = [(Z.ELGAR[s] + m['vals'][m['ph'][q]]) % 24 for q, s in enumerate(slots)]
        elif m['type'] == 'turn':
            p = [Z.ELGAR[Z.TURN[m['vals'][m['ph'][q]]][s]] for q, s in enumerate(slots)]
        else:
            p = [Z.ELGAR[slots[r]] for r in m['R']]
        out.append((Z.score_letters(p), Z.text_of(p), m['label']))
    return out


def sa_job(seq, m, restarts, iters, seed):
    slots = [Z.slot_of(x) for x in seq]
    if m['type'] == 'col':
        s, p, _ = Z.anneal('masc', [slots[r] for r in m['R']], [0] * 87, restarts, iters, seed)
    else:
        model = 'shift3' if m['type'] == 'shift' else 'turn3'
        s, p, _ = Z.anneal(model, slots, m['ph'], restarts, iters, seed, fixed=m['vals'])
    return s, Z.text_of(p), m['label']


def run_family(pool, fam, items, P):
    """items: list of (seq, members dict, seed). Returns per item the list of member results."""
    if fam == 'DE':
        args = [(seq, family_members(fam, mem)) for seq, mem, _ in items]
        return Z.run_jobs(pool, de_job, args)
    jobs, owner = [], []
    for k, (seq, mem, seed) in enumerate(items):
        for mi, m in enumerate(family_members(fam, mem)):
            jobs.append((seq, m, P['restarts'], P['iters'], Z.derive_seed(seed, fam, mi)))
            owner.append(k)
    res = Z.run_jobs(pool, sa_job, jobs)
    out = [[] for _ in items]
    for k, r in zip(owner, res):
        out[k].append(r)
    return out


def summarize(members_res):
    b = max(range(len(members_res)), key=lambda i: members_res[i][0])
    s, text, lab = members_res[b]
    return dict(stat=round(s, 5), best=lab, text=text, words=Z.words_of(text))


# ============================================================== controls
def plant(fam, idx, noise, seed):
    rng = random.Random(seed)
    texts = M.CONTROL_TEXTS
    raw = texts[idx % len(texts)] + ' ' + texts[(idx + 7) % len(texts)]
    plain = H.fold_base(''.join(ch for ch in raw.upper() if ch.isalpha()))[:87]
    L = [H.LIDX[ch] for ch in plain]
    mem = family_members(fam, DATE_MEMBERS)
    m = mem[rng.randrange(len(mem))]
    if fam == 'DE':
        def to_slot(x):
            return rng.choice(Z.ELGAR_INV[x]) if Z.ELGAR_INV[x] else rng.randrange(24)
    else:
        k = list(range(24))
        rng.shuffle(k)
        kinv = Z.inverse(k)

        def to_slot(x):
            return kinv[x]
    seq = [0] * 87
    if m['type'] == 'col':
        for i, r in enumerate(m['R']):
            seq[r] = Z.code_of_slot(to_slot(L[i]))
    elif m['type'] == 'shift':
        for q in range(87):
            seq[q] = Z.code_of_slot(to_slot((L[q] - m['vals'][m['ph'][q]]) % 24))
    else:
        for q in range(87):
            seq[q] = Z.code_of_slot(Z.TURN_INV[m['vals'][m['ph'][q]]][to_slot(L[q])])
    for q in rng.sample(range(87), noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return seq, plain, m['label']


# ============================================================== state and batches
def blank_state(args):
    return dict(version=VERSION, created=str(Z.now()), elapsed=0.0, sessions=0,
                params=dict(seed=args.seed, real_seeds=args.real_seeds, restarts=args.restarts,
                            iters=args.iters),
                sizes={f: len(family_members(f, DATE_MEMBERS)) for f in FAMILIES},
                targets=json.loads(json.dumps(DEFAULT_TARGETS)),
                fam={f: dict(real=[], shuffle=[], doublet=[], randkey=[], controls=[]) for f in FAMILIES},
                unit_secs={}, log=[])


def stages(state, only):
    out = []
    for f in (only or FAMILIES):
        st, tg = state['fam'][f], state['targets'][f]
        nreal = 1 if f == 'DE' else state['params']['real_seeds']
        if len(st['real']) < nreal:
            out.append((0, -1.0, f, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], f, 'controls'))
        for kind in ('doublet', 'shuffle', 'randkey', 'controls'):
            frac = len(st[kind]) / tg[kind]
            if frac < 1:
                out.append((2, frac, f, kind))
    return sorted(out)


def batch_size(state, f, kind, workers):
    if kind == 'real':
        return max(1, (1 if f == 'DE' else state['params']['real_seeds']) - len(state['fam'][f]['real']))
    n = workers * 4 if f == 'DE' else max(1, math.ceil(2 * workers / len(family_members(f, DATE_MEMBERS))))
    return max(1, min(n, state['targets'][f][kind] - len(state['fam'][f][kind])))


def do_batch(pool, state, SEQ, f, kind, count):
    st, P = state['fam'][f], state['params']
    seed0 = P['seed']
    items, meta = [], []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append((SEQ, DATE_MEMBERS, Z.derive_seed(seed0, 'real', f, n)))
        elif kind in ('shuffle', 'doublet'):
            rng = random.Random(Z.derive_seed(seed0, 'null', f, kind, n))
            if kind == 'shuffle':
                s = list(SEQ)
                rng.shuffle(s)
            else:
                s = M.doublet_shuffle(list(SEQ), rng)
            items.append((s, DATE_MEMBERS, Z.derive_seed(seed0, 'nullsearch', f, kind, n)))
        elif kind == 'randkey':
            items.append((SEQ, random_members(Z.derive_seed(seed0, 'randkey', n)),
                          Z.derive_seed(seed0, 'randsearch', f, n)))
        else:
            noise = 0 if n % 2 == 0 else 2
            seq, plain, lab = plant(f, n, noise, Z.derive_seed(seed0, 'plant', f, n))
            items.append((seq, DATE_MEMBERS, Z.derive_seed(seed0, 'ctrlsearch', f, n)))
            meta.append((noise, plain, lab))
    t0 = time.time()
    res = run_family(pool, f, items, P)
    dt = time.time() - t0
    for k, members_res in enumerate(res):
        sm = summarize(members_res)
        if kind == 'real':
            order = sorted(range(len(members_res)), key=lambda i: -members_res[i][0])[:5]
            sm['top'] = [dict(score=round(members_res[i][0], 5), text=members_res[i][1],
                              label=members_res[i][2]) for i in order]
            st['real'].append(sm)
        elif kind == 'controls':
            noise, plain, lab = meta[k]
            acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
            st['controls'].append(dict(stat=sm['stat'], noise=noise, found=int(sm['best'] == lab),
                                       acc=round(acc, 3), text=sm['text'], plain=plain))
        else:
            st[kind].append(dict(stat=sm['stat'], best=sm['best'], words=sm['words']))
    state['unit_secs'][f] = dt / count
    return dt


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def verdict(p_d, p_r, power):
    if p_d is None:
        return 'waiting for copies'
    if p_d <= BONF and (p_r or 1) <= BONF:
        return 'BEATS copies AND random keys past Bonferroni: read the decrypt by eye'
    if p_d <= BONF:
        return 'beats copies, but not random-number keys (the date is not special)'
    if p_d <= 0.05:
        return 'below 0.05, not past Bonferroni (expected by chance in 3 families)'
    if power is not None and power < 0.5:
        return 'chance level, LOW POWER: inconclusive'
    return 'chance level'


_SEG = []


def segment(text):
    """Words in a decrypt for display only (scores never use this): the big ~11,700-word list
    from the styles test, longer words preferred; the small built-in list if that is missing."""
    if not _SEG:
        try:
            import Dorabella_styles as S
            d = S.StyleDict(S.base_dictionary('english.txt'), 'S0')
            _SEG.append(lambda t: d.segment(S.fold(t)))
        except Exception:
            _SEG.append(lambda t: H.segment(t.translate(H.TOL)))
    return _SEG[0](text)


def build_report(state):
    L = ["DORABELLA DATE-KEY REPORT  ('July 14 97' as the key)",
         f"written {Z.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing "
         f"{state['params']['restarts']} restarts x {state['params']['iters']} steps per setting",
         f"cipher: {state.get('source', '')}",
         f"Bonferroni line for one of 3 families: p < {BONF:.4f}", "",
         "=" * 104,
         "HEADLINE (statistic = best quadgram score over the family's date settings; real = median of seeds;",
         "          p_rand = real cipher with the DATE keys vs the real cipher with RANDOM keys of the same size;",
         "          power = planted controls above the 95th percentile of the pair-preserving copies)",
         "=" * 104,
         f"{'fam':<4} {'settings':>8} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} {'p_shuf':>7} "
         f"{'n_rand':>6} {'p_rand':>7} {'power':>6}  verdict"]
    detail = []
    for f in FAMILIES:
        st = state['fam'][f]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        nd = [x['stat'] for x in st['doublet']]
        ns = [x['stat'] for x in st['shuffle']]
        nr = [x['stat'] for x in st['randkey']]
        p_d = Z.pval(med, nd) if med is not None and nd else None
        p_s = Z.pval(med, ns) if med is not None and ns else None
        p_r = Z.pval(med, nr) if med is not None and nr else None
        q95 = Z.pct(nd, 0.95)
        ctrl = st['controls']
        power = (sum(1 for c in ctrl if c['stat'] > q95) / len(ctrl)) if (ctrl and q95 is not None) else None
        L.append(f"{f:<4} {state['sizes'][f]:>8} {fmt(med):>8} {len(nd):>6} {fmt(p_d, 4):>7} {len(ns):>6} "
                 f"{fmt(p_s, 4):>7} {len(nr):>6} {fmt(p_r, 4):>7} "
                 f"{('-' if power is None else f'{power:.0%}'):>6}  {verdict(p_d, p_r, power)}")
        detail.append((f, st, med, nd, ns, nr, q95))
    L += ["", "A family matters only if the real cipher beats the pair-preserving copies AND the",
          "random-number keys past the Bonferroni line, its controls show it can find a message,",
          "and the decrypt reads as English."]
    for f, st, med, nd, ns, nr, q95 in detail:
        L += ["", "-" * 104, f"{f}: {FAMILY_DESC[f]}   ({state['sizes'][f]} settings)", "-" * 104]
        if not st['real']:
            L.append("real run not done yet")
            continue
        L.append("real statistic per seed: " + ', '.join(f"{r['stat']:.3f}" for r in st['real'])
                 + f"   median {med:.3f}")
        for name, nul in (('pair-preserving', nd), ('plain shuffle', ns), ('random-number keys', nr)):
            if nul:
                L.append(f"  vs {name:<19} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{Z.pct(nul, .95):.3f}, max {max(nul):.3f}, p = {Z.pval(med, nul):.4f}")
        best = max(st['real'], key=lambda r: r['stat'])
        L.append("best real decrypts (seed with the highest statistic):")
        for t in best['top']:
            L.append(f"  [{t['score']:.3f}] {t['label']}")
            for a in range(0, 87, 29):
                L.append("      " + t['text'][a:a + 29])
            L.append("      words: " + segment(t['text']))
        if st['controls']:
            for nz in sorted({c['noise'] for c in st['controls']}):
                cs = [c for c in st['controls'] if c['noise'] == nz]
                above = sum(1 for c in cs if q95 is not None and c['stat'] > q95)
                L.append(f"controls, {nz} misread symbols: n={len(cs)}, median stat "
                         f"{statistics.median(c['stat'] for c in cs):.3f}, true setting found "
                         f"{sum(c['found'] for c in cs)}/{len(cs)}, median letters right "
                         f"{statistics.median(c['acc'] for c in cs):.0%}, above null 95th pct "
                         f"{above}/{len(cs) if q95 is not None else '-'}")
            ex = max(st['controls'], key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:56]}")
            L.append(f"                                   found   {ex['text'][:56]}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== self-test
def selftest(args, pool, P, SEQ, workers):
    print(f"\n{Z.stamp()} SELF-TEST")
    mine = [Z.ELGAR[Z.slot_of(c)] for c in SEQ]
    pub = [H.LIDX[c] for c in H.fold_base(Z.STANDARD_DECRYPT)]
    agree = sum(a == b for a, b in zip(mine, pub))
    print(f"Elgar's key on this transcription agrees with the published decrypt at {agree}/87")
    assert agree >= 75
    for f in FAMILIES:
        print(f"  {f}: {len(family_members(f, DATE_MEMBERS)):>4} settings  ({FAMILY_DESC[f]})")
    rm = random_members(1)
    print(f"  random-key family check: {len(rm['sub'])} + {len(rm['col'])} settings "
          f"(date family {len(DATE_MEMBERS['sub'])} + {len(DATE_MEMBERS['col'])})")
    # plant / decode round trip with Elgar's key and the true setting
    for f in FAMILIES:
        seq, plain, lab = plant(f, 2, 0, Z.derive_seed('st-roundtrip', f))
        if f == 'DE':
            r = [x for x in de_job(seq, family_members(f, DATE_MEMBERS)) if x[2] == lab][0]
            ok = sum(a == b for a, b in zip(r[1], plain)) / 87
            print(f"  round trip DE ({lab}): {ok:.0%} letters right with the true setting")
    print("\nOne planted-English control per family (no misreads), through the full family search:")
    est = {}
    for f in (args.only or FAMILIES):
        seq, plain, lab = plant(f, 3, 0, Z.derive_seed('selftest', f))
        t0 = time.time()
        sm = summarize(run_family(pool, f, [(seq, DATE_MEMBERS, 11)], P)[0])
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
        s = M.doublet_shuffle(list(SEQ), random.Random(Z.derive_seed('selftest-null', f)))
        nm = summarize(run_family(pool, f, [(s, DATE_MEMBERS, 5)], P)[0])
        est[f] = dt
        print(f"  {f}: control {sm['stat']:.3f} (a pair-preserving copy: {nm['stat']:.3f}), true setting "
              f"found: {'yes' if sm['best'] == lab else 'no'}, letters right {acc:.0%}, {H.fmt_secs(dt)} per search")
        print(f"      setting: {lab}")
        print(f"      planted: {plain[:58]}")
        print(f"      decoded: {sm['text'][:58]}")
    tot = 0.0
    print("\nEstimated time to reach the targets:")
    for f in (args.only or FAMILIES):
        tg = DEFAULT_TARGETS[f]
        n = sum(tg.values()) + (1 if f == 'DE' else args.real_seeds)
        per = est[f] if f != 'DE' else est[f] / max(1, workers)
        tot += n * per
        print(f"  {f}: {n} searches x {H.fmt_secs(per)} = {H.fmt_secs(n * per)}")
    print(f"  total about {H.fmt_secs(tot)}")


# ============================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='date_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--real-seeds', type=int, default=5)
    ap.add_argument('--restarts', type=int, default=6)
    ap.add_argument('--iters', type=int, default=10000)
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
        if state.get('sizes') != {f: len(family_members(f, DATE_MEMBERS)) for f in FAMILIES}:
            raise SystemExit("saved state was made with different settings; use --fresh or another --out")
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
                    for kind in ('shuffle', 'doublet', 'randkey'):
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
