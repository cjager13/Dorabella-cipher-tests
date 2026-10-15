#!/usr/bin/env python3
"""
Dorabella CRIB CHECK: is the more-cribs C1 result (p = 0.009) about the CRIBS, or about the
cipher in general? And does a crib plus Elgar's own way of writing do any better?

Needs Dorabella_morecribs.py, Dorabella_cribs.py, Dorabella_elgarstyle.py and everything they
need, in the same folder, plus the FINISHED morecribs_out\\state.json (it is read, never
changed). Tests E1/E2 below also need your private Elgar_letters_to_1897.txt next to the
script; it is only read locally.

  python Dorabella_cribcheck.py --selftest         # checks the copies line up, times each part
  python Dorabella_cribcheck.py --stop-when-done   # results in cribcheck_out\\
  python Dorabella_cribcheck.py --report-only

Ctrl+C any time; the same command resumes.

THE IDEA
C1 fixes a crib and searches the rest of the key. Its score can be high for two reasons: a
crib really sits in the message, or the cipher simply completes well under ANY search, crib
or not. So the same search is run with NO crib, on the real cipher and on the very same
1,000 pair-preserving copies morecribs used (rebuilt from the same seeds and checked). The
gap "with crib minus without crib" is what the crib itself is worth. If a crib really
belongs, fixing it costs the real cipher almost nothing while costing the copies more.

The second half repeats the crib test with the quadgram model of Elgar's own voice (style
E3 of the elgarstyle test: Dracula rewritten with no AND and his abbreviations, plus his own
letters before 14 July 1897 at 4 times their weight), crib and no crib, the same way.

PRE-REGISTERED (frozen 14 Oct 2026, before any run)
  D1  CRIB WORTH, DRACULA ENGLISH. Statistic: median over 5 seeds of the morecribs C1 score
      minus median over 5 seeds of the no-crib score, real cipher; for each of the 1,000
      pair-preserving copies, its morecribs C1 score minus its no-crib score.
  DE  CRIB WORTH, ELGAR'S VOICE. The same with the E3 model: the morecribs crib search
      (same 238 eligible cribs, same settings) minus the no-crib search, both under E3,
      on the real cipher and the same 1,000 pair-preserving copies.
  Two tests, so the line is p < 0.05/2 = 0.025, against pair-preserving copies.
  A pass also needs a completed text that READS as English beyond the crib; without that
  it is written up as a statistical oddity, not a solution.
  Reported, not tested: N1 (no-crib score, Dracula) and NE (no-crib, E3) against their
  copies; CE (crib score under E3) against its copies; plain-shuffle versions of each;
  which cribs win in the copies; how often WELLS fits the real cipher against the copies
  (WELLS was found by looking, so that line is exploratory).
  Search: the no-crib search is the C1 annealer with nothing fixed, 30 restarts x 50,000
  steps in each of the 5 reading orders, the 9 uncertain symbols chosen inside the search.
  (Deeper than C1's deep pass, because with nothing fixed the whole key has to be found; a
  trial on the real cipher showed 10 x 20,000 falling well short of the best score.)
  Copies: all 1,000 pair-preserving copies for N1/NE/CE; 200 plain shuffles each, for
  reference only. Power: planted English with a crib (morecribs' own
  controls, rebuilt from the same seeds), counted as found when its gap beats the copies'
  95th percentile.
"""
import json
import math
import os
import statistics
import sys
import time

import Dorabella_morecribs as MC                     # noqa: E402  (applies the crib list)
import Dorabella_elgarstyle as ES                    # noqa: E402,F401  (adds the E3 model)

C = MC.C
S = C.S
H = C.H
mp = C.mp

VERSION = 1
MORE_STATE = os.path.join('morecribs_out', 'state.json')
BONF = 0.05 / 2


# ---------------------------------------------------------------- search jobs (run in workers)
def free_job(corpus, style, seq, alts, restarts, iters, seed):
    """No crib: the C1 annealer with nothing fixed, every reading order. Best (score, text, order)."""
    S.use_style(corpus, style)
    best = (-1e18, '', '')
    for k, o in enumerate(C.ORDERS):
        R = C.ORDER_PERM[o]
        c = [S.slot_of(seq[r]) for r in R]
        A = [S.slot_of(alts[r]) if r in alts else -1 for r in R]
        sc, txt = C.anneal_fixed(c, A, {}, set(), restarts, iters, S.derive_seed(seed, k))
        if sc > best[0]:
            best = (sc, txt, o)
    return best


def crib_job(corpus, style, seq, alts, items, restarts, iters, seed):
    """C.complete_job with the scoring model as a parameter."""
    S.use_style(corpus, style)
    out = []
    for k, (w, o, i, slots) in enumerate(items):
        R = C.ORDER_PERM[o]
        c = [S.slot_of(seq[r]) for r in R]
        A = [S.slot_of(alts[r]) if r in alts else -1 for r in R]
        fixed, frozen = {}, set()
        for k2, s in enumerate(slots):
            q = i + k2
            if c[q] != s:
                A[q], c[q] = c[q], s
            frozen.add(q)
            fixed[s] = C.LIDX[w[k2]]
        out.append(C.anneal_fixed(c, A, fixed, frozen, restarts, iters, S.derive_seed(seed, k)))
    return out


def crib_search(pool, corpus, style, seq, alts, cribs, P, seed):
    """C.c1_search with the scoring model as a parameter (same two passes, same settings)."""
    plist = [(w, o, i, sl) for w in cribs for (o, i, sl) in C.placements(w, seq, alts)]
    if not plist:
        return dict(stat=-99.0, crib='', text='', order='', n=0)
    chunk = max(1, math.ceil(len(plist) / (4 * P['workers'])))
    chunks = [plist[k:k + chunk] for k in range(0, len(plist), chunk)]
    r1 = C.run(pool, crib_job, [(corpus, style, seq, alts, ch, P['r1'], P['i1'], S.derive_seed(seed, 'a', k))
                                for k, ch in enumerate(chunks)])
    scored = sorted(((res[0], it) for ch, rr in zip(chunks, r1) for it, res in zip(ch, rr)),
                    key=lambda x: -x[0])
    deep = [it for _, it in scored[:P['keep']]]
    r2 = C.run(pool, crib_job, [(corpus, style, seq, alts, [it], P['r2'], P['i2'], S.derive_seed(seed, 'b', k))
                                for k, it in enumerate(deep)])
    res = sorted(((r[0][0], r[0][1], it) for r, it in zip(r2, deep)), key=lambda x: -x[0])
    b = res[0]
    return dict(stat=b[0], crib=b[2][0], text=b[1], order=b[2][1], n=len(plist),
                top=[dict(score=x[0], crib=x[2][0], order=x[2][1], text=x[1]) for x in res[:3]])


# ---------------------------------------------------------------- the copies and controls morecribs used
def copy_of(P, kind, m):
    return S.null_tokens(kind, S.SEQ, S.ALTS, S.derive_seed(P['seed'], 'null', kind, m))


def control_of(elig, P, m):
    noise = 0 if m % 2 == 0 else 2
    s, a, plain, crib, order = C.plant_c1(elig, m, noise, S.derive_seed(P['seed'], 'plant', m))
    return s, a, plain, crib, order


# ---------------------------------------------------------------- state
def load_more():
    if not os.path.exists(MORE_STATE):
        raise SystemExit(f"{MORE_STATE} not found: run this next to the finished morecribs_out folder")
    with open(MORE_STATE, encoding='utf-8') as fh:
        ms = json.load(fh)
    if ms['seq'] != S.SEQ:
        raise SystemExit("morecribs_out was made with a different transcription")
    return ms


def blank_state(args, ms):
    P = dict(ms['params'])
    P.update(workers=args.workers, free_r=args.free_restarts, free_i=args.free_steps, real_seeds=5)
    n_pair = len(ms['st']['C1']['doublet'])
    n_shuf = len(ms['st']['C1']['shuffle'])
    n_ctrl = len(ms['st']['C1']['controls'])
    tg = {'N1': dict(real=5, doublet=n_pair, shuffle=min(n_shuf, args.shuffles), controls=n_ctrl),
          'NE': dict(real=5, doublet=n_pair, shuffle=args.shuffles, controls=args.controls_e3),
          'CE': dict(real=5, doublet=n_pair, shuffle=args.shuffles, controls=args.controls_e3)}
    return dict(version=VERSION, created=str(C.now()), elapsed=0.0, sessions=0, seq=S.SEQ,
                eligible=ms['eligible'], params=P, targets=tg,
                st={t: {k: [] for k in ('real', 'doublet', 'shuffle', 'controls')} for t in tg}, log=[])


ORDER = [('N1', 'real'), ('NE', 'real'), ('CE', 'real'), ('N1', 'doublet'), ('N1', 'controls'),
         ('NE', 'doublet'), ('CE', 'doublet'), ('N1', 'shuffle'), ('NE', 'controls'), ('CE', 'controls'),
         ('NE', 'shuffle'), ('CE', 'shuffle')]


def stages(state):
    out = []
    for rank, (t, kind) in enumerate(ORDER):
        if len(state['st'][t][kind]) < state['targets'][t][kind]:
            out.append((rank, t, kind))
    return out


def do_unit(pool, state, corpus, test, kind):
    P = state['params']
    st = state['st'][test]
    elig = state['eligible']
    style = 'S0' if test == 'N1' else 'E3'
    n = len(st[kind])
    if test == 'CE':                    # one crib search at a time, parallel inside
        if kind == 'real':
            s, a, extra = S.SEQ, S.ALTS, {}
        elif kind == 'controls':
            s, a, plain, crib, order = control_of(elig, P, n)
            extra = dict(planted=crib, plain=plain, porder=order)
        else:
            s, a = copy_of(P, kind, n)
            extra = {}
        r = crib_search(pool, corpus, 'E3', s, a, elig, P, S.derive_seed(P['seed'], 'CE', kind, n))
        rec = dict(i=n, stat=round(r['stat'], 5), crib=r['crib'], n=r['n'])
        if kind == 'real':
            rec.update(text=r['text'], order=r['order'], top=r.get('top', []))
        if kind == 'controls':
            rec.update(found=r['crib'], planted=extra['planted'],
                       acc=round(sum(x == y for x, y in zip(r['text'], extra['plain'])) / 87
                                 if r['order'] == extra['porder'] else 0.0, 3))
        st[kind].append(rec)
        return
    # N1 / NE: a batch of no-crib searches in parallel
    count = 1 if kind == 'real' else max(1, min(4 * P['workers'], state['targets'][test][kind] - n))
    jobs = []
    for k in range(count):
        m = n + k
        if kind == 'real':
            s, a = S.SEQ, S.ALTS
        elif kind == 'controls':
            s, a = control_of(elig, P, m)[:2]
        else:
            s, a = copy_of(P, kind, m)
        jobs.append((corpus, style, s, a, P['free_r'], P['free_i'], S.derive_seed(P['seed'], test, kind, m)))
    res = C.run(pool, free_job, jobs)
    for k, (sc, txt, o) in enumerate(res):
        rec = dict(i=n + k, stat=round(sc, 5))
        if kind == 'real':
            rec.update(text=txt, order=o)
        st[kind].append(rec)


# ---------------------------------------------------------------- report
def med(xs):
    return statistics.median(xs) if xs else None


def gaps(crib_list, free_list):
    """Pair up by copy index: crib score minus no-crib score."""
    fr = {x['i']: x['stat'] for x in free_list}
    return [c['stat'] - fr[k] for k, c in enumerate(crib_list) if k in fr]


def line(name, real, nulls):
    if real is None or not nulls:
        return f"{name:<44} (not enough runs yet)"
    return (f"{name:<44} real {real:8.3f}   copies n={len(nulls):>4} mean {statistics.mean(nulls):8.3f} "
            f"95th {C.pct(nulls, .95):8.3f} max {max(nulls):8.3f}   p = {C.pval(real, nulls):.4f}")


def power(ctrl_vals, nulls):
    if not ctrl_vals or not nulls:
        return '-'
    q = C.pct(nulls, .95)
    return f"{sum(1 for v in ctrl_vals if v > q) / len(ctrl_vals):.0%} of {len(ctrl_vals)}"


_WELLS = {}


def wells_fits(P, n):
    if n not in _WELLS:
        _WELLS[n] = (len(C.placements('WELLS', S.SEQ, S.ALTS)),
                     [len(C.placements('WELLS', *copy_of(P, 'doublet', m))) for m in range(n)])
    return _WELLS[n]


def build_report(state, ms, sd):
    st, mst = state['st'], ms['st']['C1']
    c1_real = med([r['stat'] for r in mst['real']])
    n1_real = med([r['stat'] for r in st['N1']['real']])
    ne_real = med([r['stat'] for r in st['NE']['real']])
    ce_real = med([r['stat'] for r in st['CE']['real']])
    d1_real = None if n1_real is None else c1_real - n1_real
    de_real = None if (ne_real is None or ce_real is None) else ce_real - ne_real
    d1_pair = gaps(mst['doublet'], st['N1']['doublet'])
    d1_shuf = gaps(mst['shuffle'], st['N1']['shuffle'])
    de_pair = gaps(st['CE']['doublet'], st['NE']['doublet'])
    de_shuf = gaps(st['CE']['shuffle'], st['NE']['shuffle'])
    d1_ctrl = gaps(mst['controls'], st['N1']['controls'])
    de_ctrl = gaps(st['CE']['controls'], st['NE']['controls'])
    L = ["DORABELLA CRIB CHECK REPORT  (is the more-cribs C1 result about the cribs?  and cribs + Elgar's voice)",
         f"written {C.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}   cribs {len(state['eligible'])} (the morecribs list)", "",
         "=" * 112, "HEADLINE: pre-registered gaps (crib score minus no-crib score), against PAIR-PRESERVING copies",
         "=" * 112,
         line("D1 crib worth, Dracula English", d1_real, d1_pair) + f"   power {power(d1_ctrl, d1_pair)}",
         line("DE crib worth, Elgar's voice (E3)", de_real, de_pair) + f"   power {power(de_ctrl, de_pair)}",
         f"Line to beat: p < {BONF:.3f} for D1 or DE (two tests). A pass also needs a completed text that",
         "READS as English beyond the crib; otherwise it is a statistical oddity, not a solution.", "",
         "-" * 112, "REPORTED, NOT TESTED", "-" * 112,
         line("C1 crib score, Dracula (from morecribs)", c1_real, [x['stat'] for x in mst['doublet']]),
         line("N1 no-crib score, Dracula", n1_real, [x['stat'] for x in st['N1']['doublet']]),
         line("CE crib score, Elgar's voice", ce_real, [x['stat'] for x in st['CE']['doublet']]),
         line("NE no-crib score, Elgar's voice", ne_real, [x['stat'] for x in st['NE']['doublet']]),
         "plain-shuffle copies:",
         line("  D1 vs plain shuffles", d1_real, d1_shuf),
         line("  DE vs plain shuffles", de_real, de_shuf),
         line("  N1 vs plain shuffles", n1_real, [x['stat'] for x in st['N1']['shuffle']]),
         line("  CE vs plain shuffles", ce_real, [x['stat'] for x in st['CE']['shuffle']])]
    # which cribs win in copies; WELLS fits
    from collections import Counter
    for name, lst in (("morecribs C1", mst['doublet']), ("CE", st['CE']['doublet'])):
        if lst:
            top = Counter(x['crib'] for x in lst).most_common(6)
            L.append(f"best crib in pair-preserving copies, {name}: " +
                     ', '.join(f"{w} {k}" for w, k in top))
    real_fits, cf = wells_fits(state['params'], len(mst['doublet']))
    if cf:
        L.append(f"WELLS (found by looking, exploratory): fits the real cipher {real_fits} times; pair-preserving "
                 f"copies mean {statistics.mean(cf):.1f}, p (as many or more) = {C.pval(real_fits, cf):.4f}")
    for t, title in (('N1', "N1 best no-crib completions of the real cipher, Dracula English"),
                     ('NE', "NE best no-crib completions of the real cipher, Elgar's voice"),
                     ('CE', "CE best crib completions of the real cipher, Elgar's voice")):
        if not st[t]['real']:
            continue
        L += ["", "-" * 112, title, "-" * 112]
        L.append("per seed: " + ', '.join(f"{r['stat']:.3f}" + (f" ({r['crib']})" if r.get('crib') else '')
                                          for r in st[t]['real']))
        b = max(st[t]['real'], key=lambda r: r['stat'])
        shown = b.get('top') or [dict(score=b['stat'], crib='', order=b['order'], text=b['text'])]
        for x in shown:
            L.append(f"  [{x['score']:.3f}] {('crib ' + x['crib'] + ' in ') if x['crib'] else ''}{x['order']}")
            for a in range(0, 87, 29):
                L.append("      " + x['text'][a:a + 29])
            L.append("      words: " + sd.segment(x['text']))
    c = st['CE']['controls']
    if c:
        L += ["", f"CE controls: n={len(c)}, planted crib found {sum(1 for x in c if x['found'] == x['planted'])}/"
                  f"{len(c)}, median letters right {statistics.median(x['acc'] for x in c):.0%}"]
    L += ["", "copies are morecribs' own, rebuilt from its seeds; checked by placement counts at the start of every "
              "session", "", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ---------------------------------------------------------------- main
def check_alignment(ms, P, k=5):
    """The rebuilt copies must be the ones morecribs used: same number of crib placements."""
    elig = ms['eligible']
    for kind in ('doublet', 'shuffle'):
        for m in range(min(k, len(ms['st']['C1'][kind]))):
            s, a = copy_of(P, kind, m)
            n = sum(len(C.placements(w, s, a)) for w in elig)
            if n != ms['st']['C1'][kind][m]['n']:
                raise SystemExit(f"copy {kind} {m} does not match morecribs ({n} vs "
                                 f"{ms['st']['C1'][kind][m]['n']} placements); stopping")
    return True


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='cribcheck_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--shuffles', type=int, default=200)
    ap.add_argument('--free-restarts', type=int, default=30)
    ap.add_argument('--free-steps', type=int, default=50000)
    ap.add_argument('--controls-e3', type=int, default=30)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()
    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    if not os.path.exists(ES.LETTERS_FILE):
        raise SystemExit(f"{ES.LETTERS_FILE} is needed next to this script for the Elgar's-voice tests")
    S.SEQ, S.ALTS = S.load_cipher()
    ms = load_more()
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))
    sd = S.StyleDict(S.base_dictionary(args.corpus), 'S0')

    def log(msg):
        line_ = f"{C.stamp()} {msg}"
        print(line_, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line_ + '\n')

    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
    else:
        state = blank_state(args, ms)
    if args.report_only:
        S.write_text(rp, build_report(state, ms, sd))
        return
    state['params']['workers'] = args.workers
    print("checking the rebuilt copies match morecribs...", flush=True)
    check_alignment(ms, state['params'])
    print("  they match.", flush=True)
    S.keep_awake()
    pool = mp.Pool(args.workers, initializer=C.init_worker, initargs=(args.corpus,)) if args.workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, state, ms)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {args.workers}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stg = stages(state)
            if not stg:
                log("all targets met; stopping. Progress and report are saved.")
                break
            _, test, kind = stg[0]
            t0 = time.time()
            do_unit(pool, state, args.corpus, test, kind)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            msg = f"{test} {kind} in {H.fmt_secs(t1 - t0)}  ({len(state['st'][test][kind])}/" \
                  f"{state['targets'][test][kind]})"
            if kind == 'real':
                msg += f"  real {state['st'][test]['real'][-1]['stat']:.3f}"
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{C.stamp()} {msg}", flush=True)
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, ms, sd))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished unit is saved; run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, ms, sd))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


def selftest(args, pool, state, ms):
    P = state['params']
    print(f"\n{C.stamp()} SELF-TEST  ({len(state['eligible'])} cribs; targets {state['targets']})")
    s, a = copy_of(P, 'doublet', 0)
    t0 = time.time()
    (f0,) = C.run(pool, free_job, [(args.corpus, 'S0', s, a, P['free_r'], P['free_i'], 11)])
    dt_free = time.time() - t0
    gap = ms['st']['C1']['doublet'][0]['stat'] - f0[0]
    print(f"no-crib search, copy 0, Dracula: {f0[0]:.3f} ({H.fmt_secs(dt_free)} for one; batches run "
          f"{P['workers']} at a time); its morecribs C1 score {ms['st']['C1']['doublet'][0]['stat']:.3f}, "
          f"gap {gap:+.3f}")
    (f1,) = C.run(pool, free_job, [(args.corpus, 'S0', S.SEQ, S.ALTS, P['free_r'], P['free_i'], 12)])
    print(f"no-crib search, REAL cipher, Dracula: {f1[0]:.3f} in {f1[2]}: {f1[1]}")
    (f2,) = C.run(pool, free_job, [(args.corpus, 'E3', S.SEQ, S.ALTS, P['free_r'], P['free_i'], 13)])
    print(f"no-crib search, REAL cipher, Elgar's voice: {f2[0]:.3f} in {f2[2]}: {f2[1]}")
    cs, ca, plain, crib, order = control_of(state['eligible'], P, 0)
    t0 = time.time()
    r = crib_search(pool, args.corpus, 'E3', cs, ca, state['eligible'], P, 14)
    dt_crib = time.time() - t0
    acc = sum(x == y for x, y in zip(r['text'], plain)) / 87 if r['order'] == order else 0.0
    print(f"crib search, Elgar's voice, planted control: planted {crib}, found {r['crib']}, "
          f"{acc:.0%} letters right, {H.fmt_secs(dt_crib)}")
    print(f"  planted {plain[:60]}\n  found   {r['text'][:60]}")
    tg = state['targets']
    free_units = sum(tg[t][k] for t in ('N1', 'NE') for k in ('doublet', 'shuffle', 'controls'))
    crib_units = sum(tg['CE'][k] for k in ('real', 'doublet', 'shuffle', 'controls'))
    est = free_units * dt_free / max(1, P['workers']) + crib_units * dt_crib
    print(f"\nEstimated time to finish: about {H.fmt_secs(est)}")


if __name__ == '__main__':
    mp.freeze_support()
    if not any(a.startswith('--out') for a in sys.argv):
        sys.argv += ['--out', 'cribcheck_out']
    main()
