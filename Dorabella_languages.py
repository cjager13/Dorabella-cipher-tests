#!/usr/bin/env python3
"""
Dorabella LANGUAGES test: is the message in French, German, Latin or Italian instead of English?

Needs Dorabella_styles.py (and what it needs: Dorabella_hillclimb.py, Dorabella_marathon.py,
Dorabella_tests.py), Dorabella_transcription.csv, and one plain-text book per language in the
same folder:  english.txt  french.txt  german.txt  latin.txt  italian.txt
(any language file that is missing is simply left out).

  python Dorabella_languages.py --selftest        # checks each text, one control per language, timing
  python Dorabella_languages.py                   # start, or resume where it left off
  python Dorabella_languages.py --stop-when-done  # exit when every target is met
  python Dorabella_languages.py --report-only / --fresh / --hours N / --only FR LA

Ctrl+C any time; the same command resumes. Progress: languages_out/state.json, report.txt.

PRE-REGISTERED (frozen before any real run)
  FR French   DE German   LA Latin   IT Italian     (4 languages: Bonferroni line p < 0.0125)
  EN English  run with identical settings as the reference (already tested: expected at chance)

FOR EACH LANGUAGE
  * the book is cleaned: Project Gutenberg header and licence removed, accents removed
    (German umlauts written AE OE UE, as a German writes them without the dots), letters
    folded to Elgar's 24-letter alphabet exactly as for English (J -> I, V -> U)
  * a quadgram model is trained on the first 85% of the book; the last 15% is held back
    and used only for the planted controls, so a control is never text the model has seen
  * the search is ANY one-symbol-one-letter key, by simulated annealing, over the 5 reading
    orders, with the 9 uncertain symbols chosen inside the search (as in Dorabella_styles.py)

METHOD (unchanged): statistic = best score over the 5 orders; real = median over several
search seeds. Each plain shuffle and pair-preserving shuffle goes through the IDENTICAL search
in the SAME language model. Scores are only comparable within a language (each language's
model has its own scale), never across languages. Planted text from the held-back part of
each book (0 or 2 misread symbols, uncertain symbols as in the real cipher) measures power.
p = (1+k)/(n+1).
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import random
import re
import signal
import statistics
import time
import unicodedata
from collections import Counter

import Dorabella_styles as S

H, M = S.H, S.M
VERSION = 1
LANGS = [('EN', 'english.txt', 'English (reference)'), ('FR', 'french.txt', 'French'),
         ('DE', 'german.txt', 'German'), ('LA', 'latin.txt', 'Latin'), ('IT', 'italian.txt', 'Italian')]
CONFIRM = ('LC', 'latin_church.txt', 'Latin, Church (confirmation of LA)')
NAME = {c: n for c, _, n in LANGS + [CONFIRM]}
FILE = {c: f for c, f, _ in LANGS + [CONFIRM]}
NEW = ['FR', 'DE', 'LA', 'IT']
BONF = 0.05 / len(NEW)
HOLDOUT = 0.15
ORDERS = S.ORDERS
SEQ = ALTS = None
TABLES = {}        # per process: language -> quadgram table


# ============================================================== text preparation
def strip_gutenberg(raw):
    """Keep only the book text. Works when several Gutenberg books are pasted into one file:
    every '*** START OF ... ***' to the next '*** END OF' is kept; headers and the English
    licences between books are dropped."""
    starts = [m.end() for m in re.finditer(r'\*\*\* ?START OF[^\n]*\n', raw)]
    if not starts:
        return raw
    parts = []
    for s in starts:
        m = re.compile(r'\*\*\* ?END OF|End of (the )?Project Gutenberg').search(raw, s)
        parts.append(raw[s:m.start() if m else len(raw)])
    return '\n'.join(parts)


def normalise(raw, lang):
    t = strip_gutenberg(raw)
    if lang == 'DE':
        for a, b in (('ä', 'ae'), ('ö', 'oe'), ('ü', 'ue'), ('Ä', 'Ae'), ('Ö', 'Oe'), ('Ü', 'Ue')):
            t = t.replace(a, b)
    for a, b in (('ß', 'ss'), ('æ', 'ae'), ('Æ', 'AE'), ('œ', 'oe'), ('Œ', 'OE')):
        t = t.replace(a, b)
    t = unicodedata.normalize('NFKD', t)
    t = ''.join(ch for ch in t if not unicodedata.combining(ch))
    words = re.findall(r'[A-Z]+', t.upper())
    return [H.fold_base(w) for w in words]


def load_language(lang):
    with open(FILE[lang], encoding='utf-8', errors='ignore') as fh:
        words = normalise(fh.read(), lang)
    letters = ''.join(words)
    cut = int(len(letters) * (1 - HOLDOUT))
    return words, letters[:cut], letters[cut:]


def build_table(train):
    NL = 24
    li = [H.LIDX[c] for c in train]
    cnt = Counter()
    for i in range(len(li) - 3):
        cnt[((li[i] * NL + li[i + 1]) * NL + li[i + 2]) * NL + li[i + 3]] += 1
    tot = sum(cnt.values())
    floor = math.log(0.01 / tot)
    q = [floor] * (NL ** 4)
    for k, v in cnt.items():
        q[k] = math.log(v / tot)
    return q


def use_language(lang):
    if lang not in TABLES:
        TABLES[lang] = build_table(load_language(lang)[1])
    S.QL = TABLES[lang]


class WordList:
    """Words of a language for DISPLAY only (scores never use this)."""

    def __init__(self, words):
        c = Counter(words)
        self.words = {w for w, n in c.items() if n >= 3 and 3 <= len(w) <= 16}
        self.prefixes = {w[:k] for w in self.words for k in range(1, len(w) + 1)}
        self.maxlen = max((len(w) for w in self.words), default=0)

    def segment(self, t):
        n = len(t)
        best, choice = [0] * (n + 1), [0] * (n + 1)
        for i in range(n - 1, -1, -1):
            b, ch = best[i + 1], 0
            for L in range(1, min(self.maxlen, n - i) + 1):
                sub = t[i:i + L]
                if sub not in self.prefixes:
                    break
                if L >= 3 and sub in self.words and (L - 2) ** 2 + best[i + L] > b:
                    b, ch = (L - 2) ** 2 + best[i + L], L
            best[i], choice[i] = b, ch
        out, gap, i = [], [], 0
        while i < n:
            if choice[i]:
                if gap:
                    out.append(''.join(gap).lower())
                    gap = []
                out.append(t[i:i + choice[i]])
                i += choice[i]
            else:
                gap.append(t[i])
                i += 1
        if gap:
            out.append(''.join(gap).lower())
        return ' '.join(out)


# ============================================================== searching
def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def order_job(lang, seq, alts, order, restarts, iters, seed):
    use_language(lang)
    R = S.ORDER_PERM[order]
    c = [S.slot_of(seq[r]) for r in R]
    A = [S.slot_of(alts[r]) if r in alts else -1 for r in R]
    s, text, flips = S.anneal(c, A, restarts, iters, seed)
    return s, text, sorted(R[q] + 1 for q in flips)


def run_jobs(pool, jobs):
    if pool is None:
        return [order_job(*a) for a in jobs]
    res = pool.starmap_async(order_job, jobs, chunksize=1)
    while True:
        try:
            return res.get(60)
        except mp.TimeoutError:
            continue


def key_accuracy(text, found_order, plain, true_order):
    """Share of SYMBOLS decoded to the right letter, judged by where each symbol sits on the
    page, so a decode that only puts the three lines in a different order still counts."""
    Rf, Rt = S.ORDER_PERM[found_order], S.ORDER_PERM[true_order]
    true_at = {Rt[i]: plain[i] for i in range(87)}
    return sum(text[i] == true_at[Rf[i]] for i in range(87)) / 87


def plant(heldout, k, noise, seed):
    """87 letters of held-back text in the language, a random key and reading order,
    misread symbols, and alternates at the uncertain positions as in the real cipher."""
    rng = random.Random(seed)
    start = rng.randrange(0, max(1, len(heldout) - 87))
    plain = heldout[start:start + 87]
    slots = list(range(24))
    rng.shuffle(slots)
    order = rng.choice(ORDERS)
    R = S.ORDER_PERM[order]
    seq = [0] * 87
    for i, ch in enumerate(plain):
        seq[R[i]] = S.code_of_slot(slots[H.LIDX[ch]])
    for q in rng.sample([i for i in range(87) if i not in ALTS], noise):
        seq[q] = S.neighbour(seq[q], rng)
    alts = {}
    for q in sorted(ALTS):
        wrong = S.neighbour(seq[q], rng)
        if rng.random() < 0.5:
            alts[q] = wrong
        else:
            alts[q], seq[q] = seq[q], wrong
    return seq, alts, plain, order


# ============================================================== state and batches
def blank_state(args, langs):
    return dict(version=VERSION, created=str(S.now()), elapsed=0.0, sessions=0, seq=SEQ,
                alts={str(k): v for k, v in ALTS.items()},
                params=dict(seed=args.seed, real_seeds=args.real_seeds, restarts=args.restarts,
                            iters=args.iters, holdout=HOLDOUT),
                targets={L: dict(shuffle=args.copies, doublet=args.copies, controls=args.controls) for L in langs},
                st={L: dict(real=[], shuffle=[], doublet=[], controls=[]) for L in langs},
                files={}, log=[])


def stage_list(state, langs):
    out = []
    for L in langs:
        st, tg = state['st'][L], state['targets'][L]
        if len(st['real']) < state['params']['real_seeds']:
            out.append((0, 0.0, L, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], L, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            frac = len(st[kind]) / tg[kind]
            if frac < 1:
                out.append((2, frac, L, kind))
    return sorted(out)


def batch_size(state, L, kind, workers):
    if kind == 'real':
        return state['params']['real_seeds'] - len(state['st'][L]['real'])
    n = max(1, math.ceil(2 * workers / len(ORDERS)))
    return max(1, min(n, state['targets'][L][kind] - len(state['st'][L][kind])))


def do_batch(pool, state, held, L, kind, count):
    P, st = state['params'], state['st'][L]
    items = []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append((SEQ, ALTS, S.derive_seed(P['seed'], 'real', L, n), None))
        elif kind in ('shuffle', 'doublet'):
            s, a = S.null_tokens(kind, SEQ, ALTS, S.derive_seed(P['seed'], 'null', L, kind, n))
            items.append((s, a, S.derive_seed(P['seed'], 'nullsearch', L, kind, n), None))
        else:
            noise = 0 if n % 2 == 0 else 2
            s, a, plain, order = plant(held[L], n, noise, S.derive_seed(P['seed'], 'plant', L, n))
            items.append((s, a, S.derive_seed(P['seed'], 'ctrlsearch', L, n), (noise, plain, order)))
    jobs = [(L, s, a, o, P['restarts'], P['iters'], S.derive_seed(sd, o)) for s, a, sd, _ in items for o in ORDERS]
    t0 = time.time()
    res = run_jobs(pool, jobs)
    for k, (s, a, sd, meta) in enumerate(items):
        per = res[k * len(ORDERS):(k + 1) * len(ORDERS)]
        bi = max(range(len(ORDERS)), key=lambda i: per[i][0])
        stat, text, used = per[bi]
        if kind == 'real':
            ranked = sorted(range(len(ORDERS)), key=lambda i: -per[i][0])
            st['real'].append(dict(stat=stat, top=[dict(order=ORDERS[i], score=per[i][0], text=per[i][1],
                                                         alts=per[i][2]) for i in ranked]))
        elif kind in ('shuffle', 'doublet'):
            st[kind].append(dict(stat=round(stat, 5)))
        else:
            noise, plain, order = meta
            acc = key_accuracy(text, ORDERS[bi], plain, order)
            st['controls'].append(dict(stat=round(stat, 5), noise=noise, acc=round(acc, 3),
                                       right_order=int(ORDERS[bi] == order), text=text, plain=plain))
    return time.time() - t0


# ============================================================== report
def fmt(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def verdict(L, p_d, p_s, power):
    if p_d is None:
        return 'waiting for copies'
    line = BONF if L in NEW else 0.05
    if p_d <= line and (p_s or 1) <= line:
        return 'BEATS both nulls' + (' past Bonferroni' if L in NEW else '') + ': read the decrypts by eye'
    if p_d <= 0.05:
        return 'below 0.05, not past Bonferroni (expected by chance among 4 languages)'
    if power is not None and power < 0.5:
        return 'chance level, LOW POWER: inconclusive'
    return 'chance level'


def build_report(state, langs, wl):
    L_ = ["DORABELLA LANGUAGES REPORT",
          f"written {S.now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
          f"sessions {state['sessions']}",
          f"settings: real seeds {state['params']['real_seeds']} (median used), annealing "
          f"{state['params']['restarts']} restarts x {state['params']['iters']} steps per reading order, "
          f"5 reading orders, {len(ALTS)} uncertain symbols solved in the search",
          "texts: " + '; '.join(f"{L} {FILE[L]} ({state['files'].get(L, {}).get('letters', 0):,} letters)"
                                for L in langs),
          f"Bonferroni line for one of the 4 new languages: p < {BONF:.4f}", "",
          "=" * 104,
          "HEADLINE (statistic = best quadgram score in that language's own model; real = median over seeds;",
          "          power = planted held-back text above the 95th percentile of the pair-preserving copies)",
          "Scores are on each language's own scale: compare a language only with its own copies.",
          "=" * 104,
          f"{'language':<22} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} {'p_shuf':>7} "
          f"{'power':>6} {'ctrl med':>8}  verdict"]
    detail = []
    for L in langs:
        st = state['st'][L]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        ctrl = [c['stat'] for c in st['controls']]
        q95 = S.pct(dn, 0.95)
        power = sum(1 for c in ctrl if c > q95) / len(ctrl) if (ctrl and q95 is not None) else None
        p_d, p_s = S.pval(med, dn), S.pval(med, sn)
        L_.append(f"{L + ' ' + NAME[L]:<22} {fmt(med):>8} {len(dn):>6} {fmt(p_d, 4):>7} {len(sn):>6} "
                  f"{fmt(p_s, 4):>7} {('-' if power is None else f'{power:.0%}'):>6} "
                  f"{fmt(statistics.median(ctrl) if ctrl else None, 2):>8}  {verdict(L, p_d, p_s, power)}")
        detail.append((L, st, med, dn, sn, q95))
    L_ += ["", "A language matters only if the real cipher beats the PAIR-PRESERVING copies past the",
           "Bonferroni line, its controls show it can find a message, and the decrypt reads as that language."]
    for L, st, med, dn, sn, q95 in detail:
        L_ += ["", "-" * 104, f"{L}: {NAME[L]}  ({FILE[L]})", "-" * 104]
        if not st['real']:
            L_.append("real run not done yet")
            continue
        L_.append("real statistic per seed: " + ', '.join(f"{r['stat']:.3f}" for r in st['real'])
                  + f"   median {med:.3f}")
        for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L_.append(f"  vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                          f"{S.pct(nul, .95):.3f}, max {max(nul):.3f}, p = {S.pval(med, nul):.4f}")
        b = max(st['real'], key=lambda r: r['stat'])
        L_.append(f"best real decrypts (from the seed with the highest statistic, {b['stat']:.3f}):")
        for t in b['top'][:3]:
            alt = f"  alternate readings used at {t['alts']}" if t['alts'] else ''
            L_.append(f"  [{t['score']:.3f}] {t['order']}{alt}")
            for a in range(0, 87, 29):
                L_.append("      " + t['text'][a:a + 29])
            if wl.get(L):
                L_.append("      words: " + wl[L].segment(t['text']))
        cs = st['controls']
        if cs:
            for nz in (0, 2):
                c2 = [c for c in cs if c['noise'] == nz]
                if not c2:
                    continue
                above = sum(1 for c in c2 if q95 is not None and c['stat'] > q95)
                L_.append(f"controls, {nz} misread symbols: n={len(c2)}, median stat "
                          f"{statistics.median(c['stat'] for c in c2):.3f}, right reading order "
                          f"{sum(c['right_order'] for c in c2)}/{len(c2)}, median symbols decoded right "
                          f"{statistics.median(c['acc'] for c in c2):.0%}, above null 95th pct "
                          f"{above if q95 is not None else '-'}/{len(c2)}")
            ex = max(cs, key=lambda c: c['acc'])
            L_.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:56]}")
            L_.append(f"                                   found   {ex['text'][:56]}")
    L_ += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L_) + '\n'


# ============================================================== self-test
def selftest(args, pool, langs, held, words):
    print(f"\n{S.stamp()} SELF-TEST")
    print(f"cipher: 87 symbols, uncertain positions (1-based): {[q + 1 for q in sorted(ALTS)]}")
    for L in langs:
        w = words[L]
        letters = ''.join(w)
        top = ''.join(c for c, _ in Counter(letters).most_common(10))
        print(f"\n  {L} {NAME[L]:<20} {FILE[L]}: {len(w):,} words, {len(letters):,} letters; "
              f"commonest letters {top}")
        print(f"      start: {' '.join(w[:14])}")
        print(f"      end:   {' '.join(w[-14:])}")
        tail = ' '.join(w[-3000:])
        if 'GUTENBERG' in tail or 'LICENSE' in tail or 'EBOOK' in tail:
            print("      WARNING: the end of this text still looks like Project Gutenberg's English licence.")
        if len(letters) < 200000:
            print("      WARNING: short text (under 200,000 letters); a longer book gives a better model.")
    print("\nOne planted control per language (2 misread symbols, held-back text) and one "
          "pair-preserving copy, through the full search:")
    per = {}
    P = dict(restarts=args.restarts, iters=args.iters)
    for L in langs:
        seq, alts, plain, order = plant(held[L], 5, 2, S.derive_seed('selftest', L))
        t0 = time.time()
        res = run_jobs(pool, [(L, seq, alts, o, P['restarts'], P['iters'], S.derive_seed('st', L, o))
                              for o in ORDERS])
        dt = time.time() - t0
        bi = max(range(5), key=lambda i: res[i][0])
        acc = key_accuracy(res[bi][1], ORDERS[bi], plain, order)
        ns, na = S.null_tokens('doublet', SEQ, ALTS, S.derive_seed('selftest-null', L))
        nres = run_jobs(pool, [(L, ns, na, o, P['restarts'], P['iters'], S.derive_seed('stn', L, o))
                               for o in ORDERS])
        per[L] = dt
        print(f"  {L} {NAME[L]:<20} control {res[bi][0]:.3f} vs copy {max(r[0] for r in nres):.3f}; "
              f"order {'right' if ORDERS[bi] == order else 'different'}, symbols decoded right {acc:.0%}, "
              f"{H.fmt_secs(dt)} per search")
        print(f"      planted {plain[:60]}")
        print(f"      found   {res[bi][1][:60]}")
    tot = 0.0
    print("\nEstimated time to reach the targets:")
    for L in langs:
        n = 2 * args.copies + args.controls + args.real_seeds
        t = n * per[L] * 5 / max(5, args.workers)
        tot += t
        print(f"  {L}: {n} searches, about {H.fmt_secs(t)}")
    print(f"  total about {H.fmt_secs(tot)}")


# ============================================================== main
def main():
    global SEQ, ALTS, LANGS, NEW, BONF
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='languages_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--real-seeds', type=int, default=9)
    ap.add_argument('--restarts', type=int, default=12)
    ap.add_argument('--iters', type=int, default=20000)
    ap.add_argument('--copies', type=int, default=1000, help='copies of each kind per language')
    ap.add_argument('--controls', type=int, default=60)
    ap.add_argument('--only', nargs='+', choices=[c for c, _, _ in LANGS])
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true')
    ap.add_argument('--confirm-latin', action='store_true',
                    help='the pre-registered Latin confirmation: Church Latin (latin_church.txt) only, '
                         'its own scrambled copies and seeds, results in latin_confirm_out')
    args = ap.parse_args()
    if args.confirm_latin:
        # PRE-REGISTERED (9 Oct 2026, before running): confirms the LA result only if the real cipher
        # beats the pair-preserving copies at p < 0.05 in the Church Latin model AND the decrypt reads
        # as Latin. Clementine Vulgate Psalms + four Gospels; no text shared with latin.txt.
        LANGS, NEW, BONF = [CONFIRM], ['LC'], 0.05
        if args.out == 'languages_out':
            args.out = 'latin_confirm_out'

    SEQ, ALTS = S.load_cipher()
    S.SEQ, S.ALTS = SEQ, ALTS
    langs = [c for c, f, _ in LANGS if os.path.exists(f) and (not args.only or c in args.only)]
    missing = [f for c, f, _ in LANGS if not os.path.exists(f)]
    print(f"cipher read from Dorabella_transcription.csv ({len(ALTS)} uncertain symbols)")
    print(f"languages: {', '.join(f'{c} ({NAME[c]})' for c in langs)}"
          + (f"; missing (left out): {', '.join(missing)}" if missing else ''))
    if not langs:
        raise SystemExit("no language texts found next to this script")
    words, held, wl = {}, {}, {}
    for L in langs:
        w, train, ho = load_language(L)
        words[L], held[L], wl[L] = w, ho, WordList(w)
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
        if state['seq'] != SEQ or {int(k): v for k, v in state['alts'].items()} != ALTS:
            raise SystemExit("saved state was made with a different transcription; use --fresh or another --out")
        for L in langs:                         # a language file added later joins the run
            if L not in state['st']:
                state['st'][L] = dict(real=[], shuffle=[], doublet=[], controls=[])
                state['targets'][L] = dict(shuffle=args.copies, doublet=args.copies, controls=args.controls)
    else:
        state = blank_state(args, langs)
    for L in langs:
        size = len(''.join(words[L]))
        old = state['files'].get(L)
        if old and old['letters'] != size and not args.report_only:
            raise SystemExit(f"{FILE[L]} has changed since this run started ({old['letters']:,} -> {size:,} "
                             "letters); use --fresh or another --out")
        state['files'][L] = dict(letters=size)
    if args.report_only:
        S.write_text(rp, build_report(state, langs, wl))
        print(f"report written to {rp}")
        return
    workers = max(1, args.workers)
    S.keep_awake()
    pool = mp.Pool(workers, initializer=init_worker) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, langs, held, words)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}, languages {' '.join(langs)}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stages = stage_list(state, langs)
            if not stages:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for L in langs:
                    for kind in ('shuffle', 'doublet'):
                        state['targets'][L][kind] = int(state['targets'][L][kind] * 1.5)
                log("all targets met; raising the copy targets by 50% and carrying on")
                continue
            _, _, L, kind = stages[0]
            n = batch_size(state, L, kind, workers)
            dt = do_batch(pool, state, held, L, kind, n)
            t1 = time.time()
            state['elapsed'] += t1 - last
            last = t1
            st = state['st'][L]
            msg = f"{L} {kind} x{n} in {H.fmt_secs(dt)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][L][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{S.stamp()} {msg}", flush=True)
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, langs, wl))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            S.save_json(sp, state)
            S.write_text(rp, build_report(state, langs, wl))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


if __name__ == '__main__':
    mp.freeze_support()
    main()
