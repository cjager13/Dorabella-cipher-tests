#!/usr/bin/env python3
"""
Dorabella "private English" battery: is the message in English as Elgar might have
WRITTEN it privately, rather than as a dictionary spells it?

Needs Dorabella_hillclimb.py, Dorabella_marathon.py, Dorabella_tests.py,
Dorabella_transcription.csv and english.txt in the same folder.

  python Dorabella_styles.py --selftest        # ~20 min: checks, power preview, time estimate
  python Dorabella_styles.py                   # start, or resume where it left off
  python Dorabella_styles.py --hours 48        # stop by itself after 48 hours
  python Dorabella_styles.py --report-only     # rewrite the report from saved state
  python Dorabella_styles.py --only S1 S2      # work only on these styles
  python Dorabella_styles.py --fresh           # discard saved state and start over

Ctrl+C any time; the same command resumes. Progress is saved after every batch to
styles_out/state.json, and styles_out/report.txt is rewritten each time. When every
target is met it raises the targets by 50% and keeps going until stopped.

PRE-REGISTERED STYLES (frozen before any real run; 6 styles, so the Bonferroni line
for one style is p < 0.05/6 = 0.0083)

  S0 STANDARD   ordinary spelling (the reference; repeats earlier tests, but with the
                bigger dictionary, the uncertain symbols solved in the search, and the
                same settings as the other styles, so all six are comparable)
  S1 VOWELLESS  first letter kept, later A E I O U dropped: DEAR MISS PENNY -> DR MSS PNNY
  S2 BACKWARDS  every word spelled backwards, word order kept (Victorian back-slang):
                RAED SSIM YNNEP
  S3 PHONETIC   simple sound-spelling: PH->F, C->K or S, CK->K, QU->KW, X->KS, WH->W,
                GH silent, -TION->SHUN, KN/WR -> N/R, double letters single, silent
                final E dropped, final Y->I:  DEAR MIS PENI, THANKS FOR THE NAIT
  S4 INITIALS   first letter of each word only: 87 symbols = 87 words
  S5 PIG LATIN  "hog Latin", known in the 19th century: EARDAY ISSMAY ENNYPAY

For each style a quadgram model is trained on Dracula (english.txt) rewritten in that
style, and its dictionary (about 12,000 words: Dracula's vocabulary + 9,474 common
English words + Elgar names) is rewritten the same way. The search is ANY
one-symbol-one-letter key, by simulated annealing, over the 5 reading orders, with the
9 uncertain symbols (from the CSV) chosen inside the search.

METHOD (unchanged): statistic = best quadgram score over the 5 orders; real = median
over several search seeds. Each plain shuffle and pair-preserving shuffle (keeps every
adjacent symbol pair) goes through the IDENTICAL search. Planted English written in
each style (0 or 2 misread symbols, and the uncertain symbols' alternates set as in
the real cipher) measures power. The word score (longer words worth more; not
optimised) is a second, independent check. p = (1+k)/(n+1).
"""
import argparse
import base64
import csv
import datetime
import json
import math
import multiprocessing as mp
import os
import random
import re
import signal
import statistics
import sys
import time
import zlib
from collections import Counter

import Dorabella_hillclimb as H
import Dorabella_marathon as M
import Dorabella_tests as T

VERSION = 1
NL = 24
LET = H.LETTERS
LIDX = H.LIDX
STYLES = ['S0', 'S1', 'S2', 'S3', 'S4', 'S5']
STYLE_NAME = {'S0': 'standard spelling', 'S1': 'vowels dropped', 'S2': 'words spelled backwards',
              'S3': 'phonetic spelling', 'S4': 'first letters of words', 'S5': 'pig Latin'}
BONF = 0.05 / len(STYLES)
ORDERS = ['forward', 'fully reversed', 'each line reversed', 'boustrophedon A', 'boustrophedon B']
ELGAR_WORDS = """ELGAR EDWARD EDOO ALICE DORA DORABELLA PENNY MALVERN FORLI WORCESTER WOLVERHAMPTON
RECTORY BICYCLE CARICE JAEGER NIMROD TROYTE SHERIDAN WINIFRED NORBURY ISABEL FITTON BAKER
POWELL ENIGMA VARIATION CANTATA ORATORIO CAPTAIN KING OLAF CARACTACUS SAXON BAVARIAN
CHORUS CHOIR MUSIC SINGER CONCERT REHEARSAL VIOLIN PIANO ORGAN SCORE""".split()

QL = None          # quadgram table (list) for the active style in this process
QLS = {}           # style -> table (built lazily per process)
SEQ = None
ALTS = None        # {0-based position: alternate code}
KEEP_AWAKE_OK = False


# ============================================================== styles
VOW = set('AEIOU')


def st_vowelless(w):
    return w[0] + ''.join(c for c in w[1:] if c not in VOW)


def st_phonetic(w):
    s = w
    if s.startswith('KN'):
        s = s[1:]
    if s.startswith('WR'):
        s = s[1:]
    if s.startswith('WH'):
        s = 'W' + s[2:]
    s = s.replace('TION', 'SHUN').replace('SION', 'SHUN')
    s = s.replace('PH', 'F').replace('CK', 'K').replace('QU', 'KW').replace('X', 'KS')
    s = s[:1] + s[1:].replace('GH', '')
    out = []
    for i, c in enumerate(s):
        if c == 'C':
            nxt = s[i + 1] if i + 1 < len(s) else ''
            out.append('S' if nxt in ('E', 'I', 'Y') else ('' if nxt == 'H' else 'K'))
            if nxt == 'H':
                out.append('C')          # keep CH as CH
        else:
            out.append(c)
    s = ''.join(out)
    d = [s[0]] if s else []
    for c in s[1:]:
        if c != d[-1]:
            d.append(c)
    s = ''.join(d)
    if len(s) > 3 and s.endswith('E') and s[-2] not in VOW:
        s = s[:-1]
    if len(s) > 2 and s.endswith('Y') and s[-2] not in VOW:
        s = s[:-1] + 'I'
    return s or w


def st_piglatin(w):
    if w[0] in VOW:
        return w + 'WAY'
    i = 0
    while i < len(w) and not (w[i] in VOW or (w[i] == 'Y' and i > 0)):
        if w[i] == 'Q' and i + 1 < len(w) and w[i + 1] == 'U':
            i += 2
            continue
        i += 1
    return w[i:] + w[:i] + 'AY'


STYLE_FN = {
    'S0': lambda w: w,
    'S1': st_vowelless,
    'S2': lambda w: w[::-1],
    'S3': st_phonetic,
    'S4': lambda w: w[0],
    'S5': st_piglatin,
}


def fold(s):
    return s.replace('J', 'I').replace('V', 'U')


def words_of_text(raw):
    return re.findall(r'[A-Z]+', raw.upper())


def styled_text(words, style):
    f = STYLE_FN[style]
    return fold(''.join(f(w) for w in words))


# ============================================================== language models
CORPUS_WORDS = None


def corpus_words(corpus):
    global CORPUS_WORDS
    if CORPUS_WORDS is None:
        with open(corpus, encoding='utf-8', errors='ignore') as fh:
            CORPUS_WORDS = words_of_text(fh.read())
    return CORPUS_WORDS


def build_quad(corpus, style):
    text = styled_text(corpus_words(corpus), style)
    li = [LIDX[c] for c in text]
    cnt = Counter()
    for i in range(len(li) - 3):
        cnt[((li[i] * NL + li[i + 1]) * NL + li[i + 2]) * NL + li[i + 3]] += 1
    tot = sum(cnt.values())
    floor = math.log(0.01 / tot)
    q = [floor] * (NL ** 4)
    for k, v in cnt.items():
        q[k] = math.log(v / tot)
    return q


def use_style(corpus, style):
    global QL
    if style not in QLS:
        QLS[style] = build_quad(corpus, style)
    QL = QLS[style]


def score_letters(p):
    Q = QL
    t = 0.0
    for s in range(len(p) - 3):
        t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
    return t / (len(p) - 3)


# ============================================================== dictionaries
WORDBLOB = None    # filled in below the main code


def base_dictionary(corpus, extra=None):
    cw = Counter(corpus_words(corpus))
    words = Counter({w: c for w, c in cw.items() if c >= 2 and len(w) >= 2})
    common = zlib.decompress(base64.b64decode(WORDBLOB)).decode().split()
    for w in common + ELGAR_WORDS + H.EMBEDDED_WORDS:
        words[w] += 1
    if extra:
        with open(extra, encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                w = line.strip().upper()
                if w.isalpha() and w.isascii():
                    words[w] += 1
    return words


class StyleDict:
    """Dictionary in a style: styled form -> the original word it most likely is."""

    def __init__(self, base, style):
        self.style = style
        self.gloss = {}
        best = {}
        f = STYLE_FN[style]
        for w, c in base.items():
            if len(w) < 2:
                continue
            s = fold(f(w))
            if not (3 <= len(s) <= 16):
                continue
            if s not in best or c > best[s]:
                best[s], self.gloss[s] = c, w
        self.words = set(self.gloss)
        self.prefixes = {w[:k] for w in self.words for k in range(1, len(w) + 1)}
        self.maxlen = max((len(w) for w in self.words), default=0)
        self.usable = style != 'S4'

    def _dp(self, t):
        n = len(t)
        best, choice = [0] * (n + 1), [0] * (n + 1)
        for i in range(n - 1, -1, -1):
            b, ch = best[i + 1], 0
            for L in range(1, min(self.maxlen, n - i) + 1):
                sub = t[i:i + L]
                if sub not in self.prefixes:
                    break
                if L >= 3 and sub in self.words:
                    v = (L - 2) ** 2 + best[i + L]
                    if v > b:
                        b, ch = v, L
            best[i], choice[i] = b, ch
        return best, choice

    def score(self, t):
        if not self.usable:
            return 0
        return self._dp(t)[0][0]

    def segment(self, t):
        if not self.usable:
            return '(no word check for initials)'
        _, choice = self._dp(t)
        out, gap, i = [], [], 0
        while i < len(t):
            if choice[i]:
                if gap:
                    out.append(''.join(gap).lower())
                    gap = []
                w = t[i:i + choice[i]]
                g = self.gloss[w]
                out.append(w if fold(g) == w else f"{w}[{g.lower()}]")
                i += choice[i]
            else:
                gap.append(t[i])
                i += 1
        if gap:
            out.append(''.join(gap).lower())
        return ' '.join(out)


# ============================================================== cipher data and orders
def load_cipher(path='Dorabella_transcription.csv'):
    if not os.path.exists(path):
        raise SystemExit("Dorabella_transcription.csv not found next to this script.")
    with open(path, newline='', encoding='utf-8') as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r['position']))
    seq = [int(r['code']) for r in rows]
    alts = {int(r['position']) - 1: int(r['alternate_code'])
            for r in rows if r['sure'].strip().lower() == 'no' and r['alternate_code'].strip()}
    assert len(seq) == 87
    return seq, alts


def slot_of(code):
    return (code // 10 - 1) * 8 + (code % 10 - 1)


def code_of_slot(s):
    return (s // 8 + 1) * 10 + (s % 8 + 1)


ORDER_PERM = {o: T.apply_order(list(range(87)), o) for o in ORDERS}


def neighbour(code, rng):
    loop, d = divmod(code, 10)
    return loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)


# ============================================================== annealer
SA_T0, SA_T1, PFLIP = 10.0, 0.3, 0.05


def anneal(c, A, restarts, iters, seed):
    """c: symbol slots in reading order; A[q]: alternate slot or -1.
    Any one-to-one key; uncertain symbols are flipped between readings inside the search."""
    Q = QL
    rng = random.Random(seed)
    rand, rr, exp = rng.random, rng.randrange, math.exp
    n = len(c)
    amb = [q for q in range(n) if A[q] >= 0]
    fac = (SA_T1 / SA_T0) ** (1.0 / iters)
    best_all = (-1e18, None, None)

    def qs(p, ss):
        t = 0.0
        for s in ss:
            t += Q[((p[s] * 24 + p[s + 1]) * 24 + p[s + 2]) * 24 + p[s + 3]]
        return t

    def starts_of(lst):
        return tuple(sorted({s for q in lst for s in range(max(0, q - 3), min(q, n - 4) + 1)}))

    allst = range(n - 3)
    for _ in range(restarts):
        e = list(c)
        for q in amb:
            if rand() < 0.5:
                e[q] = A[q]
        pos = [[] for _ in range(24)]
        for q, x in enumerate(e):
            pos[x].append(q)
        starts = [starts_of(pos[x]) for x in range(24)]
        pc = {}
        key = list(range(24))
        rng.shuffle(key)
        p = [key[x] for x in e]
        cur = qs(p, allst)
        best, bp, be = cur, p[:], e[:]
        Tm = 10.0
        for _ in range(iters):
            if amb and rand() < PFLIP:
                q = amb[rr(len(amb))]
                old = e[q]
                new = c[q] if old == A[q] else A[q]
                ss = range(max(0, q - 3), min(q, n - 4) + 1)
                o = qs(p, ss)
                p[q] = key[new]
                d = qs(p, ss) - o
                if d >= 0 or rand() < exp(d / Tm):
                    e[q] = new
                    pos[old].remove(q)
                    pos[new].append(q)
                    starts[old] = starts_of(pos[old])
                    starts[new] = starts_of(pos[new])
                    pc.clear()
                    cur += d
                    if cur > best:
                        best, bp, be = cur, p[:], e[:]
                else:
                    p[q] = key[old]
                Tm *= fac
                continue
            i = e[rr(n)]
            j = rr(23)
            if j >= i:
                j += 1
            a, b = (i, j) if i < j else (j, i)
            ss = pc.get(a * 24 + b)
            if ss is None:
                ss = tuple(sorted(set(starts[a]) | set(starts[b])))
                pc[a * 24 + b] = ss
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
                    best, bp, be = cur, p[:], e[:]
            else:
                for q in pos[i]:
                    p[q] = li
                for q in pos[j]:
                    p[q] = lj
            Tm *= fac
        assert abs(best - qs(bp, allst)) < 1e-6, "annealer drifted"
        if best > best_all[0]:
            best_all = (best, bp, be)
    s, bp, be = best_all
    flips = [q for q in amb if be[q] == A[q]]
    return s / (n - 3), ''.join(LET[x] for x in bp), flips


# ============================================================== jobs
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    corpus_words(corpus)


def order_job(corpus, style, seq, alts, order, restarts, iters, seed):
    """One reading order of one search. seq: codes in written order; alts: {pos: code}."""
    use_style(corpus, style)
    R = ORDER_PERM[order]
    c = [slot_of(seq[r]) for r in R]
    A = [slot_of(alts[r]) if r in alts else -1 for r in R]
    s, text, flips = anneal(c, A, restarts, iters, seed)
    used = sorted(R[q] + 1 for q in flips)        # 1-based written positions read as alternate
    return s, text, used


def derive_seed(*parts):
    return random.Random('|'.join(str(x) for x in parts)).getrandbits(48)


def run_jobs(pool, argsets):
    if pool is None:
        return [order_job(*a) for a in argsets]
    res = pool.starmap_async(order_job, argsets, chunksize=1)
    while True:
        try:
            return res.get(60)
        except mp.TimeoutError:
            continue


def null_tokens(kind, seq, alts, seed):
    rng = random.Random(seed)
    toks = [(c, alts.get(i)) for i, c in enumerate(seq)]
    if kind == 'shuffle':
        rng.shuffle(toks)
    else:
        new = M.doublet_shuffle(list(seq), rng)
        out = [(c, None) for c in new]
        for code in {c for c, a in toks if a is not None}:
            mine = [a for c, a in toks if c == code and a is not None]
            slots = [i for i, x in enumerate(new) if x == code]
            for i, a in zip(rng.sample(slots, len(mine)), mine):
                out[i] = (code, a)
        toks = out
    return [c for c, _ in toks], {i: a for i, (_, a) in enumerate(toks) if a is not None}


def plant(style, k, noise, seed):
    """Planted English in a style: random key, random reading order, misreads, and
    alternates at the uncertain positions (the true reading is main or alternate at random)."""
    rng = random.Random(seed)
    words = []
    j = k
    while len(styled_text(words, style)) < 87:
        words += words_of_text(M.CONTROL_TEXTS[j % len(M.CONTROL_TEXTS)])
        j += 7
    plain = styled_text(words, style)[:87]
    slots = list(range(24))
    rng.shuffle(slots)
    order = rng.choice(ORDERS)
    R = ORDER_PERM[order]
    seq = [0] * 87
    for i, ch in enumerate(plain):
        seq[R[i]] = code_of_slot(slots[LIDX[ch]])
    unc = sorted(ALTS)
    for q in rng.sample([i for i in range(87) if i not in ALTS], noise):
        seq[q] = neighbour(seq[q], rng)
    alts = {}
    for q in unc:
        wrong = neighbour(seq[q], rng)
        if rng.random() < 0.5:
            alts[q] = wrong
        else:
            alts[q], seq[q] = seq[q], wrong
    return seq, alts, plain, order


# ============================================================== state and batches
def now():
    return datetime.datetime.now()


def stamp():
    return f"[{now():%Y-%m-%d %H:%M:%S}]"


def keep_awake():
    global KEEP_AWAKE_OK
    if sys.platform == 'win32':
        try:
            import ctypes
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
            KEEP_AWAKE_OK = True
        except Exception:
            pass


def save_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(obj, fh)
    os.replace(tmp, path)


def write_text(path, text):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        fh.write(text)
    os.replace(tmp, path)


def blank_state(args):
    return dict(version=VERSION, created=str(now()), elapsed=0.0, sessions=0,
                seq=SEQ, alts={str(k): v for k, v in ALTS.items()},
                params=dict(seed=args.seed, real_seeds=args.real_seeds,
                            restarts=args.restarts, iters=args.iters),
                targets={s: dict(shuffle=args.copies, doublet=args.copies, controls=args.controls)
                         for s in STYLES},
                st={s: dict(real=[], shuffle=[], doublet=[], controls=[]) for s in STYLES},
                log=[])


def stage_list(state, only):
    out = []
    for s in (only or STYLES):
        st, tg = state['st'][s], state['targets'][s]
        if len(st['real']) < state['params']['real_seeds']:
            out.append((0, 0.0, s, 'real'))
            continue
        if len(st['controls']) < min(10, tg['controls']):
            out.append((1, len(st['controls']) / tg['controls'], s, 'controls'))
        for kind in ('doublet', 'shuffle', 'controls'):
            frac = len(st[kind]) / tg[kind]
            if frac < 1:
                out.append((2, frac, s, kind))
    return sorted(out)


def do_batch(pool, state, corpus, dicts, style, kind, count):
    P = state['params']
    st = state['st'][style]
    items = []
    for k in range(count):
        n = len(st[kind]) + k
        if kind == 'real':
            items.append((SEQ, ALTS, derive_seed(P['seed'], 'real', style, n), None))
        elif kind in ('shuffle', 'doublet'):
            s, a = null_tokens(kind, SEQ, ALTS, derive_seed(P['seed'], 'null', style, kind, n))
            items.append((s, a, derive_seed(P['seed'], 'nullsearch', style, kind, n), None))
        else:
            noise = 0 if n % 2 == 0 else 2
            s, a, plain, order = plant(style, n, noise, derive_seed(P['seed'], 'plant', style, n))
            items.append((s, a, derive_seed(P['seed'], 'ctrlsearch', style, n), (noise, plain, order)))
    jobs = [(corpus, style, s, a, o, P['restarts'], P['iters'], derive_seed(sd, o))
            for s, a, sd, _ in items for o in ORDERS]
    t0 = time.time()
    res = run_jobs(pool, jobs)
    d = dicts[style]
    for k, (s, a, sd, meta) in enumerate(items):
        per = res[k * len(ORDERS):(k + 1) * len(ORDERS)]
        bi = max(range(len(ORDERS)), key=lambda i: per[i][0])
        stat, text, used = per[bi]
        w = d.score(text)
        if kind == 'real':
            ranked = sorted(range(len(ORDERS)), key=lambda i: -per[i][0])
            st['real'].append(dict(stat=stat, words=w, top=[
                dict(order=ORDERS[i], score=per[i][0], text=per[i][1], alts=per[i][2]) for i in ranked]))
        elif kind in ('shuffle', 'doublet'):
            st[kind].append(dict(stat=round(stat, 5), words=w))
        else:
            noise, plain, order = meta
            ref = ''.join(plain)
            acc = sum(x == y for x, y in zip(text, ref)) / 87
            st['controls'].append(dict(stat=round(stat, 5), words=w, noise=noise,
                                       acc=round(acc, 3), right_order=int(ORDERS[bi] == order),
                                       text=text, plain=plain))
    return time.time() - t0


def batch_size(state, style, kind, workers):
    if kind == 'real':
        return state['params']['real_seeds'] - len(state['st'][style]['real'])
    n = max(1, math.ceil(2 * workers / len(ORDERS)))
    return max(1, min(n, state['targets'][style][kind] - len(state['st'][style][kind])))


# ============================================================== report
def pval(real, nulls):
    if real is None or not nulls:
        return None
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


def pct(vals, q):
    v = sorted(vals)
    return v[min(len(v) - 1, int(q * len(v)))] if v else None


def f(x, nd=3):
    return '-' if x is None else f"{x:.{nd}f}"


def verdict(p_d, p_s, power):
    if p_d is None:
        return 'waiting for copies'
    if p_d <= BONF and (p_s or 1) <= BONF:
        return ('BEATS both nulls past Bonferroni: read the decrypts' +
                (' (low power: be careful)' if power is not None and power < 0.5 else ''))
    if p_d <= 0.05:
        return 'below 0.05, not past Bonferroni (expected by chance among 6 styles)'
    if power is not None and power < 0.5:
        return 'chance level, LOW POWER: inconclusive'
    return 'chance level'


def build_report(state, dicts):
    L = ["DORABELLA 'PRIVATE ENGLISH' STYLES REPORT",
         f"written {now():%Y-%m-%d %H:%M:%S}   compute time {H.fmt_secs(state['elapsed'])}   "
         f"sessions {state['sessions']}",
         f"settings: real seeds {state['params']['real_seeds']} (median used), annealing "
         f"{state['params']['restarts']} restarts x {state['params']['iters']} steps per reading order, "
         f"5 reading orders, {len(ALTS)} uncertain symbols solved in the search",
         f"Bonferroni line for one of 6 styles: p < {BONF:.4f}", "",
         "=" * 104,
         "HEADLINE (statistic = best quadgram score in that style's own model; real = median over seeds;",
         "          p_words = word score vs pair-preserving copies; power = planted controls above the",
         "          95th percentile of the pair-preserving copies)",
         "=" * 104,
         f"{'style':<26} {'real':>8} {'n_pair':>6} {'p_pair':>7} {'n_shuf':>6} {'p_shuf':>7} "
         f"{'p_words':>7} {'power':>6} {'ctrl med':>8}  verdict"]
    detail = []
    for s in STYLES:
        st = state['st'][s]
        real = [r['stat'] for r in st['real']]
        med = statistics.median(real) if real else None
        dn = [x['stat'] for x in st['doublet']]
        sn = [x['stat'] for x in st['shuffle']]
        dw = [x['words'] for x in st['doublet']]
        rw = statistics.median(r['words'] for r in st['real']) if real else None
        ctrl = [c['stat'] for c in st['controls']]
        q95 = pct(dn, 0.95)
        power = sum(1 for c in ctrl if c > q95) / len(ctrl) if (ctrl and q95 is not None) else None
        p_d, p_s = pval(med, dn), pval(med, sn)
        p_w = pval(rw, dw) if dicts[s].usable else None
        L.append(f"{s + ' ' + STYLE_NAME[s]:<26} {f(med):>8} {len(dn):>6} {f(p_d, 4):>7} {len(sn):>6} "
                 f"{f(p_s, 4):>7} {f(p_w, 3):>7} {('-' if power is None else f'{power:.0%}'):>6} "
                 f"{f(statistics.median(ctrl) if ctrl else None, 2):>8}  {verdict(p_d, p_s, power)}")
        detail.append((s, st, med, dn, sn, dw, rw, q95))
    L += ["", "A style matters only if the real cipher beats the PAIR-PRESERVING copies past the",
          "Bonferroni line, its controls show it can find a message (power near 100%), and the",
          "decrypt reads as English written in that style. Words are shown as STYLED[original]."]
    for s, st, med, dn, sn, dw, rw, q95 in detail:
        d = dicts[s]
        L += ["", "-" * 104, f"{s}: {STYLE_NAME[s]}", "-" * 104]
        if not st['real']:
            L.append("real run not done yet")
            continue
        L.append("real statistic per seed: " + ', '.join(f"{r['stat']:.3f}" for r in st['real'])
                 + f"   median {med:.3f};  word score per seed: "
                 + ', '.join(str(r['words']) for r in st['real']))
        for name, nul in (('pair-preserving', dn), ('plain shuffle', sn)):
            if nul:
                L.append(f"  vs {name:<16} n={len(nul)}, mean {statistics.mean(nul):.3f}, 95th pct "
                         f"{pct(nul, .95):.3f}, max {max(nul):.3f}, p = {pval(med, nul):.4f}")
        if dw and d.usable:
            L.append(f"  word score vs pair-preserving: real median {rw}, copies mean "
                     f"{statistics.mean(dw):.1f}, 95th pct {pct(dw, .95)}, p = {pval(rw, dw):.4f}")
        b = max(st['real'], key=lambda r: r['stat'])
        L.append(f"best real decrypts (from the seed with the highest statistic, {b['stat']:.3f}):")
        for t in b['top'][:3]:
            alt = f"  alternate readings used at {t['alts']}" if t['alts'] else ''
            L.append(f"  [{t['score']:.3f}] {t['order']}{alt}")
            for a in range(0, 87, 29):
                L.append("      " + t['text'][a:a + 29])
            L.append("      words: " + d.segment(t['text']))
        cs = st['controls']
        if cs:
            for nz in (0, 2):
                c2 = [c for c in cs if c['noise'] == nz]
                if not c2:
                    continue
                above = sum(1 for c in c2 if q95 is not None and c['stat'] > q95)
                L.append(f"controls, {nz} misread symbols: n={len(c2)}, median stat "
                         f"{statistics.median(c['stat'] for c in c2):.3f}, right reading order "
                         f"{sum(c['right_order'] for c in c2)}/{len(c2)}, median letters right "
                         f"{statistics.median(c['acc'] for c in c2):.0%}, above null 95th pct "
                         f"{above if q95 is not None else '-'}/{len(c2)}")
            ex = max(cs, key=lambda c: c['acc'])
            L.append(f"  example control ({ex['acc']:.0%} right): planted {ex['plain'][:44]}")
            L.append(f"                                  found   {ex['text'][:44]}")
    L += ["", "Run log (last 20):"] + ["  " + x for x in state['log'][-20:]]
    return '\n'.join(L) + '\n'


# ============================================================== selftest
def selftest(args, pool, dicts, workers):
    print(f"\n{stamp()} SELF-TEST")
    sample = words_of_text("Dear Miss Penny thanks for the night we arrived home quite safely")
    for s in STYLES:
        print(f"  {s} {STYLE_NAME[s]:<24} {fold(' '.join(STYLE_FN[s](w) for w in sample))}")
    print(f"\ncipher: {len(SEQ)} symbols, uncertain positions (1-based): {[q + 1 for q in sorted(ALTS)]}")
    for s in STYLES:
        print(f"  {s} dictionary: {len(dicts[s].words):,} styled words")
    d = M.doublet_shuffle(list(SEQ), random.Random(1))
    assert Counter(zip(d, d[1:])) == Counter(zip(SEQ, SEQ[1:]))
    print("pair-preserving shuffle check passed")
    print("\nOne planted control per style (2 misread symbols, alternates as in the real cipher), "
          "and one pair-preserving copy, through the full search:")
    P = dict(seed=args.seed, restarts=args.restarts, iters=args.iters)
    per = {}
    for s in (args.only or STYLES):
        seq, alts, plain, order = plant(s, 5, 2, derive_seed('selftest', s))
        t0 = time.time()
        jobs = [(args.corpus, s, seq, alts, o, P['restarts'], P['iters'], derive_seed('st', s, o)) for o in ORDERS]
        res = run_jobs(pool, jobs)
        dt = time.time() - t0
        bi = max(range(5), key=lambda i: res[i][0])
        acc = sum(x == y for x, y in zip(res[bi][1], plain)) / 87
        ns, na = null_tokens('doublet', SEQ, ALTS, derive_seed('selftest-null', s))
        nres = run_jobs(pool, [(args.corpus, s, ns, na, o, P['restarts'], P['iters'], derive_seed('stn', s, o))
                               for o in ORDERS])
        per[s] = dt
        print(f"  {s} {STYLE_NAME[s]:<24} control {res[bi][0]:.3f} vs copy {max(r[0] for r in nres):.3f}; "
              f"order {'right' if ORDERS[bi] == order else 'WRONG'}, letters right {acc:.0%}, "
              f"{H.fmt_secs(dt)} per search")
        print(f"      planted {plain[:60]}")
        print(f"      found   {res[bi][1][:60]}")
        print(f"      words   {dicts[s].segment(res[bi][1])[:100]}")
    tot = 0.0
    print("\nEstimated time to reach the targets on this machine (searches of one style run "
          f"{max(1, workers // 5)}+ at a time):")
    for s in (args.only or STYLES):
        n = 2 * args.copies + args.controls + args.real_seeds
        t = n * per[s] * 5 / max(5, workers)
        tot += t
        print(f"  {s}: {n} searches, about {H.fmt_secs(t)}")
    print(f"  total about {H.fmt_secs(tot)}, then it keeps raising targets until you stop it")


# ============================================================== main
def main():
    global SEQ, ALTS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--wordlist', help='optional extra word list, one word per line')
    ap.add_argument('--out', default='styles_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--real-seeds', type=int, default=9)
    ap.add_argument('--restarts', type=int, default=12)
    ap.add_argument('--iters', type=int, default=20000)
    ap.add_argument('--copies', type=int, default=2000, help='copies of each kind per style')
    ap.add_argument('--controls', type=int, default=100)
    ap.add_argument('--only', nargs='+', choices=STYLES)
    ap.add_argument('--fresh', action='store_true')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--stop-when-done', action='store_true',
                    help='exit when the copy targets are met (for running scripts back to back)')
    args = ap.parse_args()

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    SEQ, ALTS = load_cipher()
    print(f"cipher read from Dorabella_transcription.csv ({len(ALTS)} uncertain symbols)")
    os.makedirs(args.out, exist_ok=True)
    sp, rp, lp = (os.path.join(args.out, x) for x in ('state.json', 'report.txt', 'log.txt'))

    def log(msg):
        line = f"{stamp()} {msg}"
        print(line, flush=True)
        with open(lp, 'a', encoding='utf-8') as fh:
            fh.write(line + '\n')

    base = base_dictionary(args.corpus, args.wordlist)
    dicts = {s: StyleDict(base, s) for s in STYLES}

    if os.path.exists(sp) and not args.fresh and not args.selftest:
        with open(sp, encoding='utf-8') as fh:
            state = json.load(fh)
        if state['seq'] != SEQ or {int(k): v for k, v in state['alts'].items()} != ALTS:
            raise SystemExit("saved state was made with a different transcription; use --fresh "
                             "or another --out")
    else:
        state = blank_state(args)

    if args.report_only:
        write_text(rp, build_report(state, dicts))
        print(f"report written to {rp}")
        return

    workers = max(1, args.workers)
    keep_awake()
    pool = mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    try:
        if args.selftest:
            selftest(args, pool, dicts, workers)
            return
        state['sessions'] += 1
        log(f"START session {state['sessions']}, workers {workers}"
            f"{', keep-awake on' if KEEP_AWAKE_OK else ''}")
        last = time.time()
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                log("hour limit reached")
                break
            stages = stage_list(state, args.only)
            if not stages:
                if args.stop_when_done:
                    log("all targets met; stopping (--stop-when-done)")
                    break
                for s in (args.only or STYLES):
                    for kind in ('shuffle', 'doublet'):
                        state['targets'][s][kind] = int(state['targets'][s][kind] * 1.5)
                log("all targets met; raising the copy targets by 50% and carrying on")
                continue
            _, _, s, kind = stages[0]
            n = batch_size(state, s, kind, workers)
            dt = do_batch(pool, state, args.corpus, dicts, s, kind, n)
            t = time.time()
            state['elapsed'] += t - last
            last = t
            st = state['st'][s]
            msg = f"{s} {kind} x{n} in {H.fmt_secs(dt)}"
            msg += (f"  (real stat {st['real'][-1]['stat']:.3f})" if kind == 'real'
                    else f"  ({len(st[kind])}/{state['targets'][s][kind]})")
            state['log'] = (state['log'] + [msg])[-200:]
            print(f"{stamp()} {msg}", flush=True)
            save_json(sp, state)
            write_text(rp, build_report(state, dicts))
    except KeyboardInterrupt:
        print("\nStopped by Ctrl+C. Progress up to the last finished batch is saved; "
              "run the same command to resume.")
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()
        if not args.selftest:
            save_json(sp, state)
            write_text(rp, build_report(state, dicts))
            log(f"END session, {H.fmt_secs(time.time() - t_start)} this session")


# ============================================================== embedded word list
# 9,474 common English words (google-10000-english, no swears; MIT licence), zlib + base64
WORDBLOB = (
    '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'
)


if __name__ == '__main__':
    mp.freeze_support()
    main()
