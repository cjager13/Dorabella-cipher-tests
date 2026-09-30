#!/usr/bin/env python3
"""
Dorabella hill-climbing solver (NO wheel assumption).

  python Dorabella_hillclimb.py --selftest
  python Dorabella_hillclimb.py | Tee-Object hillclimb_out.txt
  python Dorabella_hillclimb.py --corpus D:\\personal\\ciphers\\english.txt
  python Dorabella_hillclimb.py --set 23=14 --set 32=18      (1-based positions)

Model: each cipher symbol maps to one letter of a 24-letter alphabet
(J->I, V->U merged). A random-restart hill climb swaps letter assignments
to maximise an English n-gram score.

Why baselines matter: with 87 symbols a hill climber can always find
gibberish that scores well. Real signal only counts if the REAL symbol order
beats SHUFFLED copies of the same symbols under the same climb. Word score is
NOT optimised, so it is a second, independent check.

--corpus: plain-text English file (a Project Gutenberg book works; 200k+
letters is better). Switches scoring to quadgrams, which is far stronger than
the small embedded bigram/trigram tables.
"""
import argparse
import datetime
import math
import multiprocessing as mp
import os
import random
import statistics
import time
from collections import Counter


# ------------------------------------------------------------------ timer
def fmt_secs(s):
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    if h:
        return f"{int(h)}h {int(m):02d}m {sec:04.1f}s"
    if m:
        return f"{int(m)}m {sec:04.1f}s"
    return f"{sec:.2f}s"


class RunTimer:
    """Start/end timestamps, per-stage laps, and throughput. Call finish() in a finally block."""

    def __init__(self, label="script"):
        self.label = label
        self.t0 = time.perf_counter()
        self.last = self.t0
        self.laps = []
        self.start_dt = datetime.datetime.now()
        print(f"[{self.start_dt:%Y-%m-%d %H:%M:%S}] START {label}", flush=True)

    def lap(self, name, units=None, unit_name="items"):
        now = time.perf_counter()
        dt = now - self.last
        self.last = now
        rate = units / dt if units and dt > 0 else None
        self.laps.append((name, dt, units, unit_name, rate))
        extra = f"  ({units:,} {unit_name} -> {rate:,.0f} {unit_name}/s)" if rate else ""
        print(f"[{datetime.datetime.now():%H:%M:%S}] {name}: {fmt_secs(dt)}{extra}", flush=True)

    def finish(self):
        total = time.perf_counter() - self.t0
        end_dt = datetime.datetime.now()
        print(f"\n[{end_dt:%Y-%m-%d %H:%M:%S}] END {self.label}")
        print(f"Total elapsed: {fmt_secs(total)}")
        rated = [l for l in self.laps if l[4]]
        if rated:
            print("Stage summary:")
            for name, dt, units, un, rate in self.laps:
                r = f"{rate:,.0f} {un}/s" if rate else "-"
                print(f"  {name:<44} {fmt_secs(dt):>10}   {r}")
            best = max(rated, key=lambda l: l[4])
            print(f"Scaling guide (from '{best[0]}'): 1,000,000 {best[3]} "
                  f"~ {fmt_secs(1_000_000 / best[4])} at this rate")


# ------------------------------------------------------------------ cipher data
LINE1 = [23,37,24,33,13,25,11,33,16,22,34,28,18,14,28,35,28,28,25,37,37,28,15,12,22,12,15,34,38]
LINE2 = [11,28,11,25,11,33,16,26,33,24,28,28,24,25,15,18,11,28,34,28,25,21,38,18,14,12,16,16,12,34,38]
LINE3 = [24,38,25,21,38,24,14,21,31,15,38,24,28,25,21,18,38,15,33,37,15,38,25,33,14,13,33]
BASE_SEQ = LINE1 + LINE2 + LINE3
assert len(BASE_SEQ) == 87
LINE_SPANS = [(0, 29), (29, 60), (60, 87)]
CODES_ALL = [loop * 10 + d for loop in (1, 2, 3) for d in range(1, 9)]

LETTERS = 'ABCDEFGHIKLMNOPQRSTUWXYZ'          # 24 letters: J->I, V->U
LIDX = {ch: i for i, ch in enumerate(LETTERS)}
NL = 24
TOL = str.maketrans({'J': 'I', 'V': 'U', 'W': 'U', 'K': 'C', 'Z': 'S'})  # word-check folding


def fold_base(s):
    return s.replace('J', 'I').replace('V', 'U')


# ------------------------------------------------------------------ language data
BIGRAM_PCT = {
    'TH':3.56,'HE':3.07,'IN':2.43,'ER':2.05,'AN':1.99,'RE':1.85,'ND':1.35,
    'ON':1.31,'EN':1.26,'AT':1.12,'OU':1.12,'ED':1.08,'HA':1.07,'TO':1.06,
    'OR':1.06,'IT':1.05,'IS':1.02,'HI':0.93,'ES':0.93,'NG':0.89,'VE':0.83,
    'CO':0.79,'DE':0.76,'RA':0.75,'AL':0.75,'RO':0.73,'LI':0.70,'IC':0.70,
    'LE':0.69,'AR':0.68,'NT':0.66,'IO':0.66,'ST':0.66,'EA':0.65,'ME':0.63,
    'SE':0.63,'FO':0.62,'TE':0.62,'NE':0.62,'WA':0.60,'CA':0.55,'AS':0.55,
    'LA':0.55,'RI':0.55,'WI':0.55,'BE':0.53,'SI':0.51,'DI':0.51,'EL':0.51,
    'UR':0.50,'CH':0.48,'AC':0.48,'TA':0.47,'EC':0.47,'PE':0.46,'SO':0.45,
    'NO':0.45,'LO':0.44,'FA':0.42,'MA':0.42,'PR':0.42,'IL':0.41,'MO':0.40,
    'GE':0.40,'AM':0.40,'PA':0.39,'SS':0.39,'NS':0.38,'US':0.38,'GA':0.37,
    'CT':0.37,'RS':0.36,'WE':0.36,'UT':0.36,'AD':0.36,'OM':0.36,'ID':0.35,
    'ET':0.35,'IE':0.35,'EE':0.33,'IV':0.33,'NC':0.33,'UN':0.33,'EM':0.33,
    'SP':0.32,'GH':0.31,'TI':0.31,'OW':0.30,'AI':0.30,'PL':0.30,'IG':0.29,
    'OS':0.29,'UL':0.28,'OP':0.28,
}
# Approximate percentages, from memory: only the most common English trigrams.
TRIGRAM_PCT = {
    'THE':1.81,'AND':0.73,'ING':0.72,'ENT':0.42,'ION':0.42,'HER':0.36,
    'FOR':0.34,'THA':0.33,'NTH':0.33,'INT':0.32,'ERE':0.31,'TIO':0.31,
    'TER':0.30,'EST':0.28,'ERS':0.28,'ATI':0.26,'HAT':0.26,'ATE':0.25,
    'ALL':0.25,'ETH':0.24,'HES':0.24,'VER':0.24,'HIS':0.24,'OFT':0.22,
    'ITH':0.21,'FTH':0.21,'STH':0.21,'OTH':0.21,'RES':0.21,'ONT':0.20,
}
BI_FLOOR, TRI_FLOOR, W_TRI = 0.02, 0.002, 0.5

EMBEDDED_WORDS = """
THE AND FOR ARE BUT NOT YOU ALL CAN HER WAS ONE OUR OUT DAY GET HAS HIM HIS HOW
MAN NEW NOW OLD SEE TWO WAY WHO BOY DID ITS LET PUT SAY SHE TOO USE ANY MAY HAD
WHY YES YET SIR MRS DOG TEA BIG END FAR FEW GOT OWN RUN SET SUN TRY WIN ASK FUN
JOY SAD WAR ART AGE AIR EAR EYE ARM BAD BED BOX CAR CUT EAT FLY GOD HOT ICE KEY
LAW LIE LOT MAP NET NOR OFF PAY RED SIT SIX TEN TOP VOW WET ZOO THAT WITH HAVE
THIS WILL YOUR FROM THEY KNOW WANT BEEN GOOD MUCH SOME TIME VERY WHEN COME HERE
JUST LIKE LONG MAKE MANY OVER SUCH TAKE THAN THEM WELL WERE WHAT ALSO ONLY THEN
EVEN BACK WORK LIFE MISS DEAR LOVE HOPE WISH SOON HOME NOTE SONG SING TUNE WALK
RIDE STAY WEEK DAYS NIGHT TELL TOLD SAID CODE WORD LAST NEXT MORE MOST MUST SEEN
GONE HEAR HEART SMILE LAUGH HAPPY GLAD KIND FOND DARK LARK SIGH SOUL MIND ABOUT
AFTER AGAIN BEFORE COULD EVERY FIRST GREAT OTHER THEIR THERE THESE THINK THOSE
UNDER WHERE WHICH WHILE WOULD WRITE PLEASE THANK THANKS DEAREST LETTER MUSIC
DANCE CHORUS CHOIR ORGAN PIANO VIOLIN TOMORROW TODAY LONDON MALVERN WORCESTER
WOLVERHAMPTON RECTORY VISIT EVENING MORNING BICYCLE FRIEND FRIENDS KINDLY MISSED
REMEMBER FORGET ALWAYS NEVER OFTEN THOUGHT MERRY SORRY PLEASURE DELIGHT LOVELY
SECRET RIDDLE PUZZLE CIPHER CYPHER MESSAGE MEANING CANNOT MIGHT SHALL SHOULD
BEING DOING GOING COMING THING THINGS NOTHING SOMETHING ELGAR EDOO EDWARD ALICE
DORA DORABELLA PENNY SINGER FOOTBALL MARCO GARDEN SUMMER JULY WEATHER SUNSHINE
RAIN DINNER LUNCH BREAKFAST WITHOUT ABOVE BELOW BEHIND AMONG TOGETHER ANOTHER
BECAUSE BETWEEN THROUGH ENOUGH WONDER WONDERFUL DELIGHTFUL BEAUTIFUL YESTERDAY
GOODBYE FAREWELL SAFELY ARRIVED AGAINST ALMOST ALREADY ENIGMA VARIATION RETURN
LEAVE WAIT WAITING HEARD SPEAK SPOKE STORY TALK TALKING TEASE
""".split()
MIN_WORD, MAX_WORD = 3, 14

PLANT = ("DEAR MISS PENNY WE ARRIVED HOME SAFELY AFTER A PLEASANT JOURNEY AND THINK OF "
         "YOU OFTEN PLEASE COME TO MALVERN SOON AND BRING YOUR BICYCLE TOO")

# per-process globals (set by build_scorer / build_words)
SCORE_FN = None
MODE = ''
BI = TRI = QUAD = None
QFLOOR = 0.0
WORDS = PREFIXES = None
MAXLEN = 0


# ------------------------------------------------------------------ scorers
def score_embedded(p):
    n = len(p)
    b = 0.0
    for i in range(n - 1):
        b += BI[p[i] * NL + p[i + 1]]
    t = 0.0
    for i in range(n - 2):
        t += TRI[(p[i] * NL + p[i + 1]) * NL + p[i + 2]]
    return b / (n - 1) + W_TRI * t / (n - 2)


def score_quad(p):
    t = 0.0
    get = QUAD.get
    for i in range(len(p) - 3):
        t += get(((p[i] * NL + p[i + 1]) * NL + p[i + 2]) * NL + p[i + 3], QFLOOR)
    return t / (len(p) - 3)


def build_scorer(corpus_path=None):
    global SCORE_FN, MODE, BI, TRI, QUAD, QFLOOR
    if corpus_path:
        with open(corpus_path, encoding='utf-8', errors='ignore') as fh:
            raw = fh.read().upper()
        text = fold_base(''.join(c for c in raw if 'A' <= c <= 'Z'))
        if len(text) < 1000:
            raise SystemExit(f"corpus too small ({len(text)} letters)")
        li = [LIDX[c] for c in text]
        counts = Counter()
        for i in range(len(li) - 3):
            counts[((li[i] * NL + li[i + 1]) * NL + li[i + 2]) * NL + li[i + 3]] += 1
        total = sum(counts.values())
        QUAD = {k: math.log(v / total) for k, v in counts.items()}
        QFLOOR = math.log(0.01 / total)
        SCORE_FN = score_quad
        MODE = f"quadgram from corpus ({len(text):,} letters)"
    else:
        bi, tri = Counter(), Counter()
        for pair, f in BIGRAM_PCT.items():
            q = fold_base(pair)
            bi[LIDX[q[0]] * NL + LIDX[q[1]]] += f
        for tg, f in TRIGRAM_PCT.items():
            q = fold_base(tg)
            tri[(LIDX[q[0]] * NL + LIDX[q[1]]) * NL + LIDX[q[2]]] += f
        BI = [math.log(max(bi.get(k, 0.0), BI_FLOOR)) for k in range(NL * NL)]
        TRI = [math.log(max(tri.get(k, 0.0), TRI_FLOOR)) for k in range(NL ** 3)]
        SCORE_FN = score_embedded
        MODE = "embedded bigram+trigram (weak; supply --corpus for real power)"


def build_words(wordlist_path=None):
    global WORDS, PREFIXES, MAXLEN
    raw = set(EMBEDDED_WORDS)
    if wordlist_path:
        with open(wordlist_path, encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                w = line.strip().upper()
                if w.isalpha() and w.isascii():
                    raw.add(w)
    WORDS = {w.translate(TOL) for w in raw if MIN_WORD <= len(w) <= MAX_WORD}
    PREFIXES = set()
    for w in WORDS:
        for k in range(1, len(w) + 1):
            PREFIXES.add(w[:k])
    MAXLEN = max(len(w) for w in WORDS)


def init_worker(corpus_path, wordlist_path):
    build_scorer(corpus_path)
    build_words(wordlist_path)


# ------------------------------------------------------------------ word tools
def word_score(tol_text):
    """Best non-overlapping word cover; a word of length L is worth (L-2)^2."""
    n = len(tol_text)
    best = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        b = best[i + 1]
        for L in range(1, min(MAXLEN, n - i) + 1):
            sub = tol_text[i:i + L]
            if sub not in PREFIXES:
                break
            if L >= MIN_WORD and sub in WORDS:
                v = (L - 2) ** 2 + best[i + L]
                if v > b:
                    b = v
        best[i] = b
    return best[0]


def segment(tol_text):
    n = len(tol_text)
    best, choice = [0] * (n + 1), [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        b, c = best[i + 1], 0
        for L in range(1, min(MAXLEN, n - i) + 1):
            sub = tol_text[i:i + L]
            if sub not in PREFIXES:
                break
            if L >= MIN_WORD and sub in WORDS:
                v = (L - 2) ** 2 + best[i + L]
                if v > b:
                    b, c = v, L
        best[i], choice[i] = b, c
    out, gap, i = [], [], 0
    while i < n:
        if choice[i]:
            if gap:
                out.append(''.join(gap).lower())
                gap = []
            out.append(tol_text[i:i + choice[i]])
            i += choice[i]
        else:
            gap.append(tol_text[i])
            i += 1
    if gap:
        out.append(''.join(gap).lower())
    return ' '.join(out)


# ------------------------------------------------------------------ search
def decrypt_text(cipher_idx, key):
    return ''.join(LETTERS[key[c]] for c in cipher_idx)


def climb(cipher_idx, nsym, rng):
    """Random-start hill climb. key[slot] = letter index; slots >= nsym are unused dummies."""
    score = SCORE_FN
    key = list(range(NL))
    rng.shuffle(key)
    evals = 1
    cur = score([key[c] for c in cipher_idx])
    improved = True
    while improved:
        improved = False
        for i in range(nsym):
            for j in range(i + 1, NL):
                key[i], key[j] = key[j], key[i]
                s = score([key[c] for c in cipher_idx])
                evals += 1
                if s > cur + 1e-12:
                    cur = s
                    improved = True
                else:
                    key[i], key[j] = key[j], key[i]
    return cur, key[:nsym], evals


def run_restarts(seed, cipher_idx, nsym, n):
    rng = random.Random(seed)
    finals, evals = [], 0
    for _ in range(n):
        s, key, e = climb(cipher_idx, nsym, rng)
        finals.append((s, tuple(key)))
        evals += e
    return finals, evals


def shuffle_job(seed, cipher_idx, nsym, n):
    """Baseline: same symbols in random order, same climb, same restarts."""
    rng = random.Random(seed)
    s = list(cipher_idx)
    rng.shuffle(s)
    finals, evals = run_restarts(seed + 7, s, nsym, n)
    best_s, best_key = max(finals, key=lambda t: t[0])
    ws = word_score(decrypt_text(s, best_key).translate(TOL))
    return best_s, ws, evals


def run_parallel(pool, fn, argsets):
    if pool is None:
        return [fn(*a) for a in argsets]
    return pool.starmap(fn, argsets)


def split_counts(total, parts):
    base, extra = divmod(total, parts)
    return [base + (1 if i < extra else 0) for i in range(parts)]


def search_parallel(pool, workers, seed, cipher_idx, nsym, restarts):
    jobs = [(seed + w * 101, cipher_idx, nsym, n)
            for w, n in enumerate(split_counts(restarts, workers)) if n > 0]
    res = run_parallel(pool, run_restarts, jobs)
    finals = [item for f, _ in res for item in f]
    return finals, sum(e for _, e in res)


def prep_cipher(seq):
    uniq = sorted(set(seq))
    idx = {c: i for i, c in enumerate(uniq)}
    return uniq, [idx[c] for c in seq]


def pvalue(real, nulls):
    return (1 + sum(1 for x in nulls if x >= real - 1e-12)) / (1 + len(nulls))


# ------------------------------------------------------------------ reporting
def print_lines(text):
    if len(text) == 87:
        for a, b in LINE_SPANS:
            print("   " + text[a:b])
    else:
        print("   " + text)


def report_best(finals, cipher_idx, uniq, top):
    best_score = max(s for s, _ in finals)
    by_text = {}
    for s, key in finals:
        t = decrypt_text(cipher_idx, key)
        if t not in by_text or s > by_text[t][0]:
            by_text[t] = (s, key)
    ranked = sorted(by_text.items(), key=lambda kv: -kv[1][0])
    n_at_best = sum(1 for s, _ in finals if s >= best_score - 1e-9)
    print(f"\nRestarts: {len(finals)}   distinct final decryptions: {len(ranked)}   "
          f"restarts reaching the best score: {n_at_best}")
    print("(Many restarts hitting one optimum suggests a sharp solution; a scatter suggests "
          "the data cannot pin one down.)")
    for rank, (t, (s, key)) in enumerate(ranked[:top], 1):
        tol = t.translate(TOL)
        print(f"\n#{rank} ngram score={s:.4f}   word score={word_score(tol)}")
        print_lines(t)
        print("   segmented (folded): " + segment(tol))
    s, key = ranked[0][1]
    freq = Counter(cipher_idx)
    print("\nBest key (code -> letter, symbol count):")
    print("   " + "  ".join(f"{uniq[i]}->{LETTERS[key[i]]}({freq[i]})" for i in range(len(uniq))))
    return ranked[0][1][0], ranked[0][0]


# ------------------------------------------------------------------ selftest
def run_selftest(args, pool, workers, timer):
    plain = fold_base(''.join(ch for ch in PLANT.upper() if ch.isalpha()))[:87]
    rng = random.Random(args.seed)
    letters = sorted(set(plain))
    codes = rng.sample(CODES_ALL, len(letters))
    code_of = dict(zip(letters, codes))
    seq = [code_of[ch] for ch in plain]
    uniq, cipher_idx = prep_cipher(seq)
    nsym = len(uniq)
    true_score = SCORE_FN([LIDX[ch] for ch in plain])
    print(f"\n=== SELF-TEST: planted English message, {nsym} distinct symbols, random key ===")
    finals, evals = search_parallel(pool, workers, args.seed, cipher_idx, nsym, args.restarts)
    timer.lap("selftest climb", evals, "score-evals")
    best_score, best_key = max(finals, key=lambda t: t[0])
    text = decrypt_text(cipher_idx, best_key)
    acc = sum(1 for a, b in zip(text, plain) if a == b) / len(plain)
    print(f"true-key score {true_score:.4f}   best found {best_score:.4f}")
    print(f"letters recovered correctly: {acc:.0%}")
    print(f"planted : {plain}")
    print(f"decoded : {text}")
    n_sh = min(args.shuffles, 20)
    sh = run_parallel(pool, shuffle_job,
                      [(args.seed + 500 + i, cipher_idx, nsym, args.restarts) for i in range(n_sh)])
    timer.lap("selftest shuffle baseline", sum(x[2] for x in sh), "score-evals")
    p = pvalue(best_score, [x[0] for x in sh])
    print(f"shuffle p-value for the best score: {p:.4f}  (floor is 1/{n_sh + 1})")
    print("\nHow to read this:")
    if best_score >= true_score - 1e-9:
        print(" * Search adequate: the climb reached at least the true key's score.")
    else:
        print(" * Search WEAK: the true key scores higher than anything found. Raise --restarts.")
    if acc >= 0.9:
        print(" * Objective strong: near-perfect recovery. Real-run results are meaningful.")
    elif acc >= 0.5:
        print(" * Objective partial: partly readable recovery. A real-run hit would be hard to "
              "trust unless the shuffle p-value is small; try --corpus.")
    else:
        print(" * Objective too weak at this length/scorer: a real solution could be missed. "
              "Use --corpus with a large English text file.")


# ------------------------------------------------------------------ main
def main():
    timer = RunTimer("dorabella_hillclimb")
    try:
        ap = argparse.ArgumentParser()
        ap.add_argument('--restarts', type=int, default=100, help='random restarts per search')
        ap.add_argument('--shuffles', type=int, default=30, help='shuffled-ciphertext baselines')
        ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
        ap.add_argument('--top', type=int, default=8)
        ap.add_argument('--seed', type=int, default=1897)
        ap.add_argument('--corpus', help='plain-text English file for quadgram scoring')
        ap.add_argument('--wordlist', help='optional extra word list (one word per line)')
        ap.add_argument('--set', action='append', default=[], metavar='POS=CODE')
        ap.add_argument('--selftest', action='store_true')
        args = ap.parse_args()

        seq = list(BASE_SEQ)
        for item in args.set:
            pos, code = item.split('=')
            assert int(code) in CODES_ALL, f"bad code {code}"
            seq[int(pos) - 1] = int(code)

        build_scorer(args.corpus)
        build_words(args.wordlist)
        workers = max(1, args.workers)
        timer.lap("setup (scorer + word tables)")
        print(f"scoring: {MODE}\nrestarts: {args.restarts}  shuffles: {args.shuffles}  "
              f"workers: {workers}", flush=True)

        pool = None
        if workers > 1:
            pool = mp.Pool(workers, initializer=init_worker,
                           initargs=(args.corpus, args.wordlist))
            timer.lap("worker pool created")

        try:
            if args.selftest:
                run_selftest(args, pool, workers, timer)
                return

            uniq, cipher_idx = prep_cipher(seq)
            nsym = len(uniq)
            print(f"\nCipher: {len(seq)} symbols, {nsym} distinct codes: {uniq}", flush=True)

            finals, evals = search_parallel(pool, workers, args.seed, cipher_idx, nsym,
                                            args.restarts)
            timer.lap("real ciphertext climb", evals, "score-evals")
            real_best, best_text = report_best(finals, cipher_idx, uniq, args.top)
            real_ws = word_score(best_text.translate(TOL))
            timer.lap("report")

            print(f"\nRunning {args.shuffles} shuffled-ciphertext baselines "
                  f"({args.restarts} restarts each)...", flush=True)
            sh = run_parallel(pool, shuffle_job,
                              [(args.seed + 1000 + i, cipher_idx, nsym, args.restarts)
                               for i in range(args.shuffles)])
            timer.lap("shuffle baseline", sum(x[2] for x in sh), "score-evals")

            sh_ng = [x[0] for x in sh]
            sh_ws = [x[1] for x in sh]
            sd = statistics.pstdev(sh_ng) or 1.0
            print("\n=== BASELINE COMPARISON ===")
            print(f"Real best n-gram score:      {real_best:.4f}")
            print(f"Shuffled best n-gram (n={len(sh)}): mean {statistics.mean(sh_ng):.4f}, "
                  f"max {max(sh_ng):.4f}, z of real = {(real_best - statistics.mean(sh_ng)) / sd:.2f}")
            print(f"  -> n-gram p-value: {pvalue(real_best, sh_ng):.4f}")
            print(f"Real word score of best decrypt: {real_ws}   "
                  f"shuffled: mean {statistics.mean(sh_ws):.1f}, max {max(sh_ws)}")
            print(f"  -> word-score p-value (not optimised, so independent): "
                  f"{pvalue(real_ws, sh_ws):.4f}")
            print("\nHow to read this:")
            print(" * Both p-values small (< ~0.02): the real symbol ORDER is more English-like than "
                  "shuffles. Read the top hits by eye before believing anything.")
            print(" * Only n-gram p small: likely overfitting to symbol frequencies, not language.")
            print(" * Both large: nothing beyond chance under this model and scorer.")
            print(" * If --selftest recovered little, a large p-value here is inconclusive, "
                  "not a refutation.")
            timer.lap("statistics")
        finally:
            if pool is not None:
                pool.close()
                pool.join()
    finally:
        timer.finish()


if __name__ == '__main__':
    mp.freeze_support()
    main()