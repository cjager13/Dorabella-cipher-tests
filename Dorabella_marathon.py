#!/usr/bin/env python3
"""
Dorabella marathon: a long-running, resumable battery of tests.

Needs Dorabella_hillclimb.py, Dorabella_tests.py and english.txt in the same folder.

  python Dorabella_marathon.py                 # start, or resume where it left off
  python Dorabella_marathon.py --hours 72      # stop by itself after 72 hours
  python Dorabella_marathon.py --report-only   # just rewrite the report from saved state
  python Dorabella_marathon.py --fresh         # discard saved state and start over

Stop it any time with Ctrl+C. Progress is saved after every unit of work to
marathon_out/state.json, and marathon_out/report.txt is rewritten each time, so
you can open the report while it runs. Running it again resumes.

WHAT IT TESTS (each has its own null, and every null repeats the IDENTICAL search):

 1. MAIN        Your verify search (32 variants x 5 orders, variant-matched z) on
                real data vs many shuffled copies. Two nulls:
                  shuffle  - symbols in random order (what verify used)
                  doublet  - random order that keeps every adjacent symbol PAIR
                             count exactly as in the Dorabella. It keeps the local
                             structure but has no language. Beating it is a much
                             stronger claim than beating plain shuffles.
                Also tests word scores (the untested second condition), each reading
                order, and each alternate reading (32=18, 34=12, ...).
 2. CONTROLS    Known English, encrypted with a random substitution key and run
                through the identical pipeline. It shows what a real solution
                scores, and whether the method can find one.
 3. HELD-OUT    Learn the key from two lines, then decode the third line with that
                frozen key. A real key generalises to unseen text, while an
                overfitted key does not. Max over all 160 configs, vs the same
                procedure on nulls.
 4. DEEP        Simulated annealing (much deeper than the hill climb) over all 160
                configs, with two key models: one-to-one (as before) and homophonic
                (several symbols may share a letter). Compared with identical SA on shuffles.
"""
import argparse
import datetime
import itertools
import json
import math
import multiprocessing as mp
import os
import random
import signal
import statistics
import sys
import time
from collections import Counter, defaultdict

import Dorabella_hillclimb as H
import Dorabella_tests as T

VERSION = 1
NL = H.NL
POSITIONS = sorted(T.AMBIG)
BITS = list(itertools.product((0, 1), repeat=len(POSITIONS)))
CONFIGS = [(vi, o) for vi in range(len(BITS)) for o in T.ORDERS]
REGISTERED = [   # frozen before this script ran (from your verify output + plain baselines)
    ('each line reversed', (0, 1, 1, 1, 1)),
    ('forward', (0, 1, 1, 1, 1)),
    ('forward', (0, 0, 0, 0, 0)),
    ('fully reversed', (0, 0, 0, 0, 0)),
]

# Period-style English for positive controls (written for this script, not in the corpus).
CONTROL_TEXTS = [
    "My dear friend I was so very glad of your letter and the news that you are all well at home "
    "we shall come over on Thursday",
    "Thank you for the music which arrived this morning quite safely I have played it through twice "
    "and like the second part best of all",
    "We walked up the hill after lunch and the view over the valley was lovely though the wind was "
    "rather cold at the top",
    "Please tell your aunt that the parcel will be sent by the early train and that she need not "
    "trouble to meet it herself",
    "I am afraid we cannot come to tea on Sunday as the doctor says I must keep indoors until this "
    "wretched cold has gone",
    "The choir sang better than I had hoped and the organ was in good order so the whole evening was "
    "a great success after all",
    "It was kind of you to remember my birthday and the little book you sent is just the thing for a "
    "long railway journey",
    "Do come and stay with us in the summer when the garden is at its best and we can sit out under "
    "the trees in the evening",
    "We reached the station with only a minute to spare and had to run for the train which made us "
    "laugh for the rest of the way",
    "I have been working hard at the new piece all week and think it is nearly finished but there is "
    "one passage I cannot get right",
    "Your father was here yesterday and told us all about the concert and how well you played in the "
    "second half of the programme",
    "The weather has been so wet that we have hardly been out at all and the children are growing "
    "restless shut up in the house",
    "When you next write do tell me whether the bicycle has been mended for I should like to borrow "
    "it if you can spare it for a day",
    "We had a merry party at the rectory last night with games and singing and did not get home until "
    "long after midnight",
    "I send you this little puzzle to see whether you are as clever as you pretend to be and I shall "
    "expect an answer very soon",
    "The river was high after the storm and we watched the boats from the bridge for an hour before "
    "walking back along the bank",
    "Mother wishes me to say that she hopes you will stay to dinner when you come on Friday and that "
    "there will be room for your brother",
    "It seems a long time since we last met and I often think of the pleasant days we spent together "
    "by the sea last autumn",
    "I found the old letters in a box in the attic and read them all again with great pleasure though "
    "some of them made me rather sad",
    "The new horse is very quiet and we drove him to the market town and back without any trouble at "
    "all which pleased us greatly",
    "Forgive this hasty note but the post goes in ten minutes and I wanted you to know that we arrived "
    "safely and are quite comfortable",
    "Our friends from the north came to stay for a week and we took them to see the cathedral and the "
    "old bridge over the river",
    "I hope that your cold is better and that you will be able to come to the rehearsal on Monday for "
    "we cannot manage without you",
    "The lamp went out just as I was writing this and I had to finish by candle light so you must "
    "excuse the dreadful hand",
    "We are all looking forward to your visit and have made the spare room ready with fresh flowers "
    "from the garden this morning",
    "He says that the piano must be tuned before the party and that he will see to it himself if the "
    "man does not come by Tuesday",
    "I was very sorry to hear of your trouble and wish there were something I could do to help you "
    "beyond sending my kindest thoughts",
    "The snow came early this year and the lane was quite blocked for two days so we had no letters "
    "and no visitors at all",
    "Tell me what you think of the enclosed and whether you can make anything of it for I confess it "
    "has me quite puzzled",
    "We saw the hounds go by this morning and the children ran out to the gate to watch them pass "
    "along the road to the wood",
    "Many thanks for the photograph which is a very good likeness though I think you look rather more "
    "serious than you really are",
    "I went to the library in the town and found the book you wanted but they would not let me take "
    "it away so I copied the page",
]
LINE_LEN = [b - a for a, b in H.LINE_SPANS]

KEEP_AWAKE_OK = False


# ============================================================== small utilities
def now():
    return datetime.datetime.now()


def stamp():
    return f"[{now():%Y-%m-%d %H:%M:%S}]"


def fmt_secs(s):
    return H.fmt_secs(s)


def derive_seed(*parts):
    return random.Random('|'.join(str(p) for p in parts)).getrandbits(48)


def pval(real, nulls, higher=True):
    if not nulls:
        return None
    if higher:
        k = sum(1 for x in nulls if x >= real - 1e-12)
    else:
        k = sum(1 for x in nulls if x <= real + 1e-12)
    return (1 + k) / (1 + len(nulls))


def keep_awake():
    """Stop Windows sleeping while this runs (does not stop lid-close sleep)."""
    global KEEP_AWAKE_OK
    if sys.platform == 'win32':
        try:
            import ctypes
            ctypes.windll.kernel32.SetThreadExecutionState(0x80000000 | 0x00000001)
            KEEP_AWAKE_OK = True
        except Exception:
            pass


def run_parallel(pool, fn, argsets):
    """Like H.run_parallel, but Ctrl+C-friendly on Windows."""
    if pool is None:
        return [fn(*a) for a in argsets]
    res = pool.starmap_async(fn, argsets, chunksize=1)
    while True:                      # short waits: Windows caps one wait at ~49 days,
        try:                         # and short waits keep Ctrl+C responsive
            return res.get(60)
        except mp.TimeoutError:
            continue


# ============================================================== tokens and nulls
def make_tokens(seq, alts):
    """alts: {1-based position: alternate code}. Token = (code, alt, ambiguity idx or -1)."""
    toks = []
    for i, c in enumerate(seq, 1):
        if i in alts:
            toks.append((c, alts[i], POSITIONS.index(i)))
        else:
            toks.append((c, c, -1))
    return toks


def realize(tokens, bits):
    return [alt if (a >= 0 and bits[a]) else code for code, alt, a in tokens]


def doublet_shuffle(seq, rng):
    """Uniform random sequence with the same first symbol and the same multiset of
    adjacent pairs (Altschul-Erickson / Kandel et al.), via a random Eulerian path.
    Wilson's algorithm picks the random last-exit arborescence."""
    n = len(seq)
    if n < 3:
        return list(seq)
    out_edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        out_edges[a].append(b)
    last = seq[-1]
    in_tree = {last}
    nxt = {}
    for v in list(out_edges):
        if v == last:
            continue
        u = v
        while u not in in_tree:
            k = rng.randrange(len(out_edges[u]))
            nxt[u] = k
            u = out_edges[u][k]
        u = v
        while u not in in_tree:
            in_tree.add(u)
            u = out_edges[u][nxt[u]]
    lists = {}
    for v, es in out_edges.items():
        if v == last:
            l = list(es)
            rng.shuffle(l)
        else:
            k = nxt[v]
            l = es[:k] + es[k + 1:]
            rng.shuffle(l)
            l.append(es[k])
        lists[v] = l
    ptr = defaultdict(int)
    out = [seq[0]]
    u = seq[0]
    for _ in range(n - 1):
        w = lists[u][ptr[u]]
        ptr[u] += 1
        out.append(w)
        u = w
    assert Counter(zip(out, out[1:])) == Counter(zip(seq, seq[1:])), "doublet shuffle broke"
    return out


def null_tokens(kind, toks, seed):
    rng = random.Random(seed)
    if kind == 'shuffle':
        t = list(toks)
        rng.shuffle(t)
        return t
    if kind == 'doublet':
        new = doublet_shuffle([c for c, _, _ in toks], rng)
        amb = defaultdict(list)
        for tk in toks:
            if tk[2] >= 0:
                amb[tk[0]].append(tk)
        out = [(c, c, -1) for c in new]
        for c, lst in amb.items():
            slots = [i for i, x in enumerate(new) if x == c]
            for i, tk in zip(rng.sample(slots, len(lst)), lst):
                out[i] = tk
        return out
    raise ValueError(kind)


def oriented_lines(seq, order):
    l = [list(seq[a:b]) for a, b in H.LINE_SPANS]
    if order == 'forward':
        return l
    if order == 'fully reversed':
        return [x[::-1] for x in l[::-1]]
    if order == 'each line reversed':
        return [x[::-1] for x in l]
    if order == 'boustrophedon A':
        return [l[0], l[1][::-1], l[2]]
    return [l[0][::-1], l[1], l[2][::-1]]


# ============================================================== worker jobs
def init_worker(corpus):
    signal.signal(signal.SIGINT, signal.SIG_IGN)   # the main process handles Ctrl+C
    H.init_worker(corpus, None)


def climb_cfg(seed, seq, restarts):
    uniq, cidx = H.prep_cipher(seq)
    finals, ev = H.run_restarts(seed, cidx, len(uniq), restarts)
    s, key = max(finals, key=lambda x: x[0])
    return s, H.decrypt_text(cidx, key), ev


def climb_score(seed, seq, restarts):
    uniq, cidx = H.prep_cipher(seq)
    finals, ev = H.run_restarts(seed, cidx, len(uniq), restarts)
    return max(s for s, _ in finals), ev


def _assign_unseen(test, mapping, free, unseen):
    """Give codes that never appeared in training the best free letters for the test
    line (one-to-one, greedy swaps). Identical in real and null runs, so it is fair."""
    if not unseen:
        return mapping
    m = dict(mapping)
    free = list(free)
    for c, l in zip(unseen, free):
        m[c] = l
    pool_letters = free[len(unseen):]
    cur = H.SCORE_FN([m[c] for c in test])
    improved = True
    while improved:
        improved = False
        for c in unseen:
            for idx in range(len(pool_letters)):
                old = m[c]
                m[c], pool_letters[idx] = pool_letters[idx], old
                s = H.SCORE_FN([m[x] for x in test])
                if s > cur + 1e-12:
                    cur, improved = s, True
                else:
                    pool_letters[idx], m[c] = m[c], old
            for d in unseen:
                if d == c:
                    continue
                m[c], m[d] = m[d], m[c]
                s = H.SCORE_FN([m[x] for x in test])
                if s > cur + 1e-12:
                    cur, improved = s, True
                else:
                    m[c], m[d] = m[d], m[c]
    return m


def heldout_cfg(seed, lines, restarts):
    """Train the key on the other lines, score the held-out line with it frozen."""
    vals = []
    for f in range(len(lines)):
        train = [c for k, l in enumerate(lines) if k != f for c in l]
        test = lines[f]
        uniq, cidx = H.prep_cipher(train)
        finals, _ = H.run_restarts(seed + 7919 * f, cidx, len(uniq), restarts)
        _, key = max(finals, key=lambda x: x[0])
        mapping = {u: key[i] for i, u in enumerate(uniq)}
        used = set(key)
        free = [l for l in range(NL) if l not in used]
        unseen = sorted(set(test) - set(uniq))
        m = _assign_unseen(test, mapping, free, unseen)
        vals.append(H.SCORE_FN([m[c] for c in test]))
    return statistics.mean(vals), vals


SA_T0, SA_T1 = 10.0, 0.3   # tuned on planted English: 96% recovery vs 88% for the climb


def sa_job(seed, seq, restarts, iters, model):
    """Simulated annealing on total quadgram log-prob, with delta scoring.
    model '1to1': each symbol gets a distinct letter. 'homo': letters may repeat."""
    Q = H.QUAD
    get = Q.get
    F = H.QFLOOR
    uniq, cidx = H.prep_cipher(seq)
    nsym, n = len(uniq), len(cidx)
    pos = [[] for _ in range(nsym)]
    for q, c in enumerate(cidx):
        pos[c].append(q)
    starts = [tuple(sorted({s for q in pos[c] for s in range(max(0, q - 3), min(q, n - 4) + 1)}))
              for c in range(nsym)]
    pair_cache = {}

    def pair_starts(i, j):
        k = (i, j)
        v = pair_cache.get(k)
        if v is None:
            v = tuple(sorted(set(starts[i]) | (set(starts[j]) if j < nsym else set())))
            pair_cache[k] = v
        return v

    def qsum(p, ss):
        t = 0.0
        for s in ss:
            t += get(((p[s] * NL + p[s + 1]) * NL + p[s + 2]) * NL + p[s + 3], F)
        return t

    allst = tuple(range(n - 3))
    rng = random.Random(seed)
    exp = math.exp
    rand = rng.random
    rr = rng.randrange
    finals = []
    evals = 0
    fac = (SA_T1 / SA_T0) ** (1.0 / max(1, iters))
    for _ in range(restarts):
        if model == '1to1':
            key = list(range(NL))
            rng.shuffle(key)
        else:
            key = [rr(NL) for _ in range(NL)]
        p = [key[c] for c in cidx]
        cur = qsum(p, allst)
        best, bestkey = cur, key[:]
        Tm = SA_T0
        for _ in range(iters):
            if model == '1to1' or rand() < 0.5:
                i = rr(nsym)
                j = rr(NL - 1)
                if j >= i:
                    j += 1
                if model == 'homo' and j >= nsym:
                    j = rr(nsym)
                    if j == i:
                        Tm *= fac
                        continue
                ss = pair_starts(i, j) if i < j else pair_starts(j, i)
                old = qsum(p, ss)
                li, lj = key[i], key[j]
                for q in pos[i]:
                    p[q] = lj
                if j < nsym:
                    for q in pos[j]:
                        p[q] = li
                d = qsum(p, ss) - old
                if d >= 0 or rand() < exp(d / Tm):
                    key[i], key[j] = lj, li
                    cur += d
                    if cur > best:
                        best, bestkey = cur, key[:]
                else:
                    for q in pos[i]:
                        p[q] = li
                    if j < nsym:
                        for q in pos[j]:
                            p[q] = lj
            else:
                i = rr(nsym)
                l = rr(NL - 1)
                li = key[i]
                if l >= li:
                    l += 1
                ss = starts[i]
                old = qsum(p, ss)
                for q in pos[i]:
                    p[q] = l
                d = qsum(p, ss) - old
                if d >= 0 or rand() < exp(d / Tm):
                    key[i] = l
                    cur += d
                    if cur > best:
                        best, bestkey = cur, key[:]
                else:
                    for q in pos[i]:
                        p[q] = li
            Tm *= fac
            evals += 1
        # guard: tracked total must equal a full rescore
        assert abs(cur - qsum(p, allst)) < 1e-6, "SA delta scoring drifted"
        finals.append((best / (n - 3), tuple(bestkey[:nsym])))
    s, key = max(finals, key=lambda x: x[0])
    hits = sum(1 for f in finals if f[0] >= s - 1e-9)
    text = H.decrypt_text(cidx, key)
    assert abs(H.SCORE_FN([key[c] for c in cidx]) - s) < 1e-9
    return s, text, hits, evals


# ============================================================== searches
def words_of(text):
    return H.word_score(text.translate(H.TOL))


def full_search(pool, toks, seed, restarts):
    jobs = [(seed + i, T.apply_order(realize(toks, BITS[vi]), o), restarts)
            for i, (vi, o) in enumerate(CONFIGS)]
    res = run_parallel(pool, climb_cfg, jobs)
    return [r[0] for r in res], [r[1] for r in res], sum(r[2] for r in res)


def calibrate(pool, toks, seed, restarts, cal):
    jobs = []
    for vi in range(len(BITS)):
        for k in range(cal):
            t = null_tokens('shuffle', toks, derive_seed(seed, 'caltok', vi, k))
            jobs.append((derive_seed(seed, 'calclimb', vi, k), realize(t, BITS[vi]), restarts))
    res = run_parallel(pool, climb_score, jobs)
    mu, sd, idx = [], [], 0
    for vi in range(len(BITS)):
        sc = [res[idx + k][0] for k in range(cal)]
        idx += cal
        mu.append(statistics.mean(sc))
        sd.append(statistics.pstdev(sc) or 1.0)
    return mu, sd, sum(r[1] for r in res)


def summarize(scores, texts, mu, sd, keep_texts=0):
    zs = [(scores[i] - mu[vi]) / sd[vi] for i, (vi, o) in enumerate(CONFIGS)]
    best = max(range(len(zs)), key=lambda i: zs[i])
    words = [words_of(t) for t in texts]
    om = {o: statistics.mean(z for z, (_, oo) in zip(zs, CONFIGS) if oo == o) for o in T.ORDERS}
    bd = []
    for b in range(len(POSITIONS)):
        on = [z for z, (vi, _) in zip(zs, CONFIGS) if BITS[vi][b]]
        off = [z for z, (vi, _) in zip(zs, CONFIGS) if not BITS[vi][b]]
        bd.append(statistics.mean(on) - statistics.mean(off))
    out = dict(best_z=zs[best], best_cfg=best, best_raw=max(scores),
               words_at_best=words[best], max_words=max(words),
               mean_words=statistics.mean(words), order_mean=om, bit_delta=bd,
               win_bits=list(BITS[CONFIGS[best][0]]))
    if keep_texts:
        order = sorted(range(len(zs)), key=lambda i: -zs[i])[:keep_texts]
        out['top'] = [dict(cfg=i, z=zs[i], raw=scores[i], words=words[i], text=texts[i])
                      for i in order]
    return out


def heldout_search(pool, toks, seed, restarts):
    jobs = [(seed + i, oriented_lines(realize(toks, BITS[vi]), o), restarts)
            for i, (vi, o) in enumerate(CONFIGS)]
    res = run_parallel(pool, heldout_cfg, jobs)
    return [r[0] for r in res]


def deep_search(pool, toks, seed, restarts, iters, model):
    jobs = [(seed + i, T.apply_order(realize(toks, BITS[vi]), o), restarts, iters, model)
            for i, (vi, o) in enumerate(CONFIGS)]
    return run_parallel(pool, sa_job, jobs)


def cfg_label(i):
    vi, o = CONFIGS[i]
    ch = [f"{p}={T.AMBIG[p][1]}" for p, b in zip(POSITIONS, BITS[vi]) if b]
    return f"{o} | {','.join(ch) or 'as transcribed'}"


def cfg_index(order, bits):
    return CONFIGS.index((BITS.index(tuple(bits)), order))


# ============================================================== positive controls
def plant_control(idx, noise, seed):
    rng = random.Random(seed)
    raw = CONTROL_TEXTS[idx % len(CONTROL_TEXTS)]
    plain = H.fold_base(''.join(c for c in raw.upper() if c.isalpha()))
    while len(plain) < 87:
        plain += H.fold_base(''.join(c for c in CONTROL_TEXTS[(idx + 1) % len(CONTROL_TEXTS)]
                                     .upper() if c.isalpha()))
    plain = plain[:87]
    letters = sorted(set(plain))
    code_of = dict(zip(letters, rng.sample(H.CODES_ALL, len(letters))))
    seq = [code_of[c] for c in plain]
    amb = set(POSITIONS)
    for q in rng.sample([i for i in range(87) if i + 1 not in amb], noise):
        loop, d = divmod(seq[q], 10)
        seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    alts = {}
    for pos in POSITIONS:
        loop, d = divmod(seq[pos - 1], 10)
        alts[pos] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
    return plain, make_tokens(seq, alts)


def accuracy(text, plain, order):
    ref = ''.join(T.apply_order(list(plain), order))
    return sum(a == b for a, b in zip(text, ref)) / len(ref)


# ============================================================== state
def blank_state(args):
    return dict(version=VERSION, created=str(now()), elapsed=0.0, sessions=0,
                params=dict(restarts=args.restarts, cal=args.cal, ctrl_cal=args.ctrl_cal,
                            seed=args.seed, sa_restarts=args.sa_restarts, sa_iters=args.sa_iters,
                            heldout_restarts=args.heldout_restarts,
                            ctrl_heldout_nulls=args.ctrl_heldout_nulls),
                targets=dict(main_shuffle=400, main_doublet=400, heldout_shuffle=100,
                             heldout_doublet=100, controls=40, deep_1to1=40, deep_homo=40),
                unit_secs={}, real=None,
                main=dict(shuffle=[], doublet=[]),
                heldout=dict(shuffle=[], doublet=[]),
                deep=dict(real={}, **{'1to1': [], 'homo': []}),
                controls=[], log=[])


def save_state(path, state):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        json.dump(state, fh)
    os.replace(tmp, path)


def write_text(path, text):
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as fh:
        fh.write(text)
    os.replace(tmp, path)


# ============================================================== report
def dist_line(vals, fmt='{:.2f}'):
    if not vals:
        return 'none yet'
    return (f"n={len(vals)}, mean {fmt.format(statistics.mean(vals))}, "
            f"max {fmt.format(max(vals))}")


def p_line(real, nulls, higher=True):
    if not nulls:
        return 'p = (no replicates yet)'
    p = pval(real, nulls, higher)
    floor = 1 / (1 + len(nulls))
    k = round(p * (1 + len(nulls))) - 1
    tag = '  <- at the floor: no null matched it yet' if k == 0 else ''
    return f"p = {p:.4f}  ({k} of {len(nulls)} nulls matched or beat it; floor {floor:.4f}){tag}"


def verdict(p, n):
    if p is None or n < 20:
        return 'too few replicates to say'
    if p <= 0.001:
        return 'STRONG: the real data beats this null clearly'
    if p <= 0.01:
        return 'notable: beats this null'
    if p <= 0.05:
        return 'weak: borderline against this null'
    return 'no evidence against this null'


def build_report(state):
    L = []
    P = state['params']
    real = state['real']
    L.append("DORABELLA MARATHON REPORT")
    L.append(f"written {now():%Y-%m-%d %H:%M:%S}   total compute time {fmt_secs(state['elapsed'])}   "
             f"sessions {state['sessions']}")
    L.append(f"settings: climb restarts {P['restarts']}, calibration {P['cal']} per variant, "
             f"SA {P['sa_restarts']} restarts x {P['sa_iters']:,} steps, seed {P['seed']}")
    counts = {k: n for k, n in unit_counts(state).items()}
    L.append("units done: " + ', '.join(f"{k} {v}/{state['targets'][k]}" for k, v in counts.items()))
    if real is None:
        L.append("\n(real-data baseline not computed yet)")
        return '\n'.join(L) + '\n'

    rm = real['main']
    ms, md = state['main']['shuffle'], state['main']['doublet']

    # ---------------- headline
    L.append("\n" + "=" * 78)
    L.append("HEADLINE: does the real Dorabella beat each null? (p floor = 1/(n+1))")
    L.append("=" * 78)
    rows = []

    def add(name, realv, nulls, higher=True):
        p = pval(realv, nulls, higher)
        rows.append((name, realv, nulls, p))

    add("MAIN best matched z  vs shuffle", rm['best_z'], [r['best_z'] for r in ms])
    add("MAIN best matched z  vs doublet", rm['best_z'], [r['best_z'] for r in md])
    add("MAIN words at best z vs shuffle", rm['words_at_best'], [r['words_at_best'] for r in ms])
    add("MAIN max words (160) vs shuffle", rm['max_words'], [r['max_words'] for r in ms])
    add("MAIN max words (160) vs doublet", rm['max_words'], [r['max_words'] for r in md])
    hs, hd = state['heldout']['shuffle'], state['heldout']['doublet']
    add("HELD-OUT max score   vs shuffle", max(real['heldout']), [max(r) for r in hs])
    add("HELD-OUT max score   vs doublet", max(real['heldout']), [max(r) for r in hd])
    for model in ('1to1', 'homo'):
        dr = state['deep']['real'].get(model)
        if dr:
            add(f"DEEP {model:<5} best raw   vs shuffle", max(dr['scores']),
                [r['best'] for r in state['deep'][model]])
    for name, rv, nulls, p in rows:
        pv = f"{p:.4f}" if p is not None else '  -   '
        L.append(f"  {name:<34} real {rv:>9.3f}   n={len(nulls):<5} p={pv:<7} "
                 f"{verdict(p, len(nulls))}")

    # ---------------- controls headline
    C = state['controls']
    L.append("\nPOSITIVE CONTROLS (known English through the identical pipeline):")
    if C:
        zc = [c['best_z'] for c in C]
        L.append(f"  control best matched z: median {statistics.median(zc):.2f}, "
                 f"min {min(zc):.2f}, max {max(zc):.2f}   (Dorabella: {rm['best_z']:.2f})")
        above = sum(1 for z in zc if z > rm['best_z'])
        L.append(f"  controls scoring above the Dorabella: {above} of {len(C)}")
        for nz in (0, 3):
            sub = [c for c in C if c['noise'] == nz]
            if sub:
                L.append(f"  noise {nz} misread symbols: n={len(sub)}, median z "
                         f"{statistics.median(c['best_z'] for c in sub):.2f}, median recovery "
                         f"(forward, as transcribed) {statistics.median(c['acc_fwd'] for c in sub):.0%}, "
                         f"picked forward order {sum(c['win_order'] == 'forward' for c in sub)}/{len(sub)}, "
                         f"median words at best {statistics.median(c['words_at_best'] for c in sub):.0f}, "
                         f"held-out p at floor {sum(c['heldout_p'] <= 1 / (1 + P['ctrl_heldout_nulls']) + 1e-9 for c in sub)}/{len(sub)}")
        L.append("  -> If controls sit far above the Dorabella and read as English, the Dorabella's "
                 "score looks\n     like non-language structure, not a hidden message under this model.")
    else:
        L.append("  none yet")

    # ---------------- main detail
    L.append("\n" + "-" * 78)
    L.append("1. MAIN SEARCH (your verify search, 160 configs, variant-matched z)")
    L.append("-" * 78)
    L.append(f"Real best matched z {rm['best_z']:.2f} at [{cfg_label(rm['best_cfg'])}], "
             f"words there {rm['words_at_best']}; highest raw score (any config) {rm['best_raw']:.4f}, "
             f"max words (any config) {rm['max_words']}")
    for kind, reps in (('shuffle', ms), ('doublet', md)):
        L.append(f"\n  vs {kind} null: best z {dist_line([r['best_z'] for r in reps])}")
        L.append(f"     {p_line(rm['best_z'], [r['best_z'] for r in reps])}")
        L.append(f"     raw best: {p_line(rm['best_raw'], [r['best_raw'] for r in reps])}")
        L.append(f"     words at best z: {p_line(rm['words_at_best'], [r['words_at_best'] for r in reps])}")
        L.append(f"     max words over 160: {p_line(rm['max_words'], [r['max_words'] for r in reps])}")
        if reps:
            L.append("     reading order (mean z over 32 variants), real vs null:")
            for o in T.ORDERS:
                nv = [r['order_mean'][o] for r in reps]
                L.append(f"       {o:<20} real {rm['order_mean'][o]:+.2f}   null mean "
                         f"{statistics.mean(nv):+.2f}   {p_line(rm['order_mean'][o], nv)}")
            L.append("     alternate readings (mean z with it ON minus OFF), real vs null:")
            for b, pos in enumerate(POSITIONS):
                nv = [r['bit_delta'][b] for r in reps]
                win = sum(r['win_bits'][b] for r in reps) / len(reps)
                L.append(f"       {pos}={T.AMBIG[pos][1]:<3} real {rm['bit_delta'][b]:+.2f}   null mean "
                         f"{statistics.mean(nv):+.2f}   ON in null winners {win:.0%}   "
                         f"{p_line(rm['bit_delta'][b], nv)}")
    L.append("\n  Real top configs:")
    for t in rm['top']:
        L.append(f"\n  [{cfg_label(t['cfg'])}]  z={t['z']:.2f}  raw={t['raw']:.4f}  words={t['words']}")
        for a, b in H.LINE_SPANS:
            L.append("     " + t['text'][a:b])
        L.append("     " + H.segment(t['text'].translate(H.TOL)))

    # ---------------- held-out detail
    L.append("\n" + "-" * 78)
    L.append("3. HELD-OUT (key learned on 2 lines, frozen, scored on the 3rd; mean of 3 folds)")
    L.append("-" * 78)
    rh = real['heldout']
    bi = max(range(len(rh)), key=lambda i: rh[i])
    L.append(f"Real best held-out score {rh[bi]:.4f} at [{cfg_label(bi)}]")
    for kind, reps in (('shuffle', hs), ('doublet', hd)):
        L.append(f"  vs {kind}: max over 160 {dist_line([max(r) for r in reps], '{:.4f}')}")
        L.append(f"     {p_line(rh[bi], [max(r) for r in reps])}")
        if reps:
            L.append("     pre-registered configs (no selection, so these p-values stand alone):")
            for o, bits in REGISTERED:
                i = cfg_index(o, bits)
                L.append(f"       [{cfg_label(i)}] real {rh[i]:.4f}   "
                         f"{p_line(rh[i], [r[i] for r in reps])}")
    if C:
        L.append(f"  Controls: held-out p (forward, as transcribed; {P['ctrl_heldout_nulls']} own "
                 f"shuffles each) median {statistics.median(c['heldout_p'] for c in C):.3f}. "
                 "Near the floor = the test has power.")

    # ---------------- deep detail
    L.append("\n" + "-" * 78)
    L.append("4. DEEP SEARCH (simulated annealing, max raw score over all 160 configs)")
    L.append("-" * 78)
    for model, desc in (('1to1', 'one symbol -> one letter'), ('homo', 'homophonic: letters may repeat')):
        dr = state['deep']['real'].get(model)
        if not dr:
            L.append(f"  {model}: real run not done yet")
            continue
        sc = dr['scores']
        bi = max(range(len(sc)), key=lambda i: sc[i])
        nulls = [r['best'] for r in state['deep'][model]]
        L.append(f"\n  {model} ({desc}): real best {sc[bi]:.4f} at [{cfg_label(bi)}], "
                 f"{dr['hits'][bi]}/{P['sa_restarts']} SA restarts reached it")
        L.append(f"     null (shuffle) best: {dist_line(nulls, '{:.4f}')}")
        L.append(f"     {p_line(sc[bi], nulls)}")
        L.append(f"     climb best for comparison: {rm['best_raw']:.4f} "
                 "(if SA is much higher, the old climb was under-searching)")
        order = sorted(range(len(sc)), key=lambda i: -sc[i])[:3]
        for i in order:
            t = dr['texts'][i]
            L.append(f"     [{cfg_label(i)}] {sc[i]:.4f}  words={words_of(t)}")
            for a, b in H.LINE_SPANS:
                L.append("        " + t[a:b])

    # ---------------- controls detail
    if C:
        L.append("\n" + "-" * 78)
        L.append("2. CONTROLS (most recent 8)")
        L.append("-" * 78)
        for c in C[-8:]:
            L.append(f"  #{c['n']} text {c['text_idx']} noise {c['noise']}: best z {c['best_z']:.2f} "
                     f"[{cfg_label(c['best_cfg'])}], fwd z {c['fwd_z']:.2f}, recovery {c['acc_fwd']:.0%}, "
                     f"words {c['words_at_best']}, held-out p {c['heldout_p']:.3f}")
            L.append(f"     planted: {c['plain']}")
            L.append(f"     decoded: {c['fwd_text']}")

    L.append("\n" + "-" * 78)
    L.append("HOW TO READ THIS")
    L.append("-" * 78)
    L.append(" * A p-value can never go below its floor 1/(n+1). 'At the floor' means no null has\n"
             "   matched the real data yet, so more replicates are needed to see how extreme it is.")
    L.append(" * Beating the SHUFFLE null only shows the symbol order is not random. Beating the\n"
             "   DOUBLET null means structure beyond adjacent-pair habits, which is far more interesting.")
    L.append(" * Word tests and held-out tests are not what the search optimises, so they are\n"
             "   independent checks. A real solution should pass them as well.")
    L.append(" * Compare the Dorabella with the CONTROLS. A real English message under this model\n"
             "   should look like them: high z, a readable decode, held-out p at the floor.")
    L.append(" * DEEP homophonic keys can make gibberish score well. Only its p-value counts, never the raw score.")
    return '\n'.join(L) + '\n'


def unit_counts(state):
    return dict(main_shuffle=len(state['main']['shuffle']),
                main_doublet=len(state['main']['doublet']),
                heldout_shuffle=len(state['heldout']['shuffle']),
                heldout_doublet=len(state['heldout']['doublet']),
                controls=len(state['controls']),
                deep_1to1=len(state['deep']['1to1']),
                deep_homo=len(state['deep']['homo']))


# ============================================================== units of work
def unit_real(pool, state, P):
    """Real-data baselines for main, held-out and deep. Runs once."""
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    t0 = time.time()
    print(f"{stamp()} calibrating matched null ({P['cal']} shuffles x 32 variants)...", flush=True)
    mu, sd, _ = calibrate(pool, toks, derive_seed(P['seed'], 'cal'), P['restarts'], P['cal'])
    print(f"{stamp()}   done in {fmt_secs(time.time() - t0)}", flush=True)
    t1 = time.time()
    scores, texts, _ = full_search(pool, toks, derive_seed(P['seed'], 'real-main'), P['restarts'])
    main = summarize(scores, texts, mu, sd, keep_texts=6)
    print(f"{stamp()} real main search: best matched z {main['best_z']:.2f} "
          f"[{cfg_label(main['best_cfg'])}]  ({fmt_secs(time.time() - t1)})", flush=True)
    t2 = time.time()
    ho = heldout_search(pool, toks, derive_seed(P['seed'], 'real-heldout'), P['heldout_restarts'])
    print(f"{stamp()} real held-out: best {max(ho):.4f}  ({fmt_secs(time.time() - t2)})", flush=True)
    state['real'] = dict(mu=mu, sd=sd, main=main, heldout=ho)
    state['unit_secs']['main'] = time.time() - t1


def unit_real_deep(pool, state, P, model):
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    t0 = time.time()
    res = deep_search(pool, toks, derive_seed(P['seed'], 'real-deep', model),
                      P['sa_restarts'], P['sa_iters'], model)
    state['deep']['real'][model] = dict(scores=[r[0] for r in res], texts=[r[1] for r in res],
                                        hits=[r[2] for r in res])
    best = max(r[0] for r in res)
    print(f"{stamp()} real deep {model}: best raw {best:.4f}  ({fmt_secs(time.time() - t0)})",
          flush=True)


def unit_main(pool, state, P, kind):
    real = state['real']
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    n = len(state['main'][kind])
    t = null_tokens(kind, toks, derive_seed(P['seed'], 'main-tok', kind, n))
    scores, texts, _ = full_search(pool, t, derive_seed(P['seed'], 'main-climb', kind, n), P['restarts'])
    s = summarize(scores, texts, real['mu'], real['sd'])
    state['main'][kind].append(s)
    reps = state['main'][kind]
    return (f"best z {s['best_z']:.2f} (real {real['main']['best_z']:.2f}) -> "
            f"p {pval(real['main']['best_z'], [r['best_z'] for r in reps]):.4f}")


def unit_heldout(pool, state, P, kind):
    real = state['real']
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    n = len(state['heldout'][kind])
    t = null_tokens(kind, toks, derive_seed(P['seed'], 'ho-tok', kind, n))
    ho = heldout_search(pool, t, derive_seed(P['seed'], 'ho-climb', kind, n), P['heldout_restarts'])
    state['heldout'][kind].append([round(x, 6) for x in ho])
    reps = state['heldout'][kind]
    return (f"max {max(ho):.4f} (real {max(real['heldout']):.4f}) -> "
            f"p {pval(max(real['heldout']), [max(r) for r in reps]):.4f}")


def unit_deep(pool, state, P, model):
    dr = state['deep']['real'][model]
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    n = len(state['deep'][model])
    t = null_tokens('shuffle', toks, derive_seed(P['seed'], 'deep-tok', model, n))
    res = deep_search(pool, t, derive_seed(P['seed'], 'deep-sa', model, n),
                      P['sa_restarts'], P['sa_iters'], model)
    best_i = max(range(len(res)), key=lambda i: res[i][0])
    state['deep'][model].append(dict(best=res[best_i][0], words=words_of(res[best_i][1])))
    nulls = [r['best'] for r in state['deep'][model]]
    rb = max(dr['scores'])
    return f"best {res[best_i][0]:.4f} (real {rb:.4f}) -> p {pval(rb, nulls):.4f}"


def unit_control(pool, state, P):
    n = len(state['controls'])
    noise = 0 if n % 2 == 0 else 3
    seed = derive_seed(P['seed'], 'ctrl', n)
    plain, toks = plant_control(n, noise, seed)
    mu, sd, _ = calibrate(pool, toks, derive_seed(seed, 'cal'), P['restarts'], P['ctrl_cal'])
    scores, texts, _ = full_search(pool, toks, derive_seed(seed, 'search'), P['restarts'])
    s = summarize(scores, texts, mu, sd)
    fwd = cfg_index('forward', (0,) * len(POSITIONS))
    fwd_z = (scores[fwd] - mu[0]) / sd[0]
    win_order = CONFIGS[s['best_cfg']][1]
    # held-out power check on the correct config, against the control's own shuffles
    lines = oriented_lines(realize(toks, BITS[0]), 'forward')
    jobs = [(derive_seed(seed, 'ho', 0), lines, P['heldout_restarts'])]
    for k in range(P['ctrl_heldout_nulls']):
        t = null_tokens('shuffle', toks, derive_seed(seed, 'hotok', k))
        jobs.append((derive_seed(seed, 'ho', k + 1),
                     oriented_lines(realize(t, BITS[0]), 'forward'), P['heldout_restarts']))
    ho = [r[0] for r in run_parallel(pool, heldout_cfg, jobs)]
    rec = dict(n=n + 1, text_idx=n % len(CONTROL_TEXTS), noise=noise, plain=plain,
               best_z=s['best_z'], best_cfg=s['best_cfg'], win_order=win_order,
               fwd_z=fwd_z, fwd_text=texts[fwd], acc_fwd=accuracy(texts[fwd], plain, 'forward'),
               acc_best=accuracy(texts[s['best_cfg']], plain, win_order),
               words_at_best=s['words_at_best'], max_words=s['max_words'],
               heldout_real=ho[0], heldout_p=pval(ho[0], ho[1:]))
    state['controls'].append(rec)
    return (f"noise {noise}: best z {rec['best_z']:.2f} [{win_order}], recovery {rec['acc_fwd']:.0%}, "
            f"words {rec['words_at_best']}, held-out p {rec['heldout_p']:.3f}")


# ============================================================== self-checks
def self_checks():
    rng = random.Random(1)
    seq = list(H.BASE_SEQ)
    for _ in range(50):
        doublet_shuffle(seq, rng)            # asserts internally
    for o in T.ORDERS:
        assert sum(oriented_lines(seq, o), []) == T.apply_order(seq, o), o
    toks = make_tokens(H.BASE_SEQ, {p: a for p, (_, a) in T.AMBIG.items()})
    for kind in ('shuffle', 'doublet'):
        t = null_tokens(kind, toks, 5)
        assert sorted(t) == sorted(toks), kind
    plain, ctoks = plant_control(0, 0, 3)
    assert len(plain) == 87 and len(ctoks) == 87
    s, text, hits, ev = sa_job(7, list(H.BASE_SEQ), 2, 3000, '1to1')
    s2, text2, _, _ = sa_job(7, list(H.BASE_SEQ), 2, 3000, 'homo')
    return True


# ============================================================== main loop
def pick_stage(state):
    c = unit_counts(state)
    tg = state['targets']
    if all(c[k] >= tg[k] for k in tg):
        for k in tg:
            tg[k] *= 2
        print(f"{stamp()} all targets reached; doubling them: {tg}", flush=True)
    return min(tg, key=lambda k: (c[k] / tg[k], k))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--out', default='marathon_out')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 1)
    ap.add_argument('--hours', type=float, default=0, help='stop after this many hours (0 = never)')
    ap.add_argument('--restarts', type=int, default=80, help='climb restarts (80 = same as verify)')
    ap.add_argument('--cal', type=int, default=64, help='calibration shuffles per variant (verify used 16)')
    ap.add_argument('--ctrl-cal', type=int, default=12)
    ap.add_argument('--heldout-restarts', type=int, default=40)
    ap.add_argument('--ctrl-heldout-nulls', type=int, default=39)
    ap.add_argument('--sa-restarts', type=int, default=16)
    ap.add_argument('--sa-iters', type=int, default=30000)
    ap.add_argument('--seed', type=int, default=1897)
    ap.add_argument('--fresh', action='store_true', help='discard saved state and start over')
    ap.add_argument('--report-only', action='store_true')
    ap.add_argument('--only', nargs='+', help='restrict to these stages, e.g. main_doublet controls')
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    spath = os.path.join(args.out, 'state.json')
    rpath = os.path.join(args.out, 'report.txt')
    lpath = os.path.join(args.out, 'log.txt')

    if os.path.exists(spath) and not args.fresh:
        with open(spath, encoding='utf-8') as fh:
            state = json.load(fh)
        print(f"{stamp()} resuming from {spath} ({fmt_secs(state['elapsed'])} done so far). "
              "Saved settings are used; command-line search settings are ignored.", flush=True)
    else:
        state = blank_state(args)
    P = state['params']

    if args.report_only:
        if state['real'] is not None:
            H.build_scorer(args.corpus)
            H.build_words(None)
        write_text(rpath, build_report(state))
        print(f"report written to {rpath}")
        return

    if not os.path.exists(args.corpus):
        raise SystemExit(f"corpus not found: {args.corpus}")
    for pos, (a, _) in T.AMBIG.items():
        assert H.BASE_SEQ[pos - 1] == a, f"position {pos} is not {a}"

    H.build_scorer(args.corpus)
    H.build_words(None)
    self_checks()
    keep_awake()
    workers = max(1, args.workers)
    print(f"{stamp()} START dorabella_marathon  workers={workers}  "
          f"keep-awake={'on' if KEEP_AWAKE_OK else 'n/a (not Windows)'}", flush=True)
    print(f"  output folder: {os.path.abspath(args.out)}  (report.txt is rewritten as it goes)",
          flush=True)
    state['sessions'] += 1
    pool = mp.Pool(workers, initializer=init_worker, initargs=(args.corpus,)) if workers > 1 else None
    t_start = time.time()
    t_mark = t_start
    logf = open(lpath, 'a', encoding='utf-8')

    def checkpoint(msg):
        nonlocal t_mark
        t = time.time()
        state['elapsed'] += t - t_mark
        t_mark = t
        line = f"{stamp()} {msg}"
        print(line, flush=True)
        logf.write(line + '\n')
        logf.flush()
        save_state(spath, state)
        write_text(rpath, build_report(state))

    try:
        if state['real'] is None:
            unit_real(pool, state, P)
            checkpoint("real-data baselines done")
        for model in ('1to1', 'homo'):
            if model not in state['deep']['real'] and (not args.only or f'deep_{model}' in args.only):
                unit_real_deep(pool, state, P, model)
                checkpoint(f"real deep {model} done")
        while True:
            if args.hours and time.time() - t_start > args.hours * 3600:
                print(f"{stamp()} --hours limit reached", flush=True)
                break
            if args.only:
                c = unit_counts(state)
                stage = min(args.only, key=lambda k: (c[k] / state['targets'][k], k))
            else:
                stage = pick_stage(state)
            t0 = time.time()
            if stage.startswith('main_'):
                msg = unit_main(pool, state, P, stage[5:])
            elif stage.startswith('heldout_'):
                msg = unit_heldout(pool, state, P, stage[8:])
            elif stage.startswith('deep_'):
                model = stage[5:]
                if model not in state['deep']['real']:
                    unit_real_deep(pool, state, P, model)
                msg = unit_deep(pool, state, P, model)
            else:
                msg = unit_control(pool, state, P)
            dt = time.time() - t0
            state['unit_secs'][stage] = dt
            n = unit_counts(state)[stage]
            checkpoint(f"{stage} #{n} ({fmt_secs(dt)}): {msg}")
    except KeyboardInterrupt:
        print(f"\n{stamp()} stopped by Ctrl+C; the unit in progress was discarded, "
              "everything before it is saved.", flush=True)
        if pool is not None:
            pool.terminate()
            pool = None
    except BaseException:
        if pool is not None:         # crash: stop workers now instead of finishing queued work
            pool.terminate()
            pool = None
        raise
    finally:
        if pool is not None:
            pool.close()
            pool.join()
        state['elapsed'] += time.time() - t_mark
        t_mark = time.time()
        save_state(spath, state)
        write_text(rpath, build_report(state))
        logf.close()
        print(f"{stamp()} saved. Report: {os.path.abspath(rpath)}", flush=True)
        print("Run the same command again to resume.", flush=True)


if __name__ == '__main__':
    mp.freeze_support()
    main()
