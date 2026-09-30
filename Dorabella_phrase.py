#!/usr/bin/env python3
"""
Phrase-key test: could a word or phrase Dora knew set the ORDER and UP/DOWN of the
8 sets of three letters?

Classic Victorian keyword method, applied to the 8 sets:
  ORDER    each letter of the phrase picks the set it belongs to (any position in the
           set). Sets already used are skipped. Sets the phrase never reaches follow in
           normal order. So DORABELLA (standard grouping) gives DEF MNO PQR ABC JKL,
           then GHI STU WXY.
  UP/DOWN  tried five ways:
             letter    : the letter that picked the set decides: first letter = up,
                         third = down, middle letter = both tried
             zigzag-up : up, down, up, down... (your notebook pattern)
             zigzag-dn : down, up, down, up...
             all-up / all-down
  PLACE    set #1 goes on each of the 8 directions, clockwise or anticlockwise (16).

Each phrase is also tried as: vowels dropped (first letter kept), and initials (for
multi-word phrases).

  python dora_phrase.py                          # built-in Elgar/Dora phrases
  python dora_phrase.py --phrases myphrases.txt  # add your own, one per line
  python dora_phrase.py --grouping standard      # or yours / both (default both)

Output: screen + phrase_out.txt (readable log) + phrase_keys.csv (every key, scored).

Scores: n-gram (English-likeness; real English ~ -9.2, gibberish -14 to -18) and word
score. Because phonetic/abbreviated text can score poorly, READ the decrypts by eye as
well: phrase_out.txt prints every phrase's best decrypts, not only the top ones.
'chance p' = how often a random set-order with the same number of tries does as well.
"""
import argparse
import csv
import itertools
import random
import re

import Dorabella_hillclimb as H

GROUPINGS = {'yours': ['ABC', 'DEF', 'GHI', 'KLM', 'NOP', 'QRS', 'TUW', 'XYZ'],
             'standard': ['ABC', 'DEF', 'GHI', 'JKL', 'MNO', 'PQR', 'STU', 'WXY']}
# letters not written in a grouping's sets are mapped to their alphabet twin
TWINS = {'yours': {'J': 'I', 'V': 'U'}, 'standard': {'V': 'U', 'Z': 'S'}}

BUILT_IN = """
DORA
DORABELLA
DORA PENNY
MISS PENNY
MISS DORA PENNY
MIS DORA PENY
MISDORAPENY
DORA BELLA
DORABELLA PENNY
MISS DORABELLA
DORIS
DORIS PENNY
PENNY
ALFRED PENNY
MRS PENNY
WOLVERHAMPTON
THE RECTORY
WOLVERHAMPTON RECTORY
EDWARD ELGAR
ELGAR
EDOO
EDU
E E
ALICE ELGAR
ALICE
CARICE
CRAEG LEA
FORLI
MALVERN
WELLS ROAD
ENIGMA
VARIATION
VARIATIONS
ENIGMA VARIATIONS
DORABELLA VARIATION
INTERMEZZO
NIMROD
JAEGER
SCHUSTER
BRAHMS
WAGNER
MENDELSSOHN
WORCESTER
CATHEDRAL
FOOTBALL
WOLVES
WOLVERHAMPTON WANDERERS
POWICK
HEREFORD
BIRCHWOOD
KEY
CIPHER
CYPHER
SECRET
PUZZLE
RIDDLE
JULY
JULY 1897
FOURTEENTH JULY
MUSIC
DEAR DORA
MY DEAR DORA
DORA DEAR
"""

OUT = None


def say(msg=''):
    print(msg, flush=True)
    OUT.write(msg + '\n')


def letters(phrase):
    return ''.join(c for c in phrase.upper() if 'A' <= c <= 'Z')


def variants(phrase):
    base = letters(phrase)
    out = {'as-is': base}
    nv = base[0] + re.sub('[AEIOU]', '', base[1:]) if base else ''
    if nv and nv != base:
        out['no-vowels'] = nv
    words = [w for w in re.split(r'[^A-Za-z]+', phrase.upper()) if w]
    if len(words) > 1:
        out['initials'] = ''.join(w[0] for w in words)
    return out


def set_of(ch, sets, twins):
    ch = twins.get(ch, ch)
    for s, t in enumerate(sets):
        if ch in t:
            return s, t.index(ch)
    return None, None


def order_and_hints(word, sets, twins):
    """Set order (by first appearance) and, for each set, the position of the letter
    that picked it (0 first, 1 middle, 2 third; None if appended)."""
    order, hint = [], {}
    for ch in word:
        s, pos = set_of(ch, sets, twins)
        if s is not None and s not in hint:
            order.append(s)
            hint[s] = pos
    for s in range(8):
        if s not in hint:
            order.append(s)
            hint[s] = None
    return order, hint


def updown_patterns(order, hint):
    """Yield (mode, tuple of down-flags aligned with order)."""
    opts = []
    for s in order:
        h = hint[s]
        opts.append((False,) if h == 0 else (True,) if h == 2 else (False, True))
    seen = set()
    for combo in itertools.product(*opts):
        if combo not in seen:
            seen.add(combo)
            yield 'letter', combo
    for mode, combo in (('zigzag-up', tuple(i % 2 == 1 for i in range(8))),
                        ('zigzag-dn', tuple(i % 2 == 0 for i in range(8))),
                        ('all-up', (False,) * 8), ('all-down', (True,) * 8)):
        if combo not in seen:
            seen.add(combo)
            yield mode, combo


def build_key(sets, order, downs, start, sense):
    key = {}
    for i, (s, down) in enumerate(zip(order, downs)):
        d = (start + sense * i) % 8 + 1
        t = H.fold_base(sets[s])
        seq = t[::-1] if down else t
        for loop in (1, 2, 3):
            key[loop * 10 + d] = seq[loop - 1]
    return key


def decode(key):
    return ''.join(key[c] for c in H.BASE_SEQ)


def score(text):
    return H.SCORE_FN([H.LIDX[c] for c in text])


def all_keys(sets, order, hint):
    for mode, downs in updown_patterns(order, hint):
        for start in range(8):
            for sense in (1, -1):
                yield mode, downs, start + 1, sense, build_key(sets, order, downs, start, sense)


def best_for(sets, order, hint):
    best = None
    n = 0
    for mode, downs, st, se, key in all_keys(sets, order, hint):
        n += 1
        text = decode(key)
        sc = score(text)
        if best is None or sc > best[0]:
            best = (sc, mode, downs, st, se, text)
    return best, n


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--phrases', action='append', default=[], help='file, one phrase per line')
    ap.add_argument('--no-builtin', action='store_true')
    ap.add_argument('--grouping', choices=['yours', 'standard', 'both'], default='both')
    ap.add_argument('--baseline', type=int, default=2000, help='random set-orders for the chance p')
    ap.add_argument('--show', type=int, default=3, help='decrypts printed per phrase')
    args = ap.parse_args()
    OUT = open('phrase_out.txt', 'w', encoding='utf-8')
    H.build_scorer(args.corpus)
    H.build_words(None)

    phrases = [] if args.no_builtin else [p.strip() for p in BUILT_IN.strip().splitlines() if p.strip()]
    for path in args.phrases:
        with open(path, encoding='utf-8', errors='ignore') as fh:
            phrases += [l.strip() for l in fh if l.strip()]
    say(f"phrases: {len(phrases)}   (each tried as-is, without vowels, and as initials)")

    cw = csv.writer(open('phrase_keys.csv', 'w', newline='', encoding='utf-8'))
    cw.writerow(['grouping', 'phrase', 'variant', 'key_word', 'set_order', 'mode', 'updown',
                 'start_dir', 'sense', 'ngram', 'words', 'decrypt'])

    for g in (['yours', 'standard'] if args.grouping == 'both' else [args.grouping]):
        sets, twins = GROUPINGS[g], TWINS[g]
        say(f"\n{'=' * 78}\nGROUPING {g}: {' '.join(sets)}\n{'=' * 78}")

        # chance baseline: random set orders, random 'letter' flags, same modes & placements
        rng = random.Random(1897)
        base = []
        for _ in range(args.baseline):
            order = list(range(8))
            rng.shuffle(order)
            hint = {s: rng.choice((0, 1, 2)) for s in order}
            base.append(best_for(sets, order, hint)[0][0])
        base.sort()
        say(f"chance baseline ({args.baseline} random set-orders): best-of-family median "
            f"{base[len(base) // 2]:.3f}, 99th pct {base[int(0.99 * len(base))]:.3f}, max {base[-1]:.3f}")

        summary = []
        done = set()
        for ph in phrases:
            for vname, word in variants(ph).items():
                order, hint = order_and_hints(word, sets, twins)
                sig = (tuple(order), tuple(hint[s] for s in order))
                rows = []
                for mode, downs, st, se, key in all_keys(sets, order, hint):
                    text = decode(key)
                    sc = score(text)
                    wd = H.word_score(text.translate(H.TOL))
                    rows.append((sc, wd, mode, downs, st, se, text))
                    cw.writerow([g, ph, vname, word, ' '.join(sets[s] for s in order), mode,
                                 ''.join('v' if d else '^' for d in downs), st, se, f"{sc:.4f}", wd, text])
                rows.sort(reverse=True)
                best = rows[0]
                p = (1 + sum(1 for b in base if b >= best[0])) / (1 + len(base))
                summary.append((best[0], ph, vname, word, p, best, order))
                if sig in done:
                    continue
                done.add(sig)
                say(f"\n{ph} [{vname}: {word}]  order {' '.join(sets[s] for s in order)}   "
                    f"best {best[0]:.3f}  chance p {p:.3f}")
                for sc, wd, mode, downs, st, se, text in rows[:args.show]:
                    say(f"   {sc:.3f} w{wd:<3} {mode:<9} "
                        f"{''.join('v' if d else '^' for d in downs)} d{st}{'cw' if se == 1 else 'acw'}  "
                        f"{text[:29]} {text[29:60]} {text[60:]}")
                    say(f"        {H.segment(text.translate(H.TOL))}")

        summary.sort(reverse=True)
        say(f"\nTOP 10 PHRASE KEYS ({g}):")
        for sc, ph, vname, word, p, best, order in summary[:10]:
            text = best[6]
            say(f"\n  {sc:.3f}  chance p {p:.3f}  {ph} [{vname}]  {best[2]} "
                f"{''.join('v' if d else '^' for d in best[3])} d{best[4]}{'cw' if best[5] == 1 else 'acw'}")
            for a, b in H.LINE_SPANS:
                say("     " + text[a:b])
            say("     " + H.segment(text.translate(H.TOL)))
        n = len(summary)
        say(f"\n  Note: {n} phrase variants were tested, so about {n / 100:.1f} of them would reach "
            "chance p <= 0.01 by luck alone.")

    say("\nHow to read: a real key gives a decrypt near -10 that reads as English, or one that")
    say("reads as Elgar-style phonetic/abbreviated English even if it scores lower. Scan the")
    say("decrypts, not only the numbers.")
    OUT.close()
    print("\nsaved: phrase_out.txt, phrase_keys.csv")


if __name__ == '__main__':
    main()
