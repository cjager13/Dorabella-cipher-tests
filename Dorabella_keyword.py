#!/usr/bin/env python3
"""
Keyword test: could an 8-letter word set the order AND up/down of the 8 sets of three?

Letters I/J and U/V are interchangeable (as in the cipher alphabet), so a keyword
containing J or V is matched as I or U.

Rule tested: the keyword's letters, in order, sit on directions 1..8 (every starting
direction and both senses are tried too). Each letter must be the FIRST or THIRD letter
of a different set. That letter is the 1-loop letter, so first = set read "up" and
third = set read "down". (Use --three-loop to say the keyword letter is the 3-loop letter.)

Needs Dorabella_hillclimb.py, Dorabella_tests.py, english.txt in the same folder.

  python dora_keyword.py                          # words from english.txt
  python dora_keyword.py --wordlist words.txt     # plus any word list (one word per line)
  python dora_keyword.py --wordlist words.txt --three-loop

Output: screen + keyword_out.txt. Each hit is decoded and scored (real English of this
length scores about -9.2 on the n-gram scale; gibberish about -14 to -18).
"""
import argparse
import itertools
import re
from collections import Counter

import Dorabella_hillclimb as H

GROUPINGS = {'yours': ['ABC', 'DEF', 'GHI', 'KLM', 'NOP', 'QRS', 'TUW', 'XYZ'],
             'standard': ['ABC', 'DEF', 'GHI', 'JKL', 'MNO', 'PQR', 'STU', 'WXY']}
OUT = None


def say(msg=''):
    print(msg, flush=True)
    OUT.write(msg + '\n')


def fold(w):
    return w.replace('J', 'I').replace('V', 'U')


def key_from(sets, combo, start, sense, three_loop):
    """combo[i] = (set, 'first'/'third') for keyword letter i. Returns {code: letter}."""
    key = {}
    for i, (s, end) in enumerate(combo):
        d = (start + sense * i) % 8 + 1
        t = fold(sets[s])
        up = (end == 'first') != three_loop          # which letter sits at 1 loop
        seq = t if up else t[::-1]
        for loop in (1, 2, 3):
            key[loop * 10 + d] = seq[loop - 1]
    return key


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--wordlist', action='append', default=[])
    ap.add_argument('--three-loop', action='store_true')
    ap.add_argument('--show', type=int, default=10)
    args = ap.parse_args()
    OUT = open('keyword_out.txt', 'w', encoding='utf-8')
    H.build_scorer(args.corpus)
    H.build_words(None)

    words = Counter()
    for path in [args.corpus] + args.wordlist:
        with open(path, encoding='utf-8', errors='ignore') as fh:
            words.update(re.findall(r"[A-Z]+", fh.read().upper()))
    words.update(H.EMBEDDED_WORDS)
    cands = {fold(w) for w in words if len(w) == 8}
    say(f"8-letter words to test: {len(cands):,}   keyword letter = "
        f"{'3-loop' if args.three_loop else '1-loop'} letter")

    for g, sets in GROUPINGS.items():
        ends = {}
        for s, t in enumerate(sets):
            t = fold(t)
            ends.setdefault(t[0], []).append((s, 'first'))
            ends.setdefault(t[2], []).append((s, 'third'))
        alt = lambda c: {'I': 'I/J', 'J': 'J/I', 'U': 'U/V', 'V': 'V/U'}.get(c, c)
        say(f"\n=== {g}: allowed letters " + ' '.join('{' + alt(t[0]) + ',' + alt(t[2]) + '}' for t in sets))
        results = []
        for w in sorted(cands):
            if not all(c in ends for c in w):
                continue
            for combo in itertools.product(*[ends[c] for c in w]):
                if len({s for s, _ in combo}) != 8:
                    continue
                for start in range(8):
                    for sense in (1, -1):
                        key = key_from(sets, combo, start, sense, args.three_loop)
                        text = ''.join(key[c] for c in H.BASE_SEQ)
                        sc = H.SCORE_FN([H.LIDX[c] for c in text])
                        results.append((sc, w, start + 1, sense, text))
        keywords = sorted({r[1] for r in results})
        say(f"keywords that fit: {len(keywords)}  {keywords[:40]}")
        say(f"keys decoded: {len(results)}")
        for sc, w, st, se, text in sorted(results, reverse=True)[:args.show]:
            say(f"\n  {w}  start d{st} {'clockwise' if se == 1 else 'anticlockwise'}  score {sc:.3f}  "
                f"words {H.word_score(text.translate(H.TOL))}")
            for a, b in H.LINE_SPANS:
                say("     " + text[a:b])
            say("     " + H.segment(text.translate(H.TOL)))
    OUT.close()
    print("\nsaved to keyword_out.txt")


if __name__ == '__main__':
    main()
