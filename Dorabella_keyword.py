#!/usr/bin/env python3
"""
Keyword test: could a word Dora knew set the order AND up/down of the 8 sets of three?

Rule: the keyword's letters, in order, sit on directions 1..8 (all 8 starting directions
and both senses are tried). Each letter must be the FIRST or THIRD letter of a different
set. It is the 1-loop letter, so first = set read "up", third = "down" (--three-loop flips
that). I/J and U/V are interchangeable, as in the cipher alphabet.

Words are reduced three ways before testing (all are tried and reported):
  as-is    : the word itself                              DORABELLA -> DORABELLA
  doubles  : doubled letters count once                   DORABELLA -> DORABELA
  keyword  : every repeated letter dropped (the classic   DORABELLA -> DORABEL
             Victorian keyword-cipher method)

Match levels:
  FULL  : 8 key letters, each from a different set
  NEAR  : 7 of 8 fit (one letter misfits, or the reduced word has only 7 letters).
          The leftover set goes on the unmatched direction, tried both up and down.

Needs Dorabella_hillclimb.py, Dorabella_tests.py, english.txt in the same folder.

  python Dorabella_keyword.py                              # words from english.txt
  python Dorabella_keyword.py --wordlist words.txt         # add any word list(s)
  python Dorabella_keyword.py --test DORABELLA PENNY       # explain specific words in detail
  python Dorabella_keyword.py --three-loop

Output: screen + keyword_out.txt (log) + keyword_keys.csv (every key tried, scored).
Real English of this length scores about -9.2 on the n-gram scale; gibberish -14 to -18.
"""
import argparse
import csv
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


def alt(c):
    return {'I': 'I/J', 'J': 'J/I', 'U': 'U/V', 'V': 'V/U'}.get(c, c)


def reductions(w):
    w = fold(w)
    doubles = re.sub(r'(.)\1+', r'\1', w)
    seen, kw = set(), []
    for c in w:
        if c not in seen:
            seen.add(c)
            kw.append(c)
    out = {}
    for name, r in (('as-is', w), ('doubles', doubles), ('keyword', ''.join(kw))):
        if len(r) in (7, 8) and r not in out.values():
            out[name] = r
    return out


def end_options(sets):
    ends = {}
    for s, t in enumerate(sets):
        t = fold(t)
        ends.setdefault(t[0], []).append((s, 'first'))
        ends.setdefault(t[2], []).append((s, 'third'))
    return ends


def fits(r, ends):
    """All ways to give the letters of r distinct sets. For 8 letters at most one letter
    may misfit (None); for 7 letters all must fit. Returns (full, near) lists of combos."""
    allow_skip = 1 if len(r) == 8 else 0
    full, near = [], []

    def rec(i, used, combo, skips):
        if i == len(r):
            if len(r) == 8 and skips == 0:
                full.append(list(combo))
            else:
                near.append(list(combo))
            return
        for s, end in ends.get(r[i], []):
            if s not in used:
                used.add(s)
                combo.append((s, end))
                rec(i + 1, used, combo, skips)
                combo.pop()
                used.discard(s)
        if skips < allow_skip:
            combo.append(None)
            rec(i + 1, used, combo, skips + 1)
            combo.pop()

    rec(0, set(), [], 0)
    return full, near


def keys_for(sets, combo, three_loop):
    """Yield (start, sense, fill, key) for a combo (length 7 or 8, may contain one None)."""
    combo = list(combo) + [None] * (8 - len(combo))
    leftover = [s for s in range(8) if s not in {c[0] for c in combo if c}]
    fills = [None]
    if leftover:
        fills = [(leftover[0], 'first'), (leftover[0], 'third')]
    for fill in fills:
        full = [c if c else fill for c in combo]
        for start in range(8):
            for sense in (1, -1):
                key = {}
                for i, (s, end) in enumerate(full):
                    d = (start + sense * i) % 8 + 1
                    t = fold(sets[s])
                    up = (end == 'first') != three_loop
                    seq = t if up else t[::-1]
                    for loop in (1, 2, 3):
                        key[loop * 10 + d] = seq[loop - 1]
                yield start + 1, sense, fill, key


def combo_str(r, combo, sets):
    parts = []
    for i, c in enumerate(r):
        if i < len(combo) and combo[i]:
            s, end = combo[i]
            parts.append(f"{c}={sets[s]}{'^' if end == 'first' else 'v'}")
        else:
            parts.append(f"{c}=(misfit)")
    return ' '.join(parts)


def explain(word, sets, ends):
    say(f"\n  {word}:")
    for name, r in {**{'as-is': fold(word)}, **reductions(word)}.items():
        marks = []
        for c in r:
            opts = ends.get(c, [])
            marks.append(f"{c}:" + ('/'.join(sets[s] + ('^' if e == 'first' else 'v') for s, e in opts)
                                    if opts else 'middle-or-missing'))
        full, near = fits(r, ends) if len(r) in (7, 8) else ([], [])
        verdict = 'FULL fit' if full else 'NEAR fit (7 of 8)' if near else 'no fit'
        say(f"    {name:<8} {r:<12} {verdict:<18} " + '  '.join(marks))


def main():
    global OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--corpus', default='english.txt')
    ap.add_argument('--wordlist', action='append', default=[])
    ap.add_argument('--test', nargs='+', default=['DORABELLA', 'DORA', 'PENNY', 'EDWARD', 'ELGAR'])
    ap.add_argument('--three-loop', action='store_true')
    ap.add_argument('--show', type=int, default=8, help='best decrypts printed per grouping')
    ap.add_argument('--list', type=int, default=150, help='max near-fit words listed per grouping')
    args = ap.parse_args()
    OUT = open('keyword_out.txt', 'w', encoding='utf-8')
    H.build_scorer(args.corpus)
    H.build_words(None)

    words = Counter()
    for path in [args.corpus] + args.wordlist:
        with open(path, encoding='utf-8', errors='ignore') as fh:
            words.update(re.findall(r"[A-Z]+", fh.read().upper()))
    words.update(H.EMBEDDED_WORDS)
    words.update(w.upper() for w in args.test)
    cand = {}
    for w, n in words.items():
        for name, r in reductions(w).items():
            cand.setdefault((r, name), []).append((n, w))
    say(f"words read: {len(words):,}   distinct 7/8-letter key strings after reduction: "
        f"{len({r for r, _ in cand}):,}")
    say(f"keyword letter = {'3-loop' if args.three_loop else '1-loop'} letter;  "
        "^ = set read up, v = set read down")

    csv_f = open('keyword_keys.csv', 'w', newline='', encoding='utf-8')
    cw = csv.writer(csv_f)
    cw.writerow(['grouping', 'level', 'score', 'words_found', 'keyword', 'reduction', 'source_words',
                 'start_dir', 'sense', 'filled_set', 'decrypt'])

    for g, sets in GROUPINGS.items():
        ends = end_options(sets)
        say(f"\n{'=' * 78}\n{g}: allowed letters " +
            ' '.join('{' + alt(fold(t)[0]) + ',' + alt(fold(t)[2]) + '}' for t in sets))
        say('=' * 78)
        say("Your test words (letter: which set it could come from):")
        for w in args.test:
            explain(w.upper(), sets, ends)

        full_hits, near_hits, scored = [], [], []
        for (r, name), src in cand.items():
            full, near = fits(r, ends)
            src_words = ', '.join(w for n, w in sorted(src, reverse=True)[:3])
            freq = sum(n for n, _ in src)
            for level, combos in (('FULL', full), ('NEAR', near)):
                if not combos:
                    continue
                (full_hits if level == 'FULL' else near_hits).append((freq, r, name, src_words, combos[0]))
                for combo in combos:
                    for st, se, fill, key in keys_for(sets, combo, args.three_loop):
                        text = ''.join(key[c] for c in H.BASE_SEQ)
                        sc = H.SCORE_FN([H.LIDX[c] for c in text])
                        wd = H.word_score(text.translate(H.TOL))
                        fs = f"{sets[fill[0]]}{'^' if fill[1] == 'first' else 'v'}" if fill else ''
                        scored.append((sc, wd, level, r, name, src_words, st, se, fs, text))
                        cw.writerow([g, level, f"{sc:.4f}", wd, r, name, src_words, st, se, fs, text])

        full_hits.sort(reverse=True)
        near_hits.sort(reverse=True)
        say(f"\nFULL fits (all 8 letters, one per set): {len(full_hits)}")
        for freq, r, name, src, combo in full_hits:
            say(f"  {r:<9} [{name}] from {src}   {combo_str(r, combo, sets)}")
        say(f"\nNEAR fits (7 of 8), most common words first: {len(near_hits)}")
        for freq, r, name, src, combo in near_hits[:args.list]:
            say(f"  {r:<9} [{name:<7}] from {src:<30} {combo_str(r, combo, sets)}")
        if len(near_hits) > args.list:
            say(f"  ... {len(near_hits) - args.list} more in keyword_keys.csv")

        say(f"\nkeys decoded: {len(scored):,}.  Best by n-gram score:")
        for sc, wd, level, r, name, src, st, se, fs, text in sorted(scored, reverse=True)[:args.show]:
            say(f"\n  {level} {r} [{name}] from {src}; start d{st} "
                f"{'clockwise' if se == 1 else 'anticlockwise'}"
                f"{'; filled ' + fs if fs else ''}   score {sc:.3f}  words {wd}")
            for a, b in H.LINE_SPANS:
                say("     " + text[a:b])
            say("     " + H.segment(text.translate(H.TOL)))
        if scored:
            best_w = max(scored, key=lambda x: (x[1], x[0]))
            say(f"\n  Most words: {best_w[1]} ({best_w[3]} [{best_w[4]}], score {best_w[0]:.3f}): "
                f"{H.segment(best_w[9].translate(H.TOL))}")

    csv_f.close()
    say("\nHow to read: a real key would give a decrypt scoring near -10 that reads as English.")
    say("Decoding many keywords will always turn up a few short words by chance; look for")
    say("whole phrases, not isolated THE/AND/SEE.")
    OUT.close()
    print("\nsaved: keyword_out.txt (log), keyword_keys.csv (every key tried)")


if __name__ == '__main__':
    main()
