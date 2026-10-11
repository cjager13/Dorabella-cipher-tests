#!/usr/bin/env python3
"""
Repair text extracted from an old-style typeset PDF (e.g. the 1662 Book of Common Prayer):

  1. words split before an "s" (the long-s glyph):  "requi site" -> "requisite",
     "ea se" -> "ease", "the se per sons" -> "these persons"
  2. a lost "ct" ligature:  "satisfaion" -> "satisfaction", "expe ed" -> "expected",
     "subjes" -> "subjects"

A dictionary decides every repair: a change is made only when it turns a non-word into a
real word. The dictionary is built from clean reference texts you already have.

  python Dorabella_fixtext.py Bible_Psalter\\BCP_1662.txt Bible_Psalter\\BCP_1662_fixed.txt ^
         --ref Bible_KJV\\Bible_KJV.txt Bible_Douay\\Bible_Douay.txt english.txt

Then keep only the _fixed file in the folder. It prints how many repairs of each kind it
made and some examples, so you can check them.
"""
import argparse
import re
from collections import Counter


def load_vocab(paths):
    v = Counter()
    for p in paths:
        with open(p, encoding='utf-8', errors='ignore') as fh:
            for w in re.findall(r'[A-Za-z]+', fh.read()):
                v[w.lower()] += 1
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('infile')
    ap.add_argument('outfile')
    ap.add_argument('--ref', nargs='+', required=True, help='clean texts used as the dictionary')
    ap.add_argument('--map', nargs=2, action='append', default=[], metavar=('HEXCODE', 'LETTERS'),
                    help='replace a special character, e.g. --map E03A ss, or --map E001 ct/st when it varies (repeatable)')
    ap.add_argument('--min-count', type=int, default=1, help='a word must appear this often in --ref')
    args = ap.parse_args()

    vocab = load_vocab(args.ref)
    try:                                       # add the 12,000-word list from the styles test
        import Dorabella_styles as S
        for w in S.base_dictionary('english.txt'):
            vocab[w.lower()] += 3
    except Exception:
        pass
    words = {w for w, c in vocab.items() if c >= args.min_count}

    def isword(w):
        return w.lower() in words

    with open(args.infile, encoding='utf-8', errors='ignore') as fh:
        text = fh.read()
    # special glyphs (ligatures the PDF reader could not name): report them, map known ones
    PUA_LO, PUA_HI = chr(0xE000), chr(0xF8FF)
    special = Counter(ch for ch in text if (PUA_LO <= ch <= PUA_HI) or ord(ch) > 0x2000)
    if special:
        print("special characters found (code, count, examples with the character shown as [#]):")
        for ch, n in special.most_common(30):
            ctx = []
            for m in re.finditer(re.escape(ch), text):
                ctx.append(text[max(0, m.start() - 12):m.end() + 12].replace('\n', ' ').replace(ch, '[#]'))
                if len(ctx) == 3:
                    break
            print(f"    U+{ord(ch):04X}  x{n:<6} " + '  |  '.join(ctx))
    fixes = Counter()
    examples = {}

    def note(kind, s):
        fixes[kind] += 1
        ex = examples.setdefault(kind, [])
        if len(ex) < 12:
            ex.append(s)

    def best(cands, min_count=1):
        cands = [c for c in cands if isword(c) and vocab[c.lower()] >= min_count]
        return max(cands, key=lambda c: vocab[c.lower()]) if cands else None

    # 0: mapped glyphs. "--map E001 ct/st" means the glyph is ct in some words and st in others:
    #    each word gets whichever makes a real word (the first choice if neither does)
    for code, letters in args.map:
        ch = chr(int(code, 16))
        alts = letters.split('/')
        pat = re.compile(r'([A-Za-z]*)( ?)' + re.escape(ch) + r'([A-Za-z]*)')

        def rep(m):
            pre, sp, post = m.groups()
            c = best(pre + L + post for L in alts)
            if c:
                if len(alts) > 1:
                    note(f'glyph U+{code.upper()} chosen per word', f"{pre}[#]{post} -> {c}")
                return c
            if sp and isword(pre):                  # glyph starts a new word
                c = best(L + post for L in alts)
                if c:
                    return pre + ' ' + c
            return pre + alts[0] + post
        text = pat.sub(rep, text)
    text = ''.join(ch for ch in text if not (PUA_LO <= ch <= PUA_HI))

    LOST = ('ct', 'st', 'ss')          # letter pairs the PDF printed as one glyph that can vanish

    # 1: join pieces of a split word: "di splea sure", "Univer sity", "afte r", "epi le" (lost st)
    tok = re.split(r'(\s+)', text)     # keep separators
    W = re.compile(r'[A-Za-z]+')
    WT = re.compile(r'([A-Za-z]+)(\W*)')

    def piece(j, last):
        """the word at token j (trailing punctuation only allowed on the last piece)"""
        if j >= len(tok):
            return None
        m = WT.fullmatch(tok[j] or '')
        if not m or (not last and m.group(2)):
            return None
        return m.group(1), m.group(2)

    def gap(j):
        return j < len(tok) and tok[j] in (' ', '  ')

    out = []
    i = 0
    while i < len(tok):
        t = tok[i]
        done = False
        a = piece(i, False)
        if a and gap(i + 1):
            b = piece(i + 2, False)
            c = piece(i + 4, True) if b and gap(i + 3) else None
            if c:                                   # three pieces, not all real words, common result
                j3 = a[0] + b[0] + c[0]
                if (isword(j3) and vocab[j3.lower()] >= 3 and not isword(a[0] + b[0])
                        and not (isword(a[0]) and isword(b[0]) and isword(c[0]))):
                    out.append(j3 + c[1])
                    note('split word (3 pieces)', f"{a[0]} {b[0]} {c[0]} -> {j3}")
                    i += 5
                    done = True
            if not done:
                b = piece(i + 2, True)
                if b:
                    joined = a[0] + b[0]
                    both = isword(a[0]) and isword(b[0])
                    if isword(joined) and ((not both and vocab[joined.lower()] >= 2)
                                           or (both and b[0][:1] in 'sS' and vocab[joined.lower()] >= 3)):
                        out.append(joined + b[1])
                        note('split word', f"{a[0]} {b[0]} -> {joined}")
                        i += 3
                        done = True
                    elif not both:
                        cand = best((a[0] + L + b[0] for L in LOST), 2)
                        if cand:
                            out.append(cand + b[1])
                            note('split + lost ct/st/ss', f"{a[0]} {b[0]} -> {cand}")
                            i += 3
                            done = True
        if not done:
            out.append(t)
            i += 1
    text = ''.join(out)

    # 2: a word with a lost "ct" inside or a lost "st"/"ss" at the end: "satisfaion", "chri" -> "christ"
    def fix_lost(m):
        w = m.group(0)
        if isword(w) or len(w) < 3:
            return w
        cands = [w + 'st', w + 'ss']
        if len(w) >= 4:
            cands += [w[:k] + 'ct' + w[k:] for k in range(1, len(w) + 1)]
        cand = best(cands, 5)
        if cand:
            note('lost ct/st/ss', f"{w} -> {cand}")
            return cand
        return w
    text = W.sub(fix_lost, text)

    with open(args.outfile, 'w', encoding='utf-8') as fh:
        fh.write(text)
    allw = re.findall(r'[A-Za-z]+', text)
    unknown = Counter(w.lower() for w in allw if not isword(w))
    print(f"dictionary: {len(words):,} words from {len(args.ref)} reference text(s)")
    for k, v in fixes.items():
        print(f"repairs ({k}): {v:,}")
        for e in examples[k]:
            print(f"    {e}")
    print(f"words still not in the dictionary: {sum(unknown.values()):,} of {len(allw):,} "
          f"({sum(unknown.values()) / max(1, len(allw)):.1%}); most common:")
    print("    " + ', '.join(f"{w} ({c})" for w, c in unknown.most_common(25)))
    print(f"-> {args.outfile}")


if __name__ == '__main__':
    main()
