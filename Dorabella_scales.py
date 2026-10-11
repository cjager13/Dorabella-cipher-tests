#!/usr/bin/env python3
"""
Dorabella SCALES and NUMBERS test: a musical KEY as the cipher key, and the personal numbers
(birthdays, ages, the line lengths 29-31-27).

Built on Dorabella_date.py (same families, same "random numbers instead" comparison, same
controls); needs it and what it needs (Dorabella_three.py and its helpers),
Dorabella_transcription.csv and english.txt in this folder.

  python Dorabella_scales.py --set scales --selftest
  python Dorabella_scales.py --set scales --stop-when-done      # results in scales_out\
  python Dorabella_scales.py --set numbers --selftest
  python Dorabella_scales.py --set numbers --stop-when-done     # results in numbers_out\
  python Dorabella_scales.py --set notes --selftest
  python Dorabella_scales.py --set notes --stop-when-done       # results in notes_out\

THE IDEA
  Elgar's alphabet is 8 directions x 3 loop counts; an octave is 8 notes, do to do. A composer
  writing for a musical young woman could use a SCALE as the key: turn each symbol by the
  steps of the scale, step its loop count, or shift each letter by the scale's numbers.
  (His 1920 notebook shows him laying his symbols out in turning order and drawing circles
  with strokes pointing in and out - the "turn" and "loop step" settings below are exactly
  that kind of move, each tried forwards and backwards.)

PRE-REGISTERED KEYS (frozen before any run)
  --set scales  (13 number patterns, each also as a column order when its digits allow)
    intervals in semitones   major 2212221, natural minor 2122122, harmonic minor 2122131,
                             melodic minor 2122221, Dorian 2122212, Phrygian 1222122,
                             Lydian 2221221, Mixolydian 2212212
    notes of the scale       major 0,2,4,5,7,9,11,12   natural minor 0,2,3,5,7,8,10,12
                             harmonic minor 0,2,3,5,7,8,11,12
    climbing the octave      degrees 0,1,2,3,4,5,6,7 and back down 7,6,5,4,3,2,1,0
  --set numbers
    birthdays (digit keys)   Elgar 2 June 1857: 261857 020657 6257 2657
                             Dora 8 Feb 1874:   821874 080274 2874 8274
                             Alice 9 Oct 1848:  9101848 091048
                             Carice 14 Aug 1890: 1481890 140890
    ages and line lengths    40,23 and 23,40 (Elgar and Dora in July 1897); 29,31,27 and 27,31,29
  --set notes  (added 9 Oct 2026, before any notes run: the NOTE NAMES as the key, CDEFGAB)
    The musical alphabet A B C D E F G = 0 1 2 3 4 5 6, wrapping G -> A. Sharps and flats drop
    out of the letter names, so every major or minor scale written as letters is one of 7
    rising orders (C major = C D E F G A B = 2345601) or, played downwards, one of 7 falling
    orders (C down = C B A G F E D = 2106543). 14 keys, each a 7-digit pattern of 0-6 and each
    a column order of width 7. Random-key comparison: random orderings of 0-6 (same kind).

  Each key is tried as a letter shift, a symbol turn (x 45 degrees) and a loop-count step,
  forwards and backwards, running straight through the message or restarting each line;
  digit keys also as a columnar transposition. Families DE (Elgar's own key, no search),
  DA (any key + shifts/turns/steps) and DT (any key + transposition); Bonferroni p < 0.0167.

  A fixed starting note (tonic) does not matter in DA/DT: with any key, adding the same amount
  to every shift or turn is absorbed by the key. Random-key comparison: every pattern replaced
  by random numbers of the same length and size (scale steps 0-12, digits 0-9, ages 1-99).
"""
import random
import sys

import Dorabella_date as D

SETS = {
    'scales': dict(
        digit=[('2212221', 'major scale steps'), ('2122122', 'natural minor steps'),
               ('2122131', 'harmonic minor steps'), ('2122221', 'melodic minor steps'),
               ('2122212', 'Dorian steps'), ('1222122', 'Phrygian steps'),
               ('2221221', 'Lydian steps'), ('2212212', 'Mixolydian steps')],
        number=[((0, 2, 4, 5, 7, 9, 11, 12), 'major scale notes'),
                ((0, 2, 3, 5, 7, 8, 10, 12), 'natural minor notes'),
                ((0, 2, 3, 5, 7, 8, 11, 12), 'harmonic minor notes'),
                ((0, 1, 2, 3, 4, 5, 6, 7), 'octave climbing'),
                ((7, 6, 5, 4, 3, 2, 1, 0), 'octave descending')],
        rand_hi=12, out='scales_out', title="DORABELLA SCALES REPORT  (a musical key as the cipher key)"),
    'numbers': dict(
        digit=[('261857', 'Elgar D M YYYY'), ('020657', 'Elgar DD MM YY'), ('6257', 'Elgar M D YY'),
               ('2657', 'Elgar D M YY'), ('821874', 'Dora D M YYYY'), ('080274', 'Dora DD MM YY'),
               ('2874', 'Dora M D YY'), ('8274', 'Dora D M YY'), ('9101848', 'Alice D M YYYY'),
               ('091048', 'Alice DD MM YY'), ('1481890', 'Carice D M YYYY'), ('140890', 'Carice DD MM YY')],
        number=[((40, 23), 'ages Elgar,Dora'), ((23, 40), 'ages Dora,Elgar'),
                ((29, 31, 27), 'line lengths'), ((27, 31, 29), 'line lengths reversed')],
        rand_hi=99, out='numbers_out', title="DORABELLA NUMBERS REPORT  (birthdays, ages, line lengths)"),
    'notes': dict(
        digit=[(''.join(str((t + i) % 7) for i in range(7)), f"{'ABCDEFG'[t]} notes up") for t in (2, 3, 4, 5, 6, 0, 1)]
             + [(''.join(str((t - i) % 7) for i in range(7)), f"{'ABCDEFG'[t]} notes down") for t in (2, 3, 4, 5, 6, 0, 1)],
        number=[], rand_hi=0, digit_perm=7,
        out='notes_out', title="DORABELLA NOTES REPORT  (note names as the key, CDEFGAB)"),
}


def pick_set():
    for i, a in enumerate(sys.argv):
        if a == '--set' and i + 1 < len(sys.argv):
            name = sys.argv[i + 1]
            del sys.argv[i:i + 2]
            return name
        if a.startswith('--set='):
            del sys.argv[i]
            return a.split('=', 1)[1]
    raise SystemExit("choose a key set:  --set scales   or   --set numbers   or   --set notes")


def configure(name):
    cfg = SETS[name]
    D.DIGIT_KEYS = cfg['digit']
    D.NUMBER_KEYS = cfg['number']
    D.DATE_MEMBERS = D.build_members(D.DIGIT_KEYS, D.NUMBER_KEYS)
    hi = cfg['rand_hi']
    perm = cfg.get('digit_perm')

    def rand_digits(rng, s):
        if perm:                       # a random ordering of 0..perm-1 (like a note-name order)
            d = list(range(perm))
            rng.shuffle(d)
            return ''.join(map(str, d))
        return ''.join(str(rng.randrange(10)) for _ in s)

    def random_members(seed):
        rng = random.Random(seed)
        dk = [(rand_digits(rng, s), lab) for s, lab in D.DIGIT_KEYS]
        nk = [(tuple(rng.randrange(0, hi + 1) for _ in t), lab) for t, lab in D.NUMBER_KEYS]
        return D.build_members(dk, nk)
    D.random_members = random_members
    old_report = D.build_report

    def build_report(state):
        t = old_report(state)
        return t.replace("DORABELLA DATE-KEY REPORT  ('July 14 97' as the key)", cfg['title'])
    D.build_report = build_report
    D.FAMILY_DESC = {k: v.replace('date-key', f'{name}-key') for k, v in D.FAMILY_DESC.items()}
    if not any(a.startswith('--out') for a in sys.argv):
        sys.argv += ['--out', cfg['out']]
    return cfg


if __name__ == '__main__':
    D.mp.freeze_support()
    name = pick_set()
    if name not in SETS:
        raise SystemExit(f"unknown set {name!r}; choose from: {', '.join(SETS)}")
    cfg = configure(name)
    print(f"key set: {len(D.DIGIT_KEYS)} digit keys, {len(D.NUMBER_KEYS)} number keys -> "
          f"{len(D.DATE_MEMBERS['sub'])} shift/turn/step settings, {len(D.DATE_MEMBERS['col'])} transpositions")
    D.main()
