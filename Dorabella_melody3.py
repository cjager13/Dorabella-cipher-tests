#!/usr/bin/env python3
"""
Dorabella MELODY test 3 (roadmap items M1 and M2): British tunes, and loops as octaves.

WHY
Dorabella_melody.py read direction as a note of the scale and loop count as note length, scored
against German folk songs (the Essen collection); all three families came out at chance. Two
things could have hidden a tune: the corpus (German songs, not the Irish and Scottish tunes of
Elgar's world), and the reading of the loops (length, not octave). This script changes one at a
time.

Needs Dorabella_melody.py (the version with LENGTH_NAMES) and everything it needs, plus folders of
ABC tunes from the music21 project (https://github.com/cuthbertLab/music21, music21/corpus/):
  M1 needs  oneills1850, ryansMammoth and airdsAirs   (copy the three folders next to the scripts)
  M2 needs  essenFolksong                              (already used by Dorabella_melody.py)
Their licences allow non-commercial use only: do not upload them, or the .txt and .tmp files
built from them, to GitHub.

  python Dorabella_melody3.py --build                 # once: writes melody_british.txt and melody_octave.txt
  python Dorabella_melody3.py BS --selftest           # then, for each family BS, BP, BK, OS, OP, OK:
  python Dorabella_melody3.py BS --stop-when-done     #   results in melody3_BS_out\\ (and so on)
  python Dorabella_melody3.py BS --report-only
  add  --reading ours  (the S and P families only) to use this project's own reading of 5 disputed symbols

PRE-REGISTERED (frozen 19 Nov 2026)
  M1  BRITISH TUNES. Exactly the three families of Dorabella_melody.py (note of the scale x note
      length), scored against about 4,200 Irish and Scottish tunes: O'Neill's Music of Ireland
      (1850 collection), Ryan's Mammoth Collection (1883) and Aird's Selection of Scotch, English,
      Irish and Foreign Airs (1782-1803).
        BS  the scale round the circle, no search (96 keys x 5 reading orders)
        BP  direction is the note, loops are the length, any assignment, no search (241,920 keys)
        BK  any key, by annealing
      Three families: line p < 0.0167.
  M2  LOOPS AS OCTAVE. Direction is the note of the scale (or a rest), as before; the loop count
      is the octave: below, in, or above the tune's commonest octave (a rest counts as "in").
      Scored against the Essen folk songs, as in Dorabella_melody.py, so only the reading of the
      loops changes.
        OS  the scale round the circle, no search (96 keys x 5 reading orders)
        OP  direction is the note, loops are the octave, any assignment, no search (241,920 keys)
        OK  any key, by annealing
      Three families: line p < 0.0167.
  Everything else as in Dorabella_melody.py: the first 85% of each tune file trains the quadgram
  model, planted tunes come from the last 15% (0 or 2 misread symbols), pair-preserving copies are
  the main comparison, plain shuffles are reported, and a pass needs p under the line AND a decoded
  tune a musician would call a tune. Transcription: consensus reading for the S and P families;
  the K families choose the CSV's 9 uncertain symbols inside the search.
"""
import glob
import os
import random
import re
import sys
from collections import Counter

FAMS = ('BS', 'BP', 'BK', 'OS', 'OP', 'OK')
WHICH = next((a for a in sys.argv[1:] if a in FAMS), None) or os.environ.get('DORABELLA_MEL3_FAM')
if WHICH:
    os.environ['DORABELLA_MEL3_FAM'] = WHICH
os.environ['DORABELLA_MEL_FAM'] = 'none'          # stops Dorabella_melody.py setting itself up on import

import itertools                                   # noqa: E402

import Dorabella_melody as MEL                     # noqa: E402

H = MEL.H
BRITISH_FOLDERS = ('oneills1850', 'ryansMammoth', 'airdsAirs')
FILES = {'B': 'melody_british.txt', 'O': 'melody_octave.txt'}
OCTAVES = ('below', 'in', 'above')
OCT_NOTE = ("NOTES: d r m f s l t = do..ti, - = rest; CAPITAL = the octave above the tune's usual one, "
            "with ' = the octave below.")


# ---------------------------------------------------------------- tunes -> tokens
def octave_tokens(key, body):
    """[(pitch 0-7, octave class 0-2)] or None: scale step from the keynote, and below / in / above
    the tune's commonest octave (a rest: pitch 7, 'in')."""
    m = re.match(r'([A-G])', key)
    if not m:
        return None
    tonic = MEL.STEP.index(m.group(1))
    body = re.sub(r'"[^"]*"', '', body)
    body = re.sub(r'![^!]*!', '', body)
    body = re.sub(r'\{[^}]*\}', '', body)
    notes = []
    for _acc, l, octv, _num, _sl, _den in MEL.NOTE_RE.findall(body):
        if l == 'z':
            notes.append((7, None))
        else:
            step = MEL.STEP.index(l.upper()) + (7 if l.islower() else 0) + 7 * octv.count("'") - 7 * octv.count(',')
            notes.append(((step - tonic) % 7, (step - tonic) // 7))
    if len(notes) < 12:
        return None
    usual = Counter(o for p, o in notes if o is not None).most_common(1)[0][0]
    return [(p, 1 if o is None or o == usual else (0 if o < usual else 2)) for p, o in notes]


def show_octave(text):
    out = []
    for ch in text:
        o, p = divmod(H.LIDX[ch], 8)
        s = MEL.SOLFA[p]
        out.append(s.upper() if o == 2 else (s + "'" if o == 0 else s))
    return ' '.join(out)


def build_file(folders, tokens, path):
    files = []
    for d in folders:
        files += sorted(glob.glob(os.path.join(d, '**', '*.abc'), recursive=True))
    if not files:
        raise SystemExit(f"no .abc files found in {', '.join(folders)} (copy the folder(s) next to this script)")
    tunes = []
    for f in files:
        for key, body in MEL.abc_tunes(f):
            tk = tokens(key, body)
            if tk:
                tunes.append(''.join(MEL.token_letter(p, c) for p, c in tk))
    random.Random(1897).shuffle(tunes)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(tunes) + '\n')
    n = sum(len(t) for t in tunes)
    cls = Counter(H.LIDX[ch] // 8 for t in tunes for ch in t)
    print(f"{path}: {len(tunes):,} tunes, {n:,} notes and rests from {len(files)} files in {', '.join(folders)}")
    print("  loop classes (" + ', '.join('short plain long'.split() if tokens is MEL.tune_tokens else OCTAVES)
          + "): " + ' '.join(f"{cls[k] / n:.0%}" for k in range(3)))


def build():
    build_file(BRITISH_FOLDERS, MEL.tune_tokens, FILES['B'])
    build_file(('essenFolksong',), octave_tokens, FILES['O'])
    print("Do not upload these .txt files, the .tmp files or the ABC folders to GitHub.")


# ---------------------------------------------------------------- set up the chosen family
def _fix(text):
    if WHICH and WHICH[0] == 'O':
        for old in ("NOTES: d r m f s l t = do..ti, - = rest; CAPITAL = long, with ' = short.",
                    "NOTES: the letters stand for notes; d r m f s l t = do..ti, - = rest; CAPITAL = long, with ' = short."):
            text = text.replace(old, OCT_NOTE)
    corpus = ("British tunes (O'Neill 1850, Ryan's Mammoth, Aird's Airs)" if WHICH[0] == 'B'
              else 'Essen folk songs, loops read as octaves')
    head, _, rest = text.partition('\n')
    return f"{head}\nCORPUS: {corpus}\n{rest}"


MOD = None
if WHICH in FAMS and os.path.exists(FILES[WHICH[0]]):
    MEL.MELODY_FILE = FILES[WHICH[0]]
    if WHICH[0] == 'O':
        MEL.LENGTH_NAMES = OCTAVES
        MEL.solfa = show_octave
    what = 'length' if WHICH[0] == 'B' else 'octave'
    if WHICH[1] == 'S':
        MOD = MEL.setup_enumerated(WHICH, None, f"{WHICH}: the scale round the circle, loops are the {what}")
    elif WHICH[1] == 'P':
        MOD = MEL.setup_enumerated(WHICH, list(itertools.permutations(range(8))),
                                   f"{WHICH}: direction is the note, loops are the {what}, any assignment")
    else:
        MOD = MEL.setup_anykey(WHICH, f"Melody (note x {what} tokens), any key")
    _old = MOD.build_report

    def _report(*a, **k):
        return _fix(_old(*a, **k))
    MOD.build_report = _report


def main():
    if '--build' in sys.argv:
        build()
        return
    if WHICH is None:
        raise SystemExit("choose a family:  python Dorabella_melody3.py BS | BP | BK | OS | OP | OK  "
                         "[--selftest | --stop-when-done | --report-only]\n(first, once:  python Dorabella_melody3.py --build)")
    if MOD is None:
        raise SystemExit(f"{FILES[WHICH[0]]} not found: run  python Dorabella_melody3.py --build  first")
    argv, skip = [], False
    for a in sys.argv:
        if skip:
            skip = False
        elif a in FAMS or a.startswith('--reading='):
            pass
        elif a == '--reading':
            skip = True
        else:
            argv.append(a)
    sys.argv[:] = argv
    if not any(a == '--out' or a.startswith('--out=') for a in sys.argv):
        ours = MEL.READING == 'ours' and WHICH[1] != 'K'
        sys.argv += ['--out', f"melody3_{WHICH}{'_ours' if ours else ''}_out"]
    if WHICH[1] == 'K':
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'MU']
    elif not any(a == '--corpus' or a.startswith('--corpus=') for a in sys.argv):
        sys.argv += ['--corpus', MEL.write_train_file()]
    MOD.main()


if __name__ == '__main__':
    MEL.Z.mp.freeze_support()
    main()
