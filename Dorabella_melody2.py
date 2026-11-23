#!/usr/bin/env python3
"""
Dorabella MELODY test 2 (roadmap item H7): direction as a melodic STEP, not a note.

WHY
Dorabella_melody.py read each direction as a note of the scale; all three families came out at
chance. A tune can also be written as its MOVES: up one step, down two, the same note again.
Read that way, two neighbouring symbols facing opposite ways are a step and its return (up one,
down one), the commonest figure in folk song. That could explain the cipher's one firm oddity:
neighbours turned half round three times as often as in enciphered English.

Needs Dorabella_melody.py and everything it needs, and the essenFolksong folder of ABC files
(see Dorabella_melody.py for where to get it; not in this repository).

  python Dorabella_melody2.py --build essenFolksong     # once: writes melody_steps.txt
  python Dorabella_melody2.py IS --selftest             # then, for each family IS, IP, IK:
  python Dorabella_melody2.py IS --stop-when-done       #   results in melody2_IS_out\\  (and _IP_, _IK_)
  python Dorabella_melody2.py IS --report-only
  add  --reading ours  (IS and IP only) to use this project's own reading of 5 disputed symbols

HOW A TUNE BECOMES 24 SYMBOLS (frozen 19 Nov 2026)
  step, 8 values      the move from the previous note, counted in steps of the scale (octaves
                      kept): =  (same note), +1, +2, +3, -1, -2, -3, or L (a leap of 4 steps or
                      more either way, a rest, or the first note of a tune)
  length, 3 values    shorter than, equal to, or longer than the tune's commonest note length
  8 x 3 = 24 tokens, written to melody_steps.txt. The first 85% trains the quadgram model; planted
  tunes come only from the last 15%.

PRE-REGISTERED FAMILIES (frozen 19 Nov 2026; three, so the line is p < 0.05/3 = 0.0167)
  IS  OPPOSITE WAYS ARE OPPOSITE MOVES, no search: each pair of opposite directions holds a move
      and its reverse (+1 and -1, +2 and -2, +3 and -3, = and L), in any arrangement round the
      circle (4! x 2^4 = 384), and the three loop counts are the three lengths in any order (6):
      2,304 keys x 5 reading orders. This is the family in which half-turned neighbours are
      returning steps.
  IP  DIRECTION IS THE MOVE, LOOPS ARE LENGTH, any assignment, no search: 8! x 3! = 241,920
      keys x 5 reading orders (IS is inside IP).
  IK  ANY one-symbol-one-token key, by annealing, exactly as Dorabella_languages.py does.
  Compared with pair-preserving copies (main) and plain shuffles through the identical search;
  power from real tunes from the held-back 15%, written with a random key of the family and a
  random reading order, with 0 or 2 misread symbols. A pass needs p under the line AND a decoded
  tune a musician would call a tune.
  Transcription: as in Dorabella_melody.py (IS and IP: consensus reading; IK: the CSV's 9
  uncertain symbols chosen inside the search).
"""
import itertools
import os
import random
import sys
from collections import Counter

FAMS = ('IS', 'IP', 'IK')
WHICH = next((a for a in sys.argv[1:] if a in FAMS), None) or os.environ.get('DORABELLA_MEL2_FAM')
if WHICH:
    os.environ['DORABELLA_MEL2_FAM'] = WHICH
os.environ['DORABELLA_MEL_FAM'] = 'none'          # stops Dorabella_melody.py setting itself up on import

import Dorabella_melody as MEL                     # noqa: E402

H = MEL.H
STEPS = ['=', '+1', '+2', '+3', 'L', '-1', '-2', '-3']      # class k and k+4 are opposite moves
MELODY_FILE = 'melody_steps.txt'
STEP_NOTE = ("NOTES: the letters stand for moves: = same note, +1 up one step, -2 down two, L leap or rest; "
             "a trailing : = long, ' = short.")


def step_tokens(key, body):
    """[(move class 0-7, length class 0-2)] or None, with octaves kept."""
    import re
    m = re.match(r'([A-G])', key)
    if not m:
        return None
    body = re.sub(r'"[^"]*"', '', body)
    body = re.sub(r'![^!]*!', '', body)
    body = re.sub(r'\{[^}]*\}', '', body)
    notes = []
    for _acc, l, octv, num, sl, den in MEL.NOTE_RE.findall(body):
        d = float(num) if num else 1.0
        if sl:
            d /= float(den) if den else 2.0 ** len(sl)
        if l == 'z':
            notes.append((None, d))
        else:
            step = MEL.STEP.index(l.upper()) + (7 if l.islower() else 0) + 7 * octv.count("'") - 7 * octv.count(',')
            notes.append((step, d))
    if len(notes) < 12:
        return None
    mode = Counter(d for p, d in notes if p is not None).most_common(1)[0][0]
    out, prev = [], None
    for p, d in notes:
        length = 0 if d < mode - 1e-9 else (1 if d <= mode + 1e-9 else 2)
        if p is None or prev is None:
            mv = 4
        else:
            iv = p - prev
            mv = {0: 0, 1: 1, 2: 2, 3: 3, -1: 5, -2: 6, -3: 7}.get(iv, 4)
        out.append((mv, length))
        prev = p                      # a rest keeps the last note, so the next move is from it
    return out


def show_steps(text):
    out = []
    for ch in text:
        length, mv = divmod(H.LIDX[ch], 8)
        out.append(STEPS[mv] + (':' if length == 2 else ("'" if length == 0 else '')))
    return ' '.join(out)


def build(folder):
    import glob
    files = sorted(glob.glob(os.path.join(folder, '**', '*.abc'), recursive=True))
    if not files:
        raise SystemExit(f"no .abc files found under {folder}")
    tunes = []
    for f in files:
        for key, body in MEL.abc_tunes(f):
            tk = step_tokens(key, body)
            if tk:
                tunes.append(''.join(MEL.token_letter(mv, ln) for mv, ln in tk))
    random.Random(1897).shuffle(tunes)
    with open(MELODY_FILE, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(tunes) + '\n')
    n = sum(len(t) for t in tunes)
    mc = Counter(H.LIDX[ch] % 8 for t in tunes for ch in t)
    print(f"{len(tunes):,} tunes, {n:,} moves from {len(files)} files -> {MELODY_FILE}")
    print("  move shares: " + ' '.join(f"{STEPS[k]} {mc[k] / n:.0%}" for k in range(8)))
    ret = sum(1 for t in tunes for a, b in zip(t, t[1:])
              if (H.LIDX[a] % 8 - H.LIDX[b] % 8) % 8 == 4)
    pairs = sum(len(t) - 1 for t in tunes)
    print(f"  neighbouring moves that are opposites (the half-turns of family IS): {ret / pairs:.1%} of pairs "
          f"(the Dorabella: 27 of 86 = 31%)")
    print(f"  example: {show_steps(tunes[0][:24])}")
    print(f"Do not upload {MELODY_FILE} or the ABC files to GitHub.")


def is_maps():
    """Directions 0-3 take one move class each from the 4 axes, with a sign; 4-7 take the opposites."""
    maps = []
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((0, 1), repeat=4):
            m = [0] * 8
            for d in range(4):
                k = perm[d] + 4 * signs[d]
                m[d], m[d + 4] = k, (k + 4) % 8
            maps.append(m)
    return maps


def _fix_notes(text):
    for old in ("NOTES: d r m f s l t = do..ti, - = rest; CAPITAL = long, with ' = short.",
                "NOTES: the letters stand for notes; d r m f s l t = do..ti, - = rest; CAPITAL = long, with ' = short."):
        text = text.replace(old, STEP_NOTE)
    return text.replace('notes: ', 'moves: ')


MOD = None
if WHICH in FAMS and os.path.exists(MELODY_FILE):
    MEL.MELODY_FILE = MELODY_FILE
    MEL.SOLFA = STEPS
    MEL.solfa = show_steps
    if WHICH == 'IS':
        MOD = MEL.setup_enumerated('IS', is_maps(), 'IS: opposite directions are opposite moves')
    elif WHICH == 'IP':
        MOD = MEL.setup_enumerated('IP', list(itertools.permutations(range(8))),
                                   'IP: direction is the move, loops are length, any assignment')
    else:
        MOD = MEL.setup_anykey('IK', 'Melody (move x length tokens), any key')
    _old = MOD.build_report

    def _report(*a, **k):
        return _fix_notes(_old(*a, **k))
    MOD.build_report = _report


def main():
    if '--build' in sys.argv:
        i = sys.argv.index('--build')
        if i + 1 >= len(sys.argv):
            raise SystemExit("give the folder of .abc files:  --build essenFolksong")
        build(sys.argv[i + 1])
        return
    if WHICH is None:
        raise SystemExit("choose a family:  python Dorabella_melody2.py IS | IP | IK  [--selftest | --stop-when-done "
                         "| --report-only]\n(first, once:  python Dorabella_melody2.py --build essenFolksong)")
    if MOD is None:
        raise SystemExit(f"{MELODY_FILE} not found: run  python Dorabella_melody2.py --build <folder of .abc files>  first")
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
        ours = MEL.READING == 'ours' and WHICH != 'IK'
        sys.argv += ['--out', f"melody2_{WHICH}{'_ours' if ours else ''}_out"]
    if WHICH == 'IK':
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'MU']
    elif not any(a == '--corpus' or a.startswith('--corpus=') for a in sys.argv):
        sys.argv += ['--corpus', MEL.write_train_file()]
    MOD.main()


if __name__ == '__main__':
    MEL.Z.mp.freeze_support()
    main()
