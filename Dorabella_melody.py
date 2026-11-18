#!/usr/bin/env python3
"""
Dorabella MELODY test: is the cipher a tune rather than a text?

WHY (audit of 16 Nov 2026)
The best-known alternative to "it is a letter" is "it is music": Elgar was a composer, the
symbols have 8 directions (a scale has 7 notes and the octave, or 7 notes and a rest) and 3 loop
counts (short, plain, long). Hauer et al. (HistoCrypt 2021; arXiv 2509.17950, 2025) tried it with
a key search and no scrambled-copy comparison. Rounds 1-2 here only used scales and note names as
KEYS FOR A TEXT; the cipher itself was never read as notes. This script does that, the same fair
way as every other test: the real cipher against scrambled copies, and real tunes planted to
measure power.

Needs Dorabella_circlekeys.py, Dorabella_languages.py and everything they need, numpy, and a
folder of folk tunes in ABC notation (not included; see below).

  python Dorabella_melody.py --build essenFolksong     # once: reads the .abc files, writes melody.txt
  python Dorabella_melody.py MS --selftest             # then, for each family MS, MP, MK:
  python Dorabella_melody.py MS --stop-when-done       #   results in melody_MS_out\\  (and _MP_, _MK_)
  python Dorabella_melody.py MS --report-only
  add  --reading ours  (MS and MP only) to use this project's own reading of 5 disputed symbols

THE TUNES (not in this repository)
  The Essen folksong collection: about 8,500 mostly German folk songs, the kind of plain sung
  tune the idea needs, in the ABC files distributed with music21:
      https://github.com/cuthbertLab/music21  ->  music21/corpus/essenFolksong
  (Code -> Download ZIP, then copy that one folder next to these scripts.) Its licence allows
  non-commercial use and its legal status is described there as unclear, so neither the ABC
  files nor melody.txt should be uploaded to GitHub.

HOW A TUNE BECOMES 24 SYMBOLS (frozen 16 Nov 2026)
  pitch, 8 values     the note's step of the scale counted from the keynote (do re mi fa sol la
                      ti, octave ignored, sharps and flats ignored), or a rest
  length, 3 values    shorter than, equal to, or longer than the tune's commonest note length
  So every note or rest is one of 8 x 3 = 24 tokens, written to melody.txt as 24 letters. The
  first 85% of melody.txt trains a quadgram model (four notes in a row); planted tunes come
  only from the last 15%.

PRE-REGISTERED FAMILIES (frozen 16 Nov 2026; three, so the line is p < 0.05/3 = 0.0167)
  MS  THE SCALE ROUND THE CIRCLE, no search: the 8 directions are do re mi fa sol la ti rest
      in order round the circle, from any starting direction, either way round (16), and the
      three loop counts are the three lengths in any order (6): 96 keys x 5 reading orders.
  MP  DIRECTION IS PITCH, LOOPS ARE LENGTH, in any assignment, no search: all 8! x 3! =
      241,920 keys x 5 reading orders. (MS is inside MP.)
  MK  ANY one-symbol-one-token key, by annealing, exactly as Dorabella_languages.py does for a
      language (this is the Hauer et al. approach, with the scrambled copies they did not use).
  Compared with pair-preserving copies (main) and plain shuffles through the identical search.
  Power: real tunes from the held-back 15%, written with a random key of the family and reading
  order, with 0 or 2 misread symbols.
  Transcription: MS and MP cannot choose readings inside a search, so they use the consensus
  reading of Hauer et al. 2025 (symbols 23, 32, 34, 85, 86 read 14, 11, 12, 14, 12), as in
  Dorabella_notebook.py; "--reading ours" is a secondary run. MK solves the CSV's 9 uncertain
  symbols inside its search, as the language tests do.
  A pass needs p under the line AND a decoded tune that a musician would call a tune.
"""
import glob
import itertools
import os
import random
import re
import sys
from collections import Counter

FAMS = ('MS', 'MP', 'MK')
WHICH = next((a for a in sys.argv[1:] if a in FAMS), None) or os.environ.get('DORABELLA_MEL_FAM')
READING = None
for _i, _a in enumerate(sys.argv):
    if _a == '--reading' and _i + 1 < len(sys.argv):
        READING = sys.argv[_i + 1]
    elif _a.startswith('--reading='):
        READING = _a.split('=', 1)[1]
READING = READING or os.environ.get('DORABELLA_MEL_READING') or 'consensus'
if READING not in ('consensus', 'ours'):
    raise SystemExit("--reading must be 'consensus' or 'ours'")
if WHICH:
    os.environ['DORABELLA_MEL_FAM'] = WHICH
os.environ['DORABELLA_MEL_READING'] = READING

import numpy as np                                              # noqa: E402
import Dorabella_three as Z                                     # noqa: E402

H = Z.H
MELODY_FILE = 'melody.txt'
HOLDOUT = 0.15
CONSENSUS = {23: 14, 32: 11, 34: 12, 85: 14, 86: 12}
SOLFA = ['d', 'r', 'm', 'f', 's', 'l', 't', '-']               # pitch 0-6 and rest
BONF = 0.05 / 3


# ---------------------------------------------------------------- tunes -> tokens
NOTE_RE = re.compile(r"(\^{1,2}|_{1,2}|=)?([A-Ga-gz])([,']*)(\d*)(/*)(\d*)")
STEP = 'CDEFGAB'


def abc_tunes(path):
    with open(path, encoding='utf-8', errors='ignore') as fh:
        txt = fh.read()
    for block in re.split(r'\n(?=X:)', '\n' + txt):
        block = block.strip()
        if not block.startswith('X:'):
            continue
        key, body = None, []
        for line in block.split('\n'):
            line = line.split('%', 1)[0]
            m = re.match(r'^\s*([A-Za-z]):(.*)$', line)
            if m:
                if m.group(1) == 'K' and key is None:
                    key = m.group(2).strip()
                continue                                         # headers, lyrics (w:), in-tune fields
            if key is not None:
                body.append(line)
        if key:
            yield key, '\n'.join(body)


def tune_tokens(key, body):
    """[(pitch 0-7, length class 0-2)] or None. Pitch = scale step from the keynote, 7 = rest."""
    m = re.match(r'([A-G])', key)
    if not m:
        return None
    tonic = STEP.index(m.group(1))
    body = re.sub(r'"[^"]*"', '', body)                          # chord names, fingerings
    body = re.sub(r'![^!]*!', '', body)                          # decorations
    body = re.sub(r'\{[^}]*\}', '', body)                        # grace notes
    notes = []
    for _acc, l, octv, num, sl, den in NOTE_RE.findall(body):
        d = float(num) if num else 1.0
        if sl:
            d /= float(den) if den else 2.0 ** len(sl)
        if l == 'z':
            notes.append((7, d))
        else:
            step = STEP.index(l.upper()) + (7 if l.islower() else 0) + 7 * octv.count("'") - 7 * octv.count(',')
            notes.append(((step - tonic) % 7, d))
    if len(notes) < 12:
        return None
    mode = Counter(d for p, d in notes if p != 7).most_common(1)[0][0]
    return [(p, 0 if d < mode - 1e-9 else (1 if d <= mode + 1e-9 else 2)) for p, d in notes]


def token_letter(pitch, length):
    return H.LETTERS[length * 8 + pitch]


def solfa(text):
    """Letters back to notes: d r m f s l t, - for a rest; CAPITAL = long, with ' = short."""
    out = []
    for ch in text:
        t = H.LIDX[ch]
        length, p = divmod(t, 8)
        s = SOLFA[p]
        out.append(s.upper() if length == 2 else (s + "'" if length == 0 else s))
    return ' '.join(out)


def build(folder):
    files = sorted(glob.glob(os.path.join(folder, '**', '*.abc'), recursive=True))
    if not files:
        raise SystemExit(f"no .abc files found under {folder}")
    tunes = []
    for f in files:
        for key, body in abc_tunes(f):
            tk = tune_tokens(key, body)
            if tk:
                tunes.append(''.join(token_letter(p, d) for p, d in tk))
    random.Random(1897).shuffle(tunes)                           # so the held-back 15% is a fair sample
    with open(MELODY_FILE, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(tunes) + '\n')
    n = sum(len(t) for t in tunes)
    pc, dc = Counter(), Counter()
    for t in tunes:
        for ch in t:
            length, p = divmod(H.LIDX[ch], 8)
            pc[p] += 1
            dc[length] += 1
    print(f"{len(tunes):,} tunes, {n:,} notes and rests from {len(files)} files -> {MELODY_FILE}")
    print("  pitch shares (do re mi fa sol la ti rest): " + ' '.join(f"{pc[p] / n:.0%}" for p in range(8)))
    print("  length shares (short plain long): " + ' '.join(f"{dc[d] / n:.0%}" for d in range(3)))
    print(f"  example: {solfa(tunes[0][:24])}")
    print("Do not upload melody.txt or the ABC files to GitHub (see the note at the top of this script).")


def held_out_tunes():
    with open(MELODY_FILE, encoding='utf-8') as fh:
        tunes = [t for t in fh.read().split('\n') if t]
    total = sum(len(t) for t in tunes)
    cut, seen, out = total * (1 - HOLDOUT), 0, []
    for t in tunes:
        if seen >= cut:
            out.append(t)
        seen += len(t)
    return out


def write_train_file():
    """The first 85% of melody.txt, as its own file for the no-search families' scorer."""
    with open(MELODY_FILE, encoding='utf-8') as fh:
        letters = ''.join(ch for ch in fh.read() if ch.isalpha())
    path = 'melody_train.tmp'
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(letters[:int(len(letters) * (1 - HOLDOUT))])
    return path


# ---------------------------------------------------------------- no-search families (MS, MP)
def setup_enumerated(fam):
    import Dorabella_circlekeys as C
    M = C.M
    dur_perms = list(itertools.permutations(range(3)))
    if fam == 'MP':
        pitch_maps = np.array(list(itertools.permutations(range(8))), dtype=np.int32)       # [k][direction-1]
    else:
        pitch_maps = np.array([[((d - s) * sense) % 8 for d in range(8)]
                               for s in range(8) for sense in (1, -1)], dtype=np.int32)
    C.NKEYS = len(pitch_maps) * len(dur_perms)
    C.BONF = BONF
    C.DEFAULT_TARGETS = {'CK': dict(shuffle=1000, doublet=2000, controls=100)}
    C.FAMILY_DESC = {'CK': ('MS: the scale round the circle' if fam == 'MS' else
                            'MP: direction is pitch, loops are length, any assignment')
                     + f" ({C.NKEYS:,} keys x 5 reading orders), no search"}

    def all_keys():
        if C.KEYS is None:
            parts = []
            for dp in dur_perms:
                K = np.empty((len(pitch_maps), 24), dtype=np.int32)
                for loop in range(3):
                    K[:, loop * 8:loop * 8 + 8] = dp[loop] * 8 + pitch_maps
                parts.append(K)
            C.KEYS = np.concatenate(parts)
        return C.KEYS

    def key_of(ix):
        return all_keys()[ix].tolist()

    def label(ix):
        di, pi = divmod(ix, len(pitch_maps))
        pm, dp = pitch_maps[pi].tolist(), dur_perms[di]
        dirs = ' '.join(f"{d + 1}={SOLFA[p]}" for d, p in enumerate(pm))
        return (f"directions {dirs} | loops 1,2,3 = "
                + ','.join(('short', 'plain', 'long')[dp[k]] for k in range(3)))

    tunes = held_out_tunes()

    def plant(idx, noise, seed):
        rng = random.Random(seed)
        k = rng.randrange(len(tunes))
        plain = ''
        while len(plain) < 87:                                   # a tune from its start, running on into the next
            plain += tunes[k % len(tunes)]
            k += 1
        plain = plain[:87]
        L = [H.LIDX[ch] for ch in plain]
        ix = rng.randrange(C.NKEYS)
        E = key_of(ix)
        inv = [0] * 24
        for s, x in enumerate(E):
            inv[x] = s
        oi = rng.randrange(len(C.ORDERS))
        R = C.order_perms()[oi]
        seq = [0] * 87
        for i, r in enumerate(R):
            seq[r] = Z.code_of_slot(inv[L[i]])
        for q in rng.sample(range(87), noise):
            loop, d = divmod(seq[q], 10)
            seq[q] = loop * 10 + ((d - 1 + rng.choice((-1, 1))) % 8 + 1)
        return seq, plain, ix, C.ORDERS[oi]

    C.all_keys, C.key_of, C.label, C.plant, C.segment = all_keys, key_of, label, plant, solfa
    old_report = C.build_report

    def build_report(state):
        t = old_report(state)
        for a, b in (("DORABELLA CIRCLE-AND-CROSS KEYS REPORT  (every 'eight sets of three' picture as the key)",
                      f"DORABELLA MELODY REPORT, family {fam}  ({C.FAMILY_DESC['CK']})"),
                     ("circle-and-cross keys", "melody keys"), ("words: ", "notes: "),
                     ("(p < 0.05)", f"(p < {BONF:.4f})")):
            t = t.replace(a, b)
        t = '\n'.join(x for x in t.split('\n') if "Elgar's own key, read forwards" not in x)
        head, _, rest = t.partition('\n')
        note = ("TRANSCRIPTION: consensus reading (Hauer et al. 2025): symbols 23, 32, 34, 85, 86 read 14, 11, 12, 14, 12"
                if READING == 'consensus' else "TRANSCRIPTION: this project's own main readings (secondary run)")
        return (f"{head}\n{note}\nNOTES: d r m f s l t = do..ti, - = rest; CAPITAL = long, with ' = short. "
                f"Line to beat: p < {BONF:.4f} (three melody families).\n{rest}")
    C.build_report = build_report

    def selftest(pool, SEQ, workers):
        import time
        print(f"\n{Z.stamp()} SELF-TEST, family {fam}: {C.NKEYS:,} keys; {len(tunes):,} held-back tunes for controls")
        seq, plain, ix, order = plant(3, 0, Z.derive_seed('selftest'))
        t0 = time.time()
        sm = C.summarize(seq, C.run_family(pool, [seq])[0])
        dt = time.time() - t0
        acc = sum(a == b for a, b in zip(sm['text'], plain)) / 87
        res = []
        for k in range(6):
            s = M.doublet_shuffle(list(SEQ), random.Random(Z.derive_seed('selftest-null', k)))
            res.append(C.summarize(s, C.run_family(pool, [s])[0])['stat'])
        real = C.summarize(SEQ, C.run_family(pool, [SEQ])[0])
        print(f"planted tune: {sm['stat']:.3f}; the real cipher: {real['stat']:.3f}; six pair-preserving copies: "
              + ' '.join(f"{x:.3f}" for x in res))
        print(f"  true key and order found: {'yes' if sm['key'] == ix and sm['order'] == order else 'no'}, "
              f"notes right {acc:.0%}, {H.fmt_secs(dt)} per copy")
        print(f"  planted: {solfa(plain[:28])}")
        print(f"  decoded: {solfa(sm['text'][:28])}")
        print(f"  real cipher, best reading: {solfa(real['text'][:28])}")
        n = sum(C.DEFAULT_TARGETS['CK'].values()) + 1
        print(f"\nEstimated time to reach the targets: {n} copies x {H.fmt_secs(dt)} = {H.fmt_secs(n * dt)}")
    C.selftest = selftest

    old_load = Z.load_seq

    def load_seq(path='Dorabella_transcription.csv'):
        seq, source = old_load(path)
        if READING == 'consensus':
            seq = list(seq)
            for pos, code in CONSENSUS.items():
                seq[pos - 1] = code
            source += ' + consensus reading at 23, 32, 34, 85, 86'
        return seq, source
    Z.load_seq = load_seq
    return C


# ---------------------------------------------------------------- any-key family (MK)
def setup_anykey():
    import Dorabella_languages as LG
    LG.LANGS = [('MU', MELODY_FILE, 'Melody (pitch x length tokens), any key')]
    LG.NAME = {'MU': LG.LANGS[0][2]}
    LG.FILE = {'MU': MELODY_FILE}
    LG.NEW = ['MU']
    LG.BONF = BONF
    class NoteList:                                             # shows each decrypt as notes, not words
        def __init__(self, words):
            pass

        def segment(self, t):
            return solfa(t)
    LG.WordList = NoteList
    old_report = LG.build_report

    def build_report(*a, **k):
        t = old_report(*a, **k).replace('DORABELLA LANGUAGES REPORT', 'DORABELLA MELODY REPORT, family MK  (any key)')
        t = t.replace('      words: ', '      notes: ')
        t = t.replace(f"Bonferroni line for one of the 4 new languages: p < {BONF:.4f}",
                      f"Line to beat: p < {BONF:.4f} (three melody families)")
        head, _, rest = t.partition('\n')
        return (f"{head}\nNOTES: the letters stand for notes; d r m f s l t = do..ti, - = rest; CAPITAL = long, "
                f"with ' = short.\n{rest}")
    LG.build_report = build_report
    return LG


MOD = None
if WHICH in ('MS', 'MP') and os.path.exists(MELODY_FILE):
    MOD = setup_enumerated(WHICH)
elif WHICH == 'MK' and os.path.exists(MELODY_FILE):
    MOD = setup_anykey()


def main():
    if '--build' in sys.argv:
        i = sys.argv.index('--build')
        if i + 1 >= len(sys.argv):
            raise SystemExit("give the folder of .abc files:  --build essenFolksong")
        build(sys.argv[i + 1])
        return
    if WHICH is None:
        raise SystemExit("choose a family:  python Dorabella_melody.py MS | MP | MK  [--selftest | --stop-when-done "
                         "| --report-only]\n(first, once:  python Dorabella_melody.py --build essenFolksong)")
    if MOD is None:
        raise SystemExit(f"{MELODY_FILE} not found: run  python Dorabella_melody.py --build <folder of .abc files>  first")
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
        sys.argv += ['--out', f"melody_{WHICH}{'_ours' if READING == 'ours' and WHICH != 'MK' else ''}_out"]
    if WHICH == 'MK':
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'MU']
    elif not any(a == '--corpus' or a.startswith('--corpus=') for a in sys.argv):
        sys.argv += ['--corpus', write_train_file()]
    MOD.main()


if __name__ == '__main__':
    Z.mp.freeze_support()
    main()
