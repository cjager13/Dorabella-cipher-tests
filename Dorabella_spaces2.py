#!/usr/bin/env python3
"""
Dorabella WORD GAPS IN OTHER LANGUAGES (roadmap G2), and the languages again with symbol 32
reopened (roadmap T1).

Needs Dorabella_spaces.py, Dorabella_languages.py, Dorabella_elgarstyle.py and everything they
need, the language books (french.txt, german.txt, italian.txt, latin.txt, latin_church.txt,
english.txt), and for family WE your private Elgar_letters_to_1897.txt (read locally only; WE is
left out if the file is missing).

  python Dorabella_spaces2.py --selftest
  python Dorabella_spaces2.py --stop-when-done            # all ten families, results in spaces2_out\\
  python Dorabella_spaces2.py --only WF WD --stop-when-done
  python Dorabella_spaces2.py --report-only

Ctrl+C any time; the same command resumes.

WHY
  G2  The word-gap model (one symbol = the gap between words) was only ever tried in English.
  T1  Symbol 32 was fixed as 18 in every run before the November audit; every published reading
      has 11. The English test with symbols 23, 32, 34, 85 and 86 all open has already been done
      (Dorabella_spaces.py, family EN: p = 0.29). The other languages, including the Latin near
      miss (p = 0.049 twice), were never rerun with 32 open. T1 does that.

PRE-REGISTERED (frozen 20 Nov 2026, before any run; ten families, so the line is p < 0.05/10 = 0.005)
  G2  WORD GAPS. For each language, the book with each gap between words written as one extra
      letter. The alphabet must stay at 24 letters, so the language's RAREST letter becomes the
      gap and its few real uses are dropped from the words (the report names the letter).
        WF French   WD German   WI Italian   WL Latin (latin.txt)
        WE Elgar's own voice: Dracula written in Elgar's shorthand (no AND, his abbreviations, as
           style E2 of Dorabella_elgarstyle.py) plus his letters before 14 July 1897 at four times
           their weight, as in style E3; planted text comes from the Dracula part only.
  T1  RUN-TOGETHER TEXT, SYMBOL 32 REOPENED, exactly the round 2 language tests otherwise:
        TF French   TD German   TI Italian   TL Latin (latin.txt)   TC Church Latin (latin_church.txt)
  Every family: ANY one-symbol-one-letter key by annealing over the 5 reading orders, real = median
  of 9 seeds, 1,000 pair-preserving copies (main) and 1,000 plain shuffles, quadgram model from the
  first 85% of the text and planted controls from the last 15% (0 or 2 misread symbols), as in
  Dorabella_languages.py. Uncertain symbols: the CSV's 9 plus position 32 (18 or 11), chosen inside
  the search. A pass also needs a decrypt that reads in that language.
"""
import os
import sys
from collections import Counter

import Dorabella_languages as LG

S, H = LG.S, LG.H
GAP_FAMS = {'WF': ('FR', 'french.txt', 'French with a word-gap symbol'),
            'WD': ('DE', 'german.txt', 'German with a word-gap symbol'),
            'WI': ('IT', 'italian.txt', 'Italian with a word-gap symbol'),
            'WL': ('LA', 'latin.txt', 'Latin with a word-gap symbol'),
            'WE': ('EN', 'english.txt', "Elgar's own voice with a word-gap symbol")}
RERUN_FAMS = {'TF': ('FR', 'french.txt', 'French, symbol 32 reopened'),
              'TD': ('DE', 'german.txt', 'German, symbol 32 reopened'),
              'TI': ('IT', 'italian.txt', 'Italian, symbol 32 reopened'),
              'TL': ('LA', 'latin.txt', 'Latin, symbol 32 reopened'),
              'TC': ('LA', 'latin_church.txt', 'Church Latin, symbol 32 reopened')}
ALL = {**GAP_FAMS, **RERUN_FAMS}
BONF = 0.05 / len(ALL)
GAPS = {}                    # family -> the letter used as the gap

_old_load_cipher = S.load_cipher


def load_cipher(path='Dorabella_transcription.csv'):
    seq, alts = _old_load_cipher(path)
    if seq[31] == 18:
        alts = dict(alts)
        alts.setdefault(31, 11)                      # position 32: 18, or 11 as published
    return seq, alts


class GapList(list):
    gap = None


def _words(fam):
    base, path, _ = ALL[fam]
    with open(path, encoding='utf-8', errors='ignore') as fh:
        return LG.normalise(fh.read(), base)


def _elgar_words():
    """Dracula in Elgar's shorthand (style E2), and his letters, as word lists (letters may be empty)."""
    import Dorabella_elgarstyle as ES
    drac = [ES.e2(w) for w in _words('WE')]
    drac = [w for w in drac if w]
    lt = ES.letters_text('')
    mine = []
    if lt:
        mine = [ES.ABBR.get(w, w) for w in S.words_of_text(lt)]
        mine = [H.fold_base(w) for w in mine if w]
    return drac, mine


def load_language(fam):
    if fam in RERUN_FAMS:
        base, path, _ = RERUN_FAMS[fam]
        with open(path, encoding='utf-8', errors='ignore') as fh:
            words = LG.normalise(fh.read(), base)
        letters = ''.join(words)
        cut = int(len(letters) * (1 - LG.HOLDOUT))
        return words, letters[:cut], letters[cut:]
    if fam == 'WE':
        drac, mine = _elgar_words()
        words = drac
    else:
        words = _words(fam)
        mine = []
    cnt = Counter(''.join(words + mine))
    gap = min(H.LETTERS, key=lambda ch: (cnt.get(ch, 0), ch))
    GAPS[fam] = gap

    def gapped(ws):
        return ''.join(w.replace(gap, '') + gap for w in ws if w.replace(gap, ''))
    text = gapped(words)
    cut = int(len(text) * (1 - LG.HOLDOUT))
    train = text[:cut] + gapped(mine) * 4               # WE: his letters at four times their weight
    out = GapList(words)
    out.gap = gap
    return out, train, text[cut:]


class GapWords(LG.WordList):
    """Display only: for the gap families, the decrypt with the gap letter shown as spaces."""

    def __init__(self, words):
        self.gap = getattr(words, 'gap', None)
        if self.gap is None:
            super().__init__(words)

    def segment(self, t):
        return t.replace(self.gap, ' ') if self.gap else super().segment(t)


_old_report = LG.build_report


def build_report(*a, **k):
    t = _old_report(*a, **k)
    t = t.replace('DORABELLA LANGUAGES REPORT',
                  'DORABELLA WORD GAPS IN OTHER LANGUAGES (G2) AND LANGUAGES WITH SYMBOL 32 REOPENED (T1)')
    t = t.replace(f"Bonferroni line for one of the 4 new languages: p < {BONF:.4f}",
                  f"Line to beat: p < {BONF:.4f} (ten pre-registered families). Uncertain symbols: the CSV's 9 plus "
                  f"position 32 (18 or 11).\nGap letters (shown as spaces on the 'words' lines): " +
                  (', '.join(f"{f} {g}" for f, g in sorted(GAPS.items())) or 'not loaded yet'))
    t = t.replace('(expected by chance among 4 languages)', '(not past the line for ten families)')
    return t


S.load_cipher = load_cipher
LG.LANGS = [(f, path, name) for f, (_, path, name) in ALL.items()]
LG.NAME = {f: name for f, (_, _, name) in ALL.items()}
LG.FILE = {f: path for f, (_, path, _) in ALL.items()}
LG.NEW = list(ALL)
LG.BONF = BONF
LG.load_language = load_language
LG.WordList = GapWords
LG.build_report = build_report

if __name__ == '__main__':
    LG.mp.freeze_support()
    if not os.path.exists('Elgar_letters_to_1897.txt'):
        print("Elgar_letters_to_1897.txt not found: family WE is left out")
        LG.LANGS = [x for x in LG.LANGS if x[0] != 'WE']
    if not any(a == '--out' or a.startswith('--out=') for a in sys.argv):
        sys.argv += ['--out', 'spaces2_out']
    LG.main()
