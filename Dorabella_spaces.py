#!/usr/bin/env python3
"""
Dorabella WORD-DIVIDER test, and the main English test again with position 32 reopened.

WHY (audit of 16 Nov 2026)
  1. Every text test so far assumed the message runs its words together, with all 24 symbols
     standing for letters. A common alternative in simple ciphers is that ONE symbol is the gap
     between words. With a gap symbol in the message, a model of run-together English is the
     wrong model, so a real message could have been missed. This was never tested.
  2. Position 32 was fixed as code 18 on 1 Oct 2026, but every published reading, including the
     five-transcription consensus of Hauer et al. (2025), has code 11 there. The any-key tests
     choose the 9 uncertain symbols inside the search; position 32 was not one of them.

Needs Dorabella_languages.py and everything it needs, Dorabella_transcription.csv, english.txt.

  python Dorabella_spaces.py --selftest
  python Dorabella_spaces.py --stop-when-done        # results in spaces_out\\
  python Dorabella_spaces.py --report-only

PRE-REGISTERED (frozen 16 Nov 2026, before any run)
  WS  ONE SYMBOL IS THE WORD GAP. English (Dracula) with its word gaps kept: each gap is
      written as one extra letter (the letter Z is used for it, and the few real Z's are
      written S, so the alphabet stays at 24). The search is ANY one-symbol-one-letter key
      over that alphabet, so it also finds WHICH symbol is the gap. Everything else is
      exactly Dorabella_languages.py: quadgram model from the first 85% of the book, planted
      controls from the last 15%, 5 reading orders, real = median of 9 seeds, 1,000
      pair-preserving copies and 1,000 plain shuffles. One new test, so the line is p < 0.05,
      and a pass also needs a decrypt that reads as English words with gaps.
  EN  REFERENCE, not a new test: ordinary run-together English with the same settings, to show
      that reopening position 32 changes nothing. Expected at chance.
  Both use 10 uncertain symbols: the CSV's 9, plus position 32 (18 or 11).
"""
import os
import sys


import Dorabella_languages as LG                                # noqa: E402

S, H = LG.S, LG.H
GAP = 'Z'
LG.LANGS = [('WS', 'english.txt', 'English with a word-gap symbol'),
            ('EN', 'english.txt', 'English, run together (reference)')]
LG.NAME = {c: n for c, _, n in LG.LANGS}
LG.FILE = {c: f for c, f, _ in LG.LANGS}
LG.NEW = ['WS']
LG.BONF = 0.05
_old_load_language = LG.load_language
_old_load_cipher = S.load_cipher
_old_report = LG.build_report


def load_cipher(path='Dorabella_transcription.csv'):
    seq, alts = _old_load_cipher(path)
    if seq[31] == 18:
        alts = dict(alts)
        alts.setdefault(31, 11)                                  # position 32: 18, or 11 as published
    return seq, alts


def load_language(lang):
    if lang != 'WS':
        return _old_load_language(lang)
    with open(LG.FILE[lang], encoding='utf-8', errors='ignore') as fh:
        words = LG.normalise(fh.read(), 'EN')
    letters = ''.join(w.replace(GAP, 'S') + GAP for w in words)
    cut = int(len(letters) * (1 - LG.HOLDOUT))
    return _GapList(words), letters[:cut], letters[cut:]


class _GapList(list):
    """The book's words, marked so the report knows this is the word-gap language."""


class GapWords(LG.WordList):
    """For WS the 'words' line shows the decrypt with the gap letter turned back into spaces."""

    def __init__(self, words):
        self.gaps = isinstance(words, _GapList)
        if not self.gaps:
            super().__init__(words)

    def segment(self, t):
        return t.replace(GAP, ' ') if self.gaps else super().segment(t)


def build_report(*a, **k):
    t = _old_report(*a, **k)
    t = t.replace('DORABELLA LANGUAGES REPORT', 'DORABELLA WORD-GAP REPORT  (one symbol is the gap between words)')
    t = t.replace(f"Bonferroni line for one of the 4 new languages: p < {LG.BONF:.4f}",
                  "Line to beat: p < 0.0500 for WS (one new test). EN is the run-together reference, not a new test.\n"
                  "In WS the letter Z stands for the gap between words; the 'words:' line shows it as a space.\n"
                  "Uncertain symbols: the CSV's 9 plus position 32 (18 or 11).")
    return t


S.load_cipher = load_cipher
LG.load_language = load_language
LG.WordList = GapWords
LG.build_report = build_report

if __name__ == '__main__':
    LG.mp.freeze_support()
    if not any(a == '--out' or a.startswith('--out=') for a in sys.argv):
        sys.argv += ['--out', 'spaces_out']
    LG.main()
