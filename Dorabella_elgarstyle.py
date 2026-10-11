#!/usr/bin/env python3
"""
Dorabella "ELGAR'S PRIVATE ENGLISH" test: was the message written the way Elgar wrote his
private letters - with "&" instead of "and" (his cipher alphabet has no "&"), and his
habitual abbreviations?

Built on Dorabella_styles.py (same search, same scrambled copies, same controls); needs it and
its helpers, Dorabella_transcription.csv and english.txt in this folder. Style E3 also needs
your private transcription Elgar_letters_to_1897.txt next to the script (it is only read
locally, never copied anywhere).

  python Dorabella_elgarstyle.py --selftest
  python Dorabella_elgarstyle.py --stop-when-done
  python Dorabella_elgarstyle.py --make-phrases     # writes the folder Elgar_letters_quotes for the quote test

THE EVIDENCE (Elgar's own letters, 1879 to July 1897, ~13,300 words):
  "&" 429 times, "and" written out 50 times (about 9 in 10 are "&"); in his one formal
  newspaper article, "and" 36 times and "&" never. He abbreviates constantly: wd. cd. shd.
  yrs. recd. Xmas wh tho' thro' 'tis 'em; Jany. Augt. Septr Novr. Decr.; Wedy Saty.

PRE-REGISTERED STYLES (3 new ones, so the Bonferroni line is p < 0.05/3 = 0.0167)
  E0 standard English (reference only; the same as styles S0)
  E1 NO "AND": every AND left out, as if he wrote "&" and the cipher had no symbol for it
  E2 ELGAR SHORTHAND: no AND, plus his abbreviations (WOULD->WD, COULD->CD, SHOULD->SHD,
     WHICH->WH, YOURS->YRS, RECEIVED->RECD, CHRISTMAS->XMAS, THOUGH->THO, THROUGH->THRO,
     ALTHOUGH->ALTHO, THEM->EM, UNTIL->TILL, months and days as he wrote them)
  E3 ELGAR'S OWN VOICE: the E2 model with Elgar's real letters (before 14 July 1897, his
     newspaper article left out, "&" left out) mixed in at 4 times their weight, so the
     model knows his own words, names and turns of phrase

Each style's quadgram model is trained on Dracula rewritten that way; its dictionary is
rewritten the same way; planted controls are English written that way. Everything else is
exactly Dorabella_styles.py: any one-symbol-one-letter key by annealing, 5 reading orders,
the 9 uncertain symbols chosen inside the search, pair-preserving and plain scrambled copies.
"""
import os
import re
import sys

import Dorabella_styles as S

LETTERS_FILE = 'Elgar_letters_to_1897.txt'
CIPHER_DATE = '1897-07-14'
ABBR = {'WOULD': 'WD', 'COULD': 'CD', 'SHOULD': 'SHD', 'WHICH': 'WH', 'YOURS': 'YRS', 'RECEIVED': 'RECD',
        'CHRISTMAS': 'XMAS', 'THOUGH': 'THO', 'THROUGH': 'THRO', 'ALTHOUGH': 'ALTHO', 'THEM': 'EM',
        'UNTIL': 'TILL', 'JANUARY': 'JANY', 'FEBRUARY': 'FEBY', 'AUGUST': 'AUGT', 'SEPTEMBER': 'SEPTR',
        'OCTOBER': 'OCTR', 'NOVEMBER': 'NOVR', 'DECEMBER': 'DECR', 'WEDNESDAY': 'WEDY', 'SATURDAY': 'SATY',
        'AFFECTIONATELY': 'AFFTLY', 'ACCOUNTS': 'ACCTS', 'SOCIETY': 'SOCY', 'BIRMINGHAM': 'BIRM'}


def e1(w):
    return '' if w == 'AND' else w


def e2(w):
    return '' if w == 'AND' else ABBR.get(w, w)


STYLES = ['E0', 'E1', 'E2', 'E3']
STYLE_NAME = {'E0': 'standard (reference)', 'E1': 'no AND (& unwritten)', 'E2': 'Elgar shorthand',
              'E3': "Elgar's own voice"}

# applied at import, so the worker processes (which import this module) get them too
S.STYLES[:] = STYLES
S.STYLE_NAME.clear()
S.STYLE_NAME.update(STYLE_NAME)
S.STYLE_FN.update({'E0': lambda w: w, 'E1': e1, 'E2': e2, 'E3': e2})
S.BONF = 0.05 / 3
_orig_build_quad = S.build_quad
_orig_init_worker = S.init_worker
_orig_report = S.build_report


def letters_text(amp_as=''):
    """Elgar's own private writing before the cipher, '&' replaced by amp_as."""
    if not os.path.exists(LETTERS_FILE):
        return None
    with open(LETTERS_FILE, encoding='utf-8') as fh:
        raw = fh.read()
    out = []
    for part in raw.split('===== ')[1:]:
        head, _, body = part.partition('\n')
        date = head.split('|')[0].strip()
        if date < CIPHER_DATE and '| article |' not in head:
            out.append(body.replace('&c', ' etc ').replace('&', f' {amp_as} '))
    return '\n'.join(out)


def build_quad(corpus, style):
    if style != 'E3':
        return _orig_build_quad(corpus, style)
    import math
    from collections import Counter
    NL = 24
    text = S.styled_text(S.corpus_words(corpus), 'E2')
    lt = letters_text('')
    if lt is None:
        raise SystemExit(f"E3 needs {LETTERS_FILE} next to this script (or run with --only E0 E1 E2)")
    mine = S.fold(''.join(ABBR.get(w, w) for w in S.words_of_text(lt)))
    li = [S.LIDX[c] for c in text]
    cnt = Counter()
    for i in range(len(li) - 3):
        cnt[((li[i] * NL + li[i + 1]) * NL + li[i + 2]) * NL + li[i + 3]] += 1
    lj = [S.LIDX[c] for c in mine]
    for i in range(len(lj) - 3):
        cnt[((lj[i] * NL + lj[i + 1]) * NL + lj[i + 2]) * NL + lj[i + 3]] += 4
    tot = sum(cnt.values())
    floor = math.log(0.01 / tot)
    q = [floor] * (NL ** 4)
    for k, v in cnt.items():
        q[k] = math.log(v / tot)
    return q


def init_worker(corpus):
    _orig_init_worker(corpus)


def build_report(state, dicts):
    t = _orig_report(state, dicts)
    t = t.replace("DORABELLA 'PRIVATE ENGLISH' STYLES REPORT", "DORABELLA 'ELGAR'S PRIVATE ENGLISH' REPORT")
    t = t.replace('one of 6 styles', 'one of the 3 new styles (E1-E3)')
    t = t.replace('expected by chance among 6 styles', 'expected by chance among 3 styles')
    return t


S.build_quad = build_quad
S.init_worker = init_worker
S.build_report = build_report


def make_phrases(folder='Elgar_letters_quotes'):
    """Elgar's own words as a quote-test source: one copy with '&' left out, one with '&' as AND."""
    os.makedirs(folder, exist_ok=True)
    for name, amp in (('elgar_letters_amp_unwritten', ''), ('elgar_letters_amp_as_and', 'and')):
        t = letters_text(amp)
        if t is None:
            raise SystemExit(f"{LETTERS_FILE} not found next to this script")
        with open(os.path.join(folder, name + '.txt'), 'w', encoding='utf-8') as fh:
            fh.write(t)
    print(f"wrote {folder}\\ (2 files). Then:\n"
          f"  python .\\Dorabella_quotes.py --source {folder} --out elgarquotes_out --selftest\n"
          f"  python .\\Dorabella_quotes.py --source {folder} --out elgarquotes_out --stop-when-done")


if __name__ == '__main__':
    S.mp.freeze_support()
    if '--make-phrases' in sys.argv:
        make_phrases()
        sys.exit()
    if not any(a.startswith('--out') for a in sys.argv):
        sys.argv += ['--out', 'elgarstyle_out']
    if not any(a.startswith('--copies') for a in sys.argv):
        sys.argv += ['--copies', '1000']
    S.main()
