#!/usr/bin/env python3
"""
Turn a typeset (not scanned) PDF into a plain-text file for Dorabella_quotes.py.

  pip install pypdf
  python Dorabella_pdf2txt.py the-book-of-common-prayer-1662.pdf Bible_Psalter\\BCP_1662.txt

Joins words split by end-of-line hyphens and replaces typographic ligatures (fi, fl...).
Prints a sample so you can check the text is clean before using it.
"""
import re
import sys

from pypdf import PdfReader

LIG = {'ﬀ': 'ff', 'ﬁ': 'fi', 'ﬂ': 'fl', 'ﬃ': 'ffi', 'ﬄ': 'ffl',
       'ﬅ': 'st', 'ﬆ': 'st',                     # long-s + t, and s + t ligatures
       'ſ': 's', '’': "'", '‘': "'"}


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: python Dorabella_pdf2txt.py input.pdf output.txt")
    reader = PdfReader(sys.argv[1])
    parts = []
    for page in reader.pages:
        parts.append(page.extract_text() or '')
    text = '\n'.join(parts)
    text = re.sub('(\\w) ([\ufb05\ufb06])', r'\1\2', text)   # 'ha \ufb05' -> 'ha\ufb05': gap before the ligature
    for a, b in LIG.items():
        text = text.replace(a, b)
    text = re.sub(r'(\w)-\n(\w)', r'\1\2', text)        # re-join hyphenated words
    with open(sys.argv[2], 'w', encoding='utf-8') as fh:
        fh.write(text)
    words = re.findall(r'[A-Za-z]+', text)
    print(f"{len(reader.pages)} pages, {len(text):,} characters, {len(words):,} words -> {sys.argv[2]}")
    mid = len(text) // 2
    print("sample from the middle:\n" + text[mid:mid + 400])


if __name__ == '__main__':
    main()
