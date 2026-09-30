# Dorabella-cipher-tests

Fair statistical tests of Edward Elgar's 1897 Dorabella cipher.

**Result: no simple substitution key turns the Dorabella into English or German.** Across about 30 tests, the cipher scores at chance level, while known English put through the same tests is recovered every time.

The full write-up is [Testing the Dorabella Cipher.pdf](Testing%20the%20Dorabella%20Cipher.pdf).

## Why these tests are different

Most proposed Dorabella solutions show a decrypt that looks a bit like English and stop there. With 87 symbols and millions of possible keys, some key will always produce a few words by chance. Every test here therefore asks one question: **does the real cipher do better than scrambled copies of itself, put through exactly the same search?**

Two kinds of scrambled copy are used:

- **Plain shuffles.** The 87 symbols in random order.
- **Pair-preserving shuffles.** A random order that keeps every adjacent symbol pair exactly as in the Dorabella. These keep the writer's habits but carry no language, so beating them is the stronger test.

Each script also runs **planted text**: known English or German, enciphered with a random key and put through the identical pipeline. That shows what a real message looks like and proves the method can find one.

## Results

| Idea tested | Script | p vs pair-preserving copies | Verdict |
| --- | --- | --- | --- |
| Any one-symbol-one-letter key, 160 reading configurations | `Dorabella_marathon.py` | 0.085 (400 copies) | Pair habits explain it |
| Loop count picks the letter group (e.g. 1 loop = A-H), direction the letter: all 24,576 keys x 160 configs | `Dorabella_loopgroups.py` | 0.30 (200 copies); words 0.995 | Chance |
| Known English through the same pipeline | `Dorabella_marathon.py` | 41 of 41 controls score above the Dorabella | Not English-like |
| All 10,321,920 "sets of three" keys from Elgar's 1920 notebook, x 160 configs | `Dorabella_allkeys.py` | 0.84 and 0.14 (vs plain shuffles, 50 copies) | Chance |
| Keypad-style keys: 287 triplet and affine layouts x 160 configs | `Dorabella_tests.py` (Part C) | 0.72 (vs plain shuffles, 300 copies) | Chance |
| Keyword or phrase alphabets (DORABELLA, MISS PENNY, TOMORROW...) | `Dorabella_phrase.py`, `Dorabella_keyword.py` | chance-level hits only | Chance |
| Rare symbols are digits (the 1896 Courage cards) | `Dorabella_numbers.py` | 0.16 (vs plain shuffles, 100 copies) | Chance |
| Rechecked transcription, 10 uncertain symbols solved in the search, English | `Dorabella_v2.py` | 0.33 (200 copies) | Chance |
| Same, German | `Dorabella_v2.py --lang de` | 0.30 (200 copies) | Chance |

The Dorabella *does* beat plain shuffles (p = 0.0025 in 400), so its symbol order is not random. But pair-preserving shuffles match it, and it never finds more real words than chance.

## Files

| File | Purpose |
| --- | --- |
| `Dorabella_transcription.csv` | The transcription: 87 symbols as loop count + direction, with alternate readings for the 10 uncertain symbols |
| `Dorabella_hillclimb.py` | Shared library: cipher, quadgram scorer, hill-climb solver. Also runs on its own. |
| `Dorabella_tests.py` | Shared library: reading orders and uncertain readings (first transcription) |
| `Dorabella.py` | The first test: all 24,576 loop-group keys ranked by a simple bigram score (no scrambled-copy comparison) |
| `Dorabella_loopgroups.py` | The same 24,576 loop-group keys, with fair comparisons against both kinds of shuffle |
| `Dorabella_verify.py` | 160-configuration search vs 20 identical searches on shuffles |
| `Dorabella_marathon.py` | The main battery: both shuffle types, known-English controls, held-out test, simulated-annealing search. Resumable; runs until stopped. |
| `Dorabella_allkeys.py` | Every "sets of three" key (8! orders x 2^8 up/down) |
| `Dorabella_phrase.py` | Keys built from words and phrases (Victorian keyword method) |
| `Dorabella_keyword.py` | 8-letter keywords using one letter from each set, plus near misses |
| `Dorabella_numbers.py` | Rarest symbols treated as digits |
| `Dorabella_v2.py` | Rechecked transcription with uncertain symbols solved inside the search; English or German |
| `*_out/` | Logs and reports from the runs quoted in the write-up |
| `early_runs/` | Outputs of the first day's scripts (27 September 2026), summarized below |

## Symbol codes

Each symbol is a two-digit code: **loop count (1-3), then direction (1-8)**. Direction 1 is the opening facing up. The common convention's direction 1 (opening to the right) is direction 3 here.

Elgar's standard key, applied to these codes, reproduces the published decrypt (BLTACEIARWUNISNF...) at 80 of 87 positions. The mismatches fall on uncertain symbols.

## Running it

**Requirements:** Python 3.9+. `Dorabella_allkeys.py` also needs numpy (`pip install numpy`).

**Corpora** (not included; public domain from Project Gutenberg):

- `english.txt`: Bram Stoker, *Dracula* (1897), [eBook #345](https://www.gutenberg.org/ebooks/345), plain text UTF-8.
- `german.txt` (for `--lang de` only): Theodor Fontane, *Effi Briest* (1895), [eBook #5323](https://www.gutenberg.org/ebooks/5323), optionally joined with other 1890s German novels. Umlauts and Gutenberg's licence text are handled by the script.

Put all scripts and the corpora in one folder. Run each script's self-test first. Times below are on an 8-core laptop.

| Command | Time |
| --- | --- |
| `python Dorabella_hillclimb.py --selftest --corpus english.txt` | 10 s |
| `python Dorabella_verify.py` | 20 min |
| `python Dorabella_numbers.py --selftest --nulls 5`, then `python Dorabella_numbers.py` | 4 min, then 1.6 h |
| `python Dorabella_v2.py --selftest --nulls 8`, then `python Dorabella_v2.py` | 2 min, then 45 min |
| `python Dorabella_v2.py --lang de` | 45 min |
| `python Dorabella_loopgroups.py` (needs numpy) | about 10 min |
| `python Dorabella_phrase.py` / `python Dorabella_keyword.py` | 30 s / seconds |
| `python Dorabella_allkeys.py --all-configs --nulls 50` | 16.8 h |
| `python Dorabella_marathon.py` (Ctrl+C to stop; the same command resumes) | 18.6 h per pass |

## Earlier tests (27 September 2026)

The first day's scripts led to the tests above; their outputs are in `early_runs/`.

| Script and output | What it did | Result |
| --- | --- | --- |
| `Dorabella.py` -> `Dorabella_out.txt` | Ranked all 24,576 loop-group keys with a simple bigram score | Best decrypts gibberish (MANETOWEZLDIJXSIJ...); no fair comparison, now done by `Dorabella_loopgroups.py` |
| `Dorabella_hillclimb.py --selftest` -> `Dorabella_hillclimb_selftest.txt` | Planted English, random key | 95% of letters recovered |
| `Dorabella_hillclimb.py` -> `Dorabella_hillclimb_real.txt` | Solver on the real cipher, 200 restarts, 100 shuffles | p = 0.059 (score), 0.41 (words) |
| `Dorabella_tests.py` -> `Dorabella_tests_out.txt`, Part C | 287 keypad-style keys | p = 0.85 (as transcribed), 0.72 (all 160 configs) |
| same, Part A | Planted English with 0, 2, 4, 6, 8 misread symbols | 68%, 87%, 48%, 51%, 30% of letters recovered: a few misreads hurt a lot, which motivated the transcription recheck |
| same, Part B | 160 configs vs a single-search null | z = 4.98, but that null was mismatched; corrected by `Dorabella_verify.py` |
| `Dorabella_verify.py` -> `Dorabella_verify_out.txt` | 160 configs vs 20 identical searches on shuffles | 0 of 20 matched (p = 0.048, the floor); led to the 400-copy marathon |

`Dorabella_scorer.py` (loop-group keys with a word-based score) was also run that day, but its output was not saved; `Dorabella_loopgroups.py` replaces it.

## Notes on the logs

- **Two transcriptions were used.** The marathon, allkeys, loopgroups, numbers, phrase and keyword runs used the first transcription: five uncertain symbols, and position 85 read as 14. The `Dorabella_v2.py` runs used the rechecked one in `Dorabella_transcription.csv`. That file's `first_transcription_code` column shows the only difference.
- **Some early entries in `numbers_out/log.txt` should be ignored.** The runs before 14:15 on 2026-09-29 used an `english.txt` that had extra notes appended. They are superseded by the later entries in the same log.
- **The p-values have a floor.** A p-value can never go below 1 / (copies + 1). "At the floor" means no copy matched the real cipher yet.

## Open threads

- **Position 32 read as 18.** This reading helps more than chance against pair-preserving shuffles (p = 0.0075), and both the English and German searches chose it independently. It is worth rechecking on the original first.
- **Reversed reading orders** beat pair-preserving shuffles modestly (p about 0.015). That is expected by chance among the roughly ten such comparisons made.
- **A scorer for Elgar's playful spelling** (phonetic forms, abbreviations) is the natural next test.

## Sources

- [Cipher Mysteries: Dorabella Cipher](https://ciphermysteries.com/other-ciphers/the-dorabella-cipher)
- [Cipher Mysteries: The Dorabella cipher and Elgar's other little ciphertexts](https://ciphermysteries.com/2013/10/09/elgars-little-cipher)
- [Wikipedia: Dorabella Cipher](https://en.wikipedia.org/wiki/Dorabella_Cipher)
