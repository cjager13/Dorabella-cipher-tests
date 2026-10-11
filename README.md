# Dorabella-cipher-tests

Fair statistical tests of Edward Elgar's 1897 Dorabella cipher.

**Result: no key tested turns the Dorabella into readable text.** About 30 tests in round 1 and 55 pre-registered test families in round 2 all land at chance. Known text put through the same tests is recovered almost every time.

The cipher is not random: it has strong habits in which symbol follows which. But scrambled copies that keep those habits do as well as the real cipher in every test, in every language and writing style tried.

The full write-up is [Testing the Dorabella Cipher.pdf](Testing%20the%20Dorabella%20Cipher.pdf).

## Why these tests are different

Most proposed Dorabella solutions show a decrypt that looks a bit like English and stop there. With 87 symbols and millions of possible keys, some key will always produce a few words by chance. Every test here therefore asks one question: **does the real cipher do better than scrambled copies of itself, put through exactly the same search?**

Two kinds of scrambled copy are used:

- **Plain shuffles.** The 87 symbols in random order.
- **Pair-preserving shuffles.** A random order that keeps every adjacent symbol pair exactly as in the Dorabella. These keep the writer's habits but carry no language, so beating them is the stronger test.

Each script also runs **planted text**: known text, enciphered the way that test assumes and put through the identical pipeline. That shows what a real message looks like and proves the method can find one. **Power** is the share of planted messages that stand out from the copies.

Round 2 added four rules:

- **Pre-registration.** Each script's ideas, keys and texts were written into the script before any real run.
- **Bonferroni.** With several families in a script, a result counts only below 0.05 divided by the number of families.
- **Median over seeds.** The real cipher's score is the median over 5 to 9 search seeds, never the luckiest.
- **Random-key comparisons and confirmation runs.** Meaningful numbers (a date, birthdays, scales) must beat random numbers of the same size, and anything near the line gets a fresh, pre-registered confirmation run.

p-values are (1 + copies at least as good) / (copies + 1), so they can never go below 1 / (copies + 1).

## Results, round 2 (1-11 October 2026)

| Idea tested | Script | p vs pair-preserving copies | Power | Verdict |
| --- | --- | --- | --- | --- |
| Transpositions built on 3, then any key (T3, TB) | `Dorabella_three.py` | 0.35, 0.84 | 98% | Chance |
| Elgar's key through transpositions or period-3 shifts and turns (TE, PA, PB) | `Dorabella_three.py` | 0.99, 0.46, 0.47 | 100% | Chance |
| Any key plus a period-3 shift or symbol turn (PC, PD) | `Dorabella_three.py` | 0.41, 0.64 | 85% | Chance |
| Three alphabets used in turn (PE) | `Dorabella_three.py` | 0.87 | 10% | Inconclusive: too much freedom |
| Loop counts alone, read in threes: repeats, counting keys, any key (L1, L2, L3) | `Dorabella_looptrits.py` | 0.048, 0.26, 0.31 | 83%, 100%, 48% | Chance; any-key test low power |
| Is it language? Nine visual habits vs English through Elgar's key (no search) | `Dorabella_language.py` | 0.18 vs English; 0.0002 vs random symbols | - | Looks like enciphered English, not random |
| 63 cribs (greetings, names, places), any key | `Dorabella_cribs.py` | 0.37 | 98% | Chance |
| Cribs with sets-of-three or loop-group keys | `Dorabella_cribs.py` | 1.00, 1.00 | 97-100% | Chance |
| Quotation from words Elgar set (18 songs and cantatas) | `Dorabella_quotes.py` | 0.75 | 98% | Chance |
| Quotation from Elgar's own letters to 1897 | `Dorabella_quotes.py` | 0.48 | 97% | Chance |
| Quotation from the Psalms: Coverdale, King James, Douay | `Dorabella_quotes.py` | 0.50, 0.47, 0.51 | 95-100% | Chance |
| Quotation from the 1662 Prayer Book, full KJV, full Douay | `Dorabella_quotes.py` | 0.19, 0.66, 0.58 | 95-100% | Chance |
| Quotation from the texts appointed for 14 July 1897 | `Dorabella_quotes.py` | 0.73 | 95% | Chance |
| Playful spellings: no vowels, backwards, phonetic, initials, pig Latin | `Dorabella_styles.py` | 0.26-0.88 | 92-100% | Chance |
| Elgar's private English: no AND, his abbreviations, his own letters | `Dorabella_elgarstyle.py` | 0.26-0.36 | 99-100% | Chance |
| The date, July 14 97, as the key (Elgar's key; any key; transposition) | `Dorabella_date.py` | 0.15, 0.38, 0.77 | 88-100% | Chance |
| French, German, Italian, Latin | `Dorabella_languages.py` | 0.14, 0.53, 0.075, 0.049 | 98-100% | Chance |
| Latin confirmation (Church Latin, pre-registered) | `Dorabella_languages.py --confirm-latin` | 0.049 | 100% | Not confirmed: decrypt unreadable |
| Musical scales as numbers | `Dorabella_scales.py --set scales` | 0.41, 0.86, 0.72 | 90-100% | Chance |
| Birthdays, ages, line lengths 29-31-27 | `Dorabella_scales.py --set numbers` | 0.63, 0.048, 0.12 | 90-100% | Chance: random numbers do as well |
| Note names CDEFGAB as the key | `Dorabella_scales.py --set notes` | 0.66, 0.91, 0.41 | 95-100% | Chance |
| Every circle-and-cross key (589,824) | `Dorabella_circlekeys.py` | 0.35 | 100% | Chance |
| Homophones: several symbols for one letter | `Dorabella_homophones.py` | 0.37, 0.63, 0.30 | 62-78% | No evidence; low power |

Where a row gives several p-values, they are the row's families in the order named. The "is it language?" row is not a key search: it asks whether the cipher's visual habits (loop changes, turns, mirror pairs, doubled symbols, repeats) look like real English put through Elgar's key. They do (p = 0.18 against 200,000 English passages), and they clearly differ from random symbols. So the cipher is consistent with an enciphered message; what no test found is the key or text that makes it one.

**The near-misses, and why none counts.**

- **Latin, p = 0.049 twice.** The confirmation rule was written into the script first: confirmed only if p < 0.05 *and* the decrypt reads as Latin. The decrypt ("AORAUIC FACIS QUETICHQCTI QUI TIMOLLORSI QUEM ... ENIM ... DAT") does not read, while planted Latin came back 95-99% right.
- **Birthdays with any key, p = 0.048.** It misses the corrected line (0.0167), and random numbers of the same size did as well (p = 0.35).
- **Loop-count repeats, p = 0.048.** The 29 loop-count triples repeat a little more than copies' do (232 appears five times), but this misses the line (0.0167), and no counting key turns them into text (p = 0.26, power 100%).
- **MY DEAR DORA, p = 0.0047.** DEAR DORA fits the start of line 3 read backwards, but it was noticed by looking, so it was reported as exploratory from the start. The text completed around it does not read.
- **Homophones and three alphabets** are the two tests too weak to rule their ideas out: at 87 symbols, that much freedom lets many keys produce equally English-looking text.

## Results, round 1 (27 September - 1 October 2026)

| Idea tested | Script | p vs pair-preserving copies | Verdict |
| --- | --- | --- | --- |
| Any one-symbol-one-letter key, 160 reading configurations | `Dorabella_marathon.py` | 0.085 (1,016 copies) | Pair habits explain it |
| Loop count picks the letter group (e.g. 1 loop = A-H), direction the letter: all 24,576 keys x 160 configs | `Dorabella_loopgroups.py` | 0.30 (200 copies); words 0.995 | Chance |
| Known English through the same pipeline | `Dorabella_marathon.py` | 102 of 102 controls score above the Dorabella | Not English-like |
| All 10,321,920 "sets of three" keys from Elgar's 1920 notebook, x 160 configs | `Dorabella_allkeys.py` | 0.84 and 0.14 (vs plain shuffles, 50 copies) | Chance |
| Keypad-style keys: 287 triplet and affine layouts x 160 configs | `Dorabella_tests.py` (Part C) | 0.72 (vs plain shuffles, 300 copies) | Chance |
| Keyword or phrase alphabets (DORABELLA, MISS PENNY, TOMORROW...) | `Dorabella_phrase.py`, `Dorabella_keyword.py` | chance-level hits only | Chance |
| Rare symbols are digits (the 1896 Courage cards) | `Dorabella_numbers.py` | 0.16 (vs plain shuffles, 100 copies) | Chance |
| Rechecked transcription, 10 uncertain symbols solved in the search, English | `Dorabella_v2.py` | 0.33 (200 copies) | Chance |
| Same, German | `Dorabella_v2.py --lang de` | 0.30 (200 copies) | Chance |
| Key learned on 2 lines, scored on the 3rd, one pre-registered config, typical of 24 seeds | `Dorabella_heldout_check.py` | 0.20 (240 copies) | Chance |

The Dorabella *does* beat plain shuffles (p = 0.0020 in 1,016), so its symbol order is not random. But pair-preserving shuffles match it, and it never finds more real words than chance.

## Files

### Shared libraries and data

| File | Purpose |
| --- | --- |
| `Dorabella_transcription.csv` | The transcription: 87 symbols as loop count + direction, with alternate readings for the uncertain symbols |
| `Dorabella_hillclimb.py` | Shared library: cipher, quadgram scorer, hill-climb solver. Also runs on its own. |
| `Dorabella_tests.py` | Shared library: reading orders and uncertain readings |
| `Dorabella_marathon.py` | Round 1's main battery; also a shared library for round 2 (control texts, pair-preserving shuffle) |
| `Dorabella_three.py` | Round 2's shared engine (symbols, Elgar's key, symbol turns, annealer, resumable runs); also its own test |

### Round 2 scripts

| Script | Purpose |
| --- | --- |
| `Dorabella_three.py` | Transpositions built on 3 and keys that repeat every 3 symbols: 8 families |
| `Dorabella_cribs.py` | 63 cribs fixed in place, with any key, sets-of-three keys or loop-group keys; MY DEAR DORA reported separately |
| `Dorabella_quotes.py` | Quotation test: phrases from a folder of texts placed where they fit, the rest completed |
| `Dorabella_styles.py` | Five playful spellings and standard English, each with its own scoring model |
| `Dorabella_elgarstyle.py` | Elgar's private English: no AND, his abbreviations, his own letters (built on `Dorabella_styles.py`) |
| `Dorabella_date.py` | The date July 14 97 as the key: shifts, symbol turns, loop steps, transpositions; random-number comparison |
| `Dorabella_scales.py` | Built on `Dorabella_date.py`: `--set scales`, `--set numbers` (birthdays, ages, line lengths) or `--set notes` (note names) |
| `Dorabella_languages.py` | French, German, Italian and Latin; `--confirm-latin` runs the pre-registered Church Latin confirmation |
| `Dorabella_circlekeys.py` | All 589,824 circle-and-cross keys (no search, needs numpy) |
| `Dorabella_homophones.py` | Keys where several symbols may share one letter |
| `Dorabella_fixtext.py`, `Dorabella_pdf2txt.py` | Helpers: extract and repair text from PDFs (ligatures, long s, split words) |
| `Dorabella_daytexts.py` | Helper: builds the folder of texts appointed for 14 July 1897 from a King James Bible and the Coverdale Psalter |

### Round 1 scripts

| File | Purpose |
| --- | --- |
| `Dorabella.py` | The first test: all 24,576 loop-group keys ranked by a simple bigram score (no scrambled-copy comparison) |
| `Dorabella_loopgroups.py` | The same 24,576 loop-group keys, with fair comparisons against both kinds of shuffle |
| `Dorabella_verify.py` | 160-configuration search vs 20 identical searches on shuffles |
| `Dorabella_allkeys.py` | Every "sets of three" key (8! orders x 2^8 up/down) |
| `Dorabella_phrase.py` | Keys built from words and phrases (Victorian keyword method) |
| `Dorabella_keyword.py` | 8-letter keywords using one letter from each set, plus near misses |
| `Dorabella_numbers.py` | Rarest symbols treated as digits |
| `Dorabella_v2.py` | Rechecked transcription with uncertain symbols solved inside the search; English or German |
| `Dorabella_heldout_check.py` | Reruns one held-out configuration over many seeds and copies |
| `Dorabella_keyword2.py` | Keyword test with direction order and up/down set by the keyword's first or third letters (output `keyword_out.txt`) |
| `Dorabella_language.py` | "Is it language?": were the symbols chosen by a message or by the hand? Compares the cipher's statistics with English passages and hand-made fakes (output `language_out/`) |
| `Dorabella_looptrits.py` | Do the loop counts alone, read three at a time, spell a message? Run in round 2 (output `looptrits_out/`) |

### Outputs

| Folder | Contents |
| --- | --- |
| `*_out/` | Each run's `report.txt` and `log.txt`, with final reports also saved as `*_final_report.txt`. `state.json` (every copy's score, for auditing or resuming) is included where it fits GitHub's size limit. |
| `*_final_report_old.txt` | Reports as first written, kept where a report was later regenerated (see Notes on the logs) |
| `early_runs/` | Outputs of the first day's scripts (27 September 2026) |

## Symbol codes

Each symbol is a two-digit code: **loop count (1-3), then direction (1-8)**. Direction 1 is the opening facing up, and directions step round in 45-degree turns. The common convention's direction 1 (opening to the right) is direction 3 here.

Elgar's standard key, applied to these codes, reproduces the published decrypt (BLTACEIARWUNISNF...) at 78 to 80 of 87 positions, depending on the transcription; the mismatches fall on uncertain symbols.

## Running it

**Requirements:** Python 3.9+ and numpy (`pip install numpy`).

Put all scripts, `Dorabella_transcription.csv` and the corpora in one folder. **Run each script's `--selftest` first:** it checks the setup, runs one planted message through the full search, and estimates the run time. Every round-2 script can be stopped with Ctrl+C and resumed with the same command; `--stop-when-done` makes it end when its pre-set number of copies is reached, so runs can be chained; `--report-only` rewrites the report from saved state.

### Corpora and text folders (included)

| File or folder | Used by | Source |
| --- | --- | --- |
| `english.txt` | all scripts | Bram Stoker, *Dracula* (1897), [Project Gutenberg #345](https://www.gutenberg.org/ebooks/345), the current release, plain text UTF-8 |
| `german.txt` | `Dorabella_v2.py --lang de`, `Dorabella_languages.py` | Fontane, *Effi Briest*, and Goethe, *Die Leiden des jungen Werthers* (both volumes), from Project Gutenberg |
| `french.txt` | `Dorabella_languages.py` | Hugo, *Les Misérables*, tome I, in French, from Project Gutenberg |
| `italian.txt` | `Dorabella_languages.py` | Manzoni, *I promessi sposi*, from Project Gutenberg |
| `latin.txt` | `Dorabella_languages.py` | Vulgate Daniel, Caesar, Cicero and Horace |
| `latin_church.txt` | `Dorabella_languages.py --confirm-latin` | Psalms and Gospels of the Clementine Vulgate ([eBible.org latVUC](https://ebible.org/find/details.php?id=latVUC)) |
| `quotes/` | `Dorabella_quotes.py` | Words of 18 songs and cantatas Elgar set |
| `Bible_KJV/`, `Bible_KJV_Psalms/` | `Dorabella_quotes.py`, `Dorabella_daytexts.py` | King James Bible, [Project Gutenberg #10](https://www.gutenberg.org/ebooks/10) |
| `Bible_Douay/`, `Bible_Douay_Psalms/` | `Dorabella_quotes.py` | Douay-Rheims Bible |
| `Bible_Psalter/`, `Bible_Psalter_Psalms/` | `Dorabella_quotes.py` | 1662 Book of Common Prayer and its Coverdale Psalter, from a scanned PDF repaired with `Dorabella_fixtext.py` |
| `Day_14July1897/` | `Dorabella_quotes.py` | Built by `Dorabella_daytexts.py` |

These are public-domain texts, included so every run can be reproduced exactly. `english.txt` is the exact file all round-2 scores came from. `Dorabella_languages.py` keeps every text between Project Gutenberg's START and END markers and drops the licence text. The English scorer reads `english.txt` as it is, so use the same file throughout.

### Commands and run times (8-core laptop)

| Command | Time |
| --- | --- |
| `python Dorabella_three.py --stop-when-done` | 78 h |
| `python Dorabella_cribs.py --stop-when-done` | 18.5 h |
| `python Dorabella_quotes.py --source <folder> --out <name>_out --stop-when-done` (add `--min-len 14` for the BCP and the full Bibles) | 2.7-11.9 h per source |
| `python Dorabella_styles.py --stop-when-done` | 18.5 h |
| `python Dorabella_elgarstyle.py --stop-when-done` | 8.3 h |
| `python Dorabella_date.py --stop-when-done` | 8.3 h |
| `python Dorabella_languages.py --stop-when-done`, then `--confirm-latin --stop-when-done` | 12.5 h, then 1.5 h |
| `python Dorabella_scales.py --set scales` (or `numbers`, `notes`) `--stop-when-done` | 6-7.5 h each |
| `python Dorabella_circlekeys.py --stop-when-done` | 38 min |
| `python Dorabella_homophones.py --stop-when-done` | 3.1 h |
| `python Dorabella_looptrits.py --workers 8` (not resumable; runs a fixed number of copies) | 6 min |
| `python Dorabella_hillclimb.py --selftest --corpus english.txt` | 10 s |
| `python Dorabella_verify.py` | 20 min |
| `python Dorabella_numbers.py` | 1.6 h |
| `python Dorabella_v2.py` (and `--lang de`) | 45 min each |
| `python Dorabella_loopgroups.py` | about 10 min |
| `python Dorabella_phrase.py` / `python Dorabella_keyword.py` | 30 s / seconds |
| `python Dorabella_allkeys.py --all-configs --nulls 50` | 16.8 h |
| `python Dorabella_heldout_check.py` | about 25 min |
| `python Dorabella_marathon.py` (Ctrl+C to stop; the same command resumes) | 18.6 h first pass; stopped at 1,016 copies after 46.6 h |

## Private material

Two round-2 tests used Elgar's own letters before 14 July 1897, transcribed from Jerrold Northrop Moore, *Edward Elgar: Letters of a Lifetime* (Oxford, 1990), pages 1-59. **That transcription is not in this repository, for copyright reasons.** The write-up shares only derived counts (for example, "&" 429 times against "and" 50 times) and the few short phrases that appear in the reports. Without the file, `Dorabella_elgarstyle.py` runs its first three styles (`--only E0 E1 E2`), and the letters quotation test cannot be rerun.

## Notes on the logs

- **Two transcriptions were used in round 1.** The marathon, allkeys, loopgroups, numbers, phrase and keyword runs used the first transcription. The `Dorabella_v2.py` runs and all of round 2 used the rechecked one in `Dorabella_transcription.csv`, whose `first_transcription_code` column shows the differences.
- **Position 32 was settled as 18 on 1 October 2026**, after a look at the highest-resolution scan; round 2 leaves 9 symbols uncertain and chooses them inside the search where the script allows it.
- **One `english.txt` for all of round 2.** Every round-2 report was scored with the same file, the current Gutenberg #345 release of *Dracula*. A different edition shifts raw scores by a few hundredths but not the p-values, since each run scores the cipher and its copies with the same file.
- **Reports regenerated for display only.** On 10 October 2026 the reports were rewritten with `--report-only` so that every "words:" line uses the large display dictionary (about 11,700 words) instead of a 309-word list. The display never affects scores. Earlier versions are kept as `*_final_report_old.txt`, and every changed line was checked by eye: only the timestamp and word lines differ.
- **Some early entries in `numbers_out/log.txt` should be ignored.** The runs before 14:15 on 2026-09-29 used an `english.txt` that had extra notes appended.
- **`numbers_out/` holds two different tests.** Round 1's `Dorabella_numbers.py` (rare symbols as digits) and round 2's `Dorabella_scales.py --set numbers` (birthdays, ages, line lengths) both wrote to this folder; the log shows which entries are which, and the round-2 report is `numbers_out_final_report.txt`.

## Open threads

- **Homophones, three alternating alphabets and any-key loop triples** cannot be ruled out by statistical search at this length (power 62-78%, 10% and 48%).
- **A blind re-reading of the uncertain symbols**, by someone who has not seen these results, would settle positions 32 and 67.
- **Not yet tested:** names of God and Jesus as cribs, street names and addresses as cribs, and the cipher's single dot as "&".

## Images

`Dorabella_cipher.jpg` is the 1897 cipher. The redrawings of Elgar's 1920 notebook page in the write-up are the investigator's own; the notebook itself is not reproduced here.

## Sources

- [Cipher Mysteries: Dorabella Cipher](https://ciphermysteries.com/other-ciphers/the-dorabella-cipher)
- [Cipher Mysteries: The Dorabella cipher and Elgar's other little ciphertexts](https://ciphermysteries.com/2013/10/09/elgars-little-cipher)
- [Wikipedia: Dorabella Cipher](https://en.wikipedia.org/wiki/Dorabella_Cipher)
- Jerrold Northrop Moore, *Edward Elgar: Letters of a Lifetime* (Oxford, 1990)
- [Project Gutenberg #345: Dracula](https://www.gutenberg.org/ebooks/345) and [#10: The King James Bible](https://www.gutenberg.org/ebooks/10)
- [eBible.org: Latin Clementine Vulgate](https://ebible.org/find/details.php?id=latVUC)
