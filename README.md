# Dorabella-cipher-tests

Fair statistical tests of Edward Elgar's 1897 Dorabella cipher.

**Result: no key tested turns the Dorabella into readable text.** About 30 tests in round 1, 58 pre-registered test families in round 2 and 17 more in the November 2026 audit: all but three land at chance, and none of the three reads. Known text put through the same tests is recovered almost every time. The three exceptions trace to local quirks of the symbols (a cluster of repeats in line 1, and one symbol spaced like word gaps in line 3), not to a message. About 70 families at their own lines would be expected to give two or three such passes by chance.

The cipher is not random: it has strong habits in which symbol follows which. But scrambled copies that keep those habits do as well as the real cipher in almost every test, in every language, writing style and key tried. One habit is now firmly established: neighbouring symbols are far more often the same shape turned half round than in any enciphered English (p = 0.001), which replicates Hauer et al. (2021).

**Audit, November 2026.** A review of every file found that the tests built on "Elgar's own key" had used the wrong key (a modern reconstruction, not the key in Elgar's notebook), that symbol 32 had been read differently from every published transcription, and that two ideas (a word-gap symbol, and the cipher as a tune) had never been tested. All of these were rerun or tested; see [Audit, November 2026](#audit-november-2026). Nothing changed the conclusion.

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

## Audit, November 2026

### What was wrong

- **"Elgar's key" was not Elgar's key.** Rounds 1 and 2 took the key from the decrypt BLTACEIARWUNISNF... printed on Cipher Mysteries. That comes from a tidy pigpen-style reconstruction (ABC DEF GHI JKL on one cross, MNO PQR STU WXY on the other, some groups running backwards). Elgar's notebook, and the key table on `Dorabella_cipher.jpg`, write the alphabet straight round the eight directions: ABC opening right, DEF down-right, GHI down, KLM down-left, NOP left, QRS up-left, TUW up, XYZ up-right, 1 loop = first letter. That key gives the widely printed reading BPECAHTCKYFRQDRIRR... (reproduced here at 85 of 87 symbols). Every fixed-key test was affected: TE, PA, PB, the DE family of the date, scales, numbers and notes tests, tick-shift, and the main comparison of the "is it language?" test. Tests that search for ANY key were not affected, and the sets-of-three key spaces (`Dorabella_allkeys.py`, cribs C2) already contained the notebook layout.
- **Symbol 32.** It was fixed as code 18 on 1 October 2026; every published reading, including the five-transcription consensus of Hauer et al. (2025), has code 11 (one loop, opening up), and no search ever tried 11 there. The reruns use the consensus reading (symbols 23, 32, 34, 85, 86 read 14, 11, 12, 14, 12); the word-gap test treats 32 as uncertain.
- **The "is it language?" test missed the known oddity.** Its mirror features did not count neighbours that are the same shape turned half round, the excess Hauer et al. (2021) and Massey reported, and its main comparison used the wrong key.
- **Never tested:** one symbol as the gap between words, and the cipher read directly as a tune.
- **Checked and sound:** the pair-preserving shuffle, p-values, seeds, controls, and the privacy of the letters (no letters state file is in the repository). A calibration run showed planted English beats all of its own pair-preserving copies (10 of 10).

### Results

| Idea tested | Script | p vs pair-preserving copies | Power | Verdict |
| --- | --- | --- | --- | --- |
| Notebook key through every transposition, + period-3 letter shift, + period-3 symbol turns (TE, PA, PB) | `Dorabella_notebook.py three` | 0.82, 0.41, 0.29 | 100% | Chance |
| Notebook key + the date, scales, birthdays/ages/line lengths, note names (DE family of each) | `Dorabella_notebook.py date` / `scales` / `numbers` / `notes` | 0.94, 0.27, 0.50, 0.78 | 100% | Chance; random numbers do as well |
| Ticks as shifts on the notebook layout (512 keys) | `Dorabella_notebook.py tickshift` | 0.93 | 100% | Chance |
| Every walk of the alphabet round the eight directions (4,096 keys; the notebook key is one) | `Dorabella_notebook.py spokes` | 0.39 | 100% | Chance |
| Is it language? F1-F9 against English through the notebook key (P1) | `Dorabella_language2.py` | 0.033 (line 0.0167) | - | Not past the line |
| Half-turn pairs, same loops (P2) / any loops (P3), against English through the notebook key | `Dorabella_language2.py` | 0.0064 / 0.0014 (vs any key: 0.001 / 0.0002) | - | **Real habit, replicates Hauer et al. 2021** |
| The cipher as a tune: scale round the circle; any note assignment; any key (MS, MP, MK) | `Dorabella_melody.py` | 0.41, 0.19, 0.64 | 95-100% | Chance |
| One symbol is the gap between words (WS) | `Dorabella_spaces.py` | 0.005 (line 0.05) | 100% | Passes, but unreadable: one symbol in line 3 |
| Ordinary English with symbol 32 reopened (EN, reference) | `Dorabella_spaces.py` | 0.29 | 95% | Chance, as before |

Fixed-key reruns use the consensus reading; any-key searches choose the uncertain symbols inside the search.

**The word-gap pass.** The search chose symbol 38 (three loops, opening up-left) as the gap. It ends line 1 and line 2, and in line 3 it recurs every few symbols, cutting that line into "words" of 1, 2, 5, 5, 4 and 5 letters, a very English rhythm that pair-preserving copies rarely have. But lines 1 and 2 have no 38 inside them, so they would be single 28- and 30-letter words, and the decrypt does not read: "...SAUCUTELLNEEMETHESCRIPNFIDLY SURRUTH ONESEPHANDEEDIGRIPN... ALIT HONED ABOAD ON D". By the rule written into the script before the run, it is a statistical oddity, not a solution. The 38s were looked at after the run, so that explanation is exploratory.

**Half-turn pairs.** 13 neighbouring pairs have the same loop count and face opposite ways (27 face opposite ways whatever their loops), against 5.8 and 15.3 on average for English through the notebook key. These features were counted during the audit, before the script was written, so this is a replication of the published result rather than a new discovery. It is a habit of the hand or a rule of the cipher that no test here models; it is not a decryption.

## Results, round 2 (1-11 October 2026)

Rows marked * used the wrong base key ("Elgar's key" was the Cipher Mysteries reconstruction, not the notebook key); they were rerun in the audit above. The key-based verdicts stayed at chance; the "is it language?" verdict changed (see the audit).

| Idea tested | Script | p vs pair-preserving copies | Power | Verdict |
| --- | --- | --- | --- | --- |
| Transpositions built on 3, then any key (T3, TB) | `Dorabella_three.py` | 0.35, 0.84 | 98% | Chance |
| Elgar's key through transpositions or period-3 shifts and turns (TE, PA, PB) * | `Dorabella_three.py` | 0.99, 0.46, 0.47 | 100% | Chance |
| Any key plus a period-3 shift or symbol turn (PC, PD) | `Dorabella_three.py` | 0.41, 0.64 | 85% | Chance |
| Three alphabets used in turn (PE) | `Dorabella_three.py` | 0.87 | 10% | Inconclusive: too much freedom |
| Loop counts alone, read in threes: repeats, counting keys, any key (L1, L2, L3) | `Dorabella_looptrits.py` | 0.048, 0.26, 0.31 | 83%, 100%, 48% | Chance; any-key test low power |
| Is it language? Nine visual habits vs English through Elgar's key (no search) * | `Dorabella_language.py` | 0.18 vs English; 0.0002 vs random symbols | - | Looks like enciphered English, not random |
| 63 cribs (greetings, names, places), any key | `Dorabella_cribs.py` | 0.37 | 98% | Chance |
| Cribs with sets-of-three or loop-group keys | `Dorabella_cribs.py` | 1.00, 1.00 | 97-100% | Chance |
| 238 more cribs: God, Jesus, prayers, churches, homes and streets (C1, C2, C3) | `Dorabella_morecribs.py` | 0.009, 0.77, 1.00 | 95-100% | C1 passes the line; decrypt unreadable |
| Crib check: what the crib adds over no crib, Dracula English and Elgar's voice (D1, DE) | `Dorabella_cribcheck.py` | 0.034, 0.002 | 72%, 60% | DE passes; traced to one symbol cluster, decrypt unreadable |
| Quotation from words Elgar set (18 songs and cantatas) | `Dorabella_quotes.py` | 0.75 | 98% | Chance |
| Quotation from Elgar's own letters to 1897 | `Dorabella_quotes.py` | 0.48 | 97% | Chance |
| Quotation from the Psalms: Coverdale, King James, Douay | `Dorabella_quotes.py` | 0.50, 0.47, 0.51 | 95-100% | Chance |
| Quotation from the 1662 Prayer Book, full KJV, full Douay | `Dorabella_quotes.py` | 0.19, 0.66, 0.58 | 95-100% | Chance |
| Quotation from the texts appointed for 14 July 1897 | `Dorabella_quotes.py` | 0.73 | 95% | Chance |
| Playful spellings: no vowels, backwards, phonetic, initials, pig Latin | `Dorabella_styles.py` | 0.26-0.88 | 92-100% | Chance |
| Elgar's private English: no AND, his abbreviations, his own letters | `Dorabella_elgarstyle.py` | 0.26-0.36 | 99-100% | Chance |
| The date, July 14 97, as the key (Elgar's key; any key; transposition) * | `Dorabella_date.py` | 0.15, 0.38, 0.77 | 88-100% | Chance |
| French, German, Italian, Latin | `Dorabella_languages.py` | 0.14, 0.53, 0.075, 0.049 | 98-100% | Chance |
| Latin confirmation (Church Latin, pre-registered) | `Dorabella_languages.py --confirm-latin` | 0.049 | 100% | Not confirmed: decrypt unreadable |
| Musical scales as numbers * | `Dorabella_scales.py --set scales` | 0.41, 0.86, 0.72 | 90-100% | Chance |
| Birthdays, ages, line lengths 29-31-27 * | `Dorabella_scales.py --set numbers` | 0.63, 0.048, 0.12 | 90-100% | Chance: random numbers do as well |
| Note names CDEFGAB as the key * | `Dorabella_scales.py --set notes` | 0.66, 0.91, 0.41 | 95-100% | Chance |
| Every circle-and-cross key (589,824; does not include the notebook's walk-round layout, tested in the audit) | `Dorabella_circlekeys.py` | 0.35 | 100% | Chance |
| Ticks as shifts: Elgar's key with each group moved up or down by 1 or by the loop count (1,024 keys) * | `Dorabella_tickshift.py` | 0.10 | 100% | Chance |
| Homophones: several symbols for one letter | `Dorabella_homophones.py` | 0.37, 0.63, 0.30 | 62-78% | No evidence; low power |

Where a row gives several p-values, they are the row's families in the order named. The "is it language?" row is not a key search: it asks whether the cipher's visual habits (loop changes, turns, mirror pairs, doubled symbols, repeats) look like real English put through Elgar's key. They do (p = 0.18 against 200,000 English passages), and they clearly differ from random symbols. So the cipher is consistent with an enciphered message; what no test found is the key or text that makes it one. *(Audit: with Elgar's notebook key this test gives p = 0.033, and the half-turn pairs it never measured are far too many for enciphered English; see the audit.)*

**The near-misses and the one pass, and why none counts.**

- **Crib WELLS, p = 0.009, then 0.002 in the crib check.** The second crib list (238 eligible cribs) passed its line (0.0167) with WELLS: "...DORA F FAT HOUSE ... DEE DIG FIN SHE NORWAY A WELLS EE MET HER...". The follow-up `Dorabella_cribcheck.py` ran the same search with no crib on the same 1,000 copies. Without a crib the real cipher has no edge (p = 0.13); the crib simply costs it almost nothing (0.05 against 0.36 for copies). With Elgar's own voice as the scoring model that gap passes (p = 0.002, line 0.025), but the decrypt still does not read, and planted messages come back about 97% right. The cause, found by looking afterwards: line 1, symbols 15-23, has the pattern a b a a c d d a e (one symbol four times, two doubled pairs), which read backwards falls into English as WELLSEEME or EDDSEEME. Only 1 of the 1,000 pair-preserving copies has that pattern. Any search is drawn to "...SEE ME..." there, and WELLS fits it for free.
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
| `Dorabella_morecribs.py` | Built on `Dorabella_cribs.py`: 343 religious, place and address cribs plus single words from multi-word cribs (238 eligible) |
| `Dorabella_cribcheck.py` | What each crib is worth over the same search with no crib, in Dracula English and in Elgar's voice; reads `morecribs_out/state.json` and rebuilds its copies |
| `Dorabella_quotes.py` | Quotation test: phrases from a folder of texts placed where they fit, the rest completed |
| `Dorabella_styles.py` | Five playful spellings and standard English, each with its own scoring model |
| `Dorabella_elgarstyle.py` | Elgar's private English: no AND, his abbreviations, his own letters (built on `Dorabella_styles.py`) |
| `Dorabella_date.py` | The date July 14 97 as the key: shifts, symbol turns, loop steps, transpositions; random-number comparison |
| `Dorabella_scales.py` | Built on `Dorabella_date.py`: `--set scales`, `--set numbers` (birthdays, ages, line lengths) or `--set notes` (note names) |
| `Dorabella_languages.py` | French, German, Italian and Latin; `--confirm-latin` runs the pre-registered Church Latin confirmation |
| `Dorabella_circlekeys.py` | All 589,824 circle-and-cross keys (no search, needs numpy) |
| `Dorabella_tickshift.py` | Elgar's key with each group's letters moved forward (ticks out) or back (ticks in), by 1 or by the loop count: 1,024 keys, no search |
| `Dorabella_homophones.py` | Keys where several symbols may share one letter |
| `Dorabella_fixtext.py`, `Dorabella_pdf2txt.py` | Helpers: extract and repair text from PDFs (ligatures, long s, split words) |
| `Dorabella_daytexts.py` | Helper: builds the folder of texts appointed for 14 July 1897 from a King James Bible and the Coverdale Psalter |

### Audit scripts (November 2026)

| Script | Purpose |
| --- | --- |
| `Dorabella_notebook.py` | Reruns every fixed-key test with Elgar's notebook key and the consensus reading: `three`, `date`, `scales`, `numbers`, `notes`, `tickshift`, plus the new `spokes` family. `--reading ours` repeats a test with this project's own readings |
| `Dorabella_language2.py` | "Is it language?" with the notebook key, half-turn pairs added, symbol 32 treated as uncertain |
| `Dorabella_spaces.py` | One symbol as the gap between words (any key), and run-together English again with symbol 32 reopened |
| `Dorabella_melody.py` | The cipher read as a tune: direction = scale step or rest, loops = note length; three families against about 8,500 folk songs |

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

**Elgar's notebook key** in these codes: ABC on direction 3 (opening right), DEF 4, GHI 5, KLM 6, NOP 7, QRS 8, TUW 1, XYZ 2, with 1 loop = the first letter of the group (I/J and U/V share symbols). On the consensus reading it reproduces the widely printed decrypt BPECAHTCKYFRQDRIRR... at 85 of 87 positions. The key the round 1-2 scripts call `ELGAR_SETS` (giving BLTACEIARWUNISNF...) is a modern reconstruction, not Elgar's; see the audit.

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

**Folk tunes for `Dorabella_melody.py` (not included).** The Essen folksong collection in ABC notation, from the music21 project: download https://github.com/cuthbertLab/music21 and copy its `music21/corpus/essenFolksong` folder next to the scripts. Its licence allows non-commercial use, so neither that folder nor the `melody.txt` built from it is in this repository.

These are public-domain texts, included so every run can be reproduced exactly. `english.txt` is the exact file all round-2 scores came from. `Dorabella_languages.py` keeps every text between Project Gutenberg's START and END markers and drops the licence text. The English scorer reads `english.txt` as it is, so use the same file throughout.

### Commands and run times (8-core laptop)

| Command | Time |
| --- | --- |
| `python Dorabella_three.py --stop-when-done` | 78 h |
| `python Dorabella_cribs.py --stop-when-done` | 18.5 h |
| `python Dorabella_morecribs.py --stop-when-done` | 34 h |
| `python Dorabella_cribcheck.py --stop-when-done` (after morecribs; needs the private letters file) | 16.4 h |
| `python Dorabella_quotes.py --source <folder> --out <name>_out --stop-when-done` (add `--min-len 14` for the BCP and the full Bibles) | 2.7-11.9 h per source |
| `python Dorabella_styles.py --stop-when-done` | 18.5 h |
| `python Dorabella_elgarstyle.py --stop-when-done` | 8.3 h |
| `python Dorabella_date.py --stop-when-done` | 8.3 h |
| `python Dorabella_languages.py --stop-when-done`, then `--confirm-latin --stop-when-done` | 12.5 h, then 1.5 h |
| `python Dorabella_scales.py --set scales` (or `numbers`, `notes`) `--stop-when-done` | 6-7.5 h each |
| `python Dorabella_circlekeys.py --stop-when-done` | 38 min |
| `python Dorabella_tickshift.py --stop-when-done` | under 2 min |
| `python Dorabella_homophones.py --stop-when-done` | 3.1 h |
| `python Dorabella_looptrits.py --workers 8` (not resumable; runs a fixed number of copies) | 6 min |
| `python Dorabella_notebook.py three --stop-when-done` (or `date`, `scales`, `numbers`, `notes`, `tickshift`, `spokes`) | 20 min for `three`; about a minute each for the others |
| `python Dorabella_language2.py` | a few minutes |
| `python Dorabella_spaces.py --stop-when-done` | 3.5 h |
| `python Dorabella_melody.py --build essenFolksong`, then `python Dorabella_melody.py MS --stop-when-done` (or `MP`, `MK`) | seconds, 16 min, 1.6 h |
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

Three round-2 tests used Elgar's own letters before 14 July 1897, transcribed from Jerrold Northrop Moore, *Edward Elgar: Letters of a Lifetime* (Oxford, 1990), pages 1-59. **That transcription is not in this repository, for copyright reasons.** The write-up shares only derived counts (for example, "&" 429 times against "and" 50 times) and the few short phrases that appear in the reports. Without the file, `Dorabella_elgarstyle.py` runs its first three styles (`--only E0 E1 E2`), and the letters quotation test and `Dorabella_cribcheck.py` cannot be rerun.

## Notes on the logs

- **Two transcriptions were used in round 1.** The marathon, allkeys, loopgroups, numbers, phrase and keyword runs used the first transcription. The `Dorabella_v2.py` runs and all of round 2 used the rechecked one in `Dorabella_transcription.csv`, whose `first_transcription_code` column shows the differences.
- **Position 32 was settled as 18 on 1 October 2026**, after a look at the highest-resolution scan; round 2 leaves 9 symbols uncertain and chooses them inside the search where the script allows it. Every published transcription reads 11 there; the audit reruns use 11, and reopening it changed nothing (run-together English p = 0.29).
- **`nb_*_out/` reports show only the families that were rerun.** The other families in those reports say "waiting for copies"; they search for any key and did not need rerunning.
- **One `english.txt` for all of round 2.** Every round-2 report was scored with the same file, the current Gutenberg #345 release of *Dracula*. A different edition shifts raw scores by a few hundredths but not the p-values, since each run scores the cipher and its copies with the same file.
- **Reports regenerated for display only.** On 10 October 2026 the reports were rewritten with `--report-only` so that every "words:" line uses the large display dictionary (about 11,700 words) instead of a 309-word list. The display never affects scores. Earlier versions are kept as `*_final_report_old.txt`, and every changed line was checked by eye: only the timestamp and word lines differ.
- **Some early entries in `numbers_out/log.txt` should be ignored.** The runs before 14:15 on 2026-09-29 used an `english.txt` that had extra notes appended.
- **`numbers_out/` holds two different tests.** Round 1's `Dorabella_numbers.py` (rare symbols as digits) and round 2's `Dorabella_scales.py --set numbers` (birthdays, ages, line lengths) both wrote to this folder; the log shows which entries are which, and the round-2 report is `numbers_out_final_report.txt`.

## Open threads

- **Half-turn pairs.** The cipher's one firm oddity. Whether it is a habit of Elgar's hand or a rule of the cipher (Hauer et al. suggest it might mark doubled letters) is untested.
- **Symbol 38 in line 3**, spaced like word gaps, but absent inside lines 1 and 2.
- **Homophones, three alternating alphabets and any-key loop triples** cannot be ruled out by statistical search at this length (power 62-78%, 10% and 48%).
- **A blind re-reading of the uncertain symbols**, by someone who has not seen these results, would settle positions 32 and 67.
- **Not yet tested:** the cipher's single dot as "&". (Names of God and Jesus, prayers, churches, and street names and addresses were tested as cribs in `Dorabella_morecribs.py`.)

## Images

`Dorabella_cipher.jpg` is the 1897 cipher. The redrawings of Elgar's 1920 notebook page in the write-up are the investigator's own; the notebook itself is not reproduced here.

## Sources

- [Cipher Mysteries: Dorabella Cipher](https://ciphermysteries.com/other-ciphers/the-dorabella-cipher)
- [Cipher Mysteries: The Dorabella cipher and Elgar's other little ciphertexts](https://ciphermysteries.com/2013/10/09/elgars-little-cipher)
- [Wikipedia: Dorabella Cipher](https://en.wikipedia.org/wiki/Dorabella_Cipher)
- Hauer, Choi, Sundar, Hindle, Smallwood and Kondrak, [Experimental Analysis of the Dorabella Cipher with Statistical Language Models](https://ecp.ep.liu.se/index.php/histocrypt/article/download/159/115) (HistoCrypt 2021)
- Hauer, Choi, Hindle, Smallwood and Kondrak, [Dorabella Cipher as Musical Inspiration](https://arxiv.org/abs/2509.17950) (2025): the consensus transcription
- [Benzedrine: Clues to the Dorabella Cipher](https://www.benzedrine.ch/dorabella.html) and [dCode](https://www.dcode.fr/dorabella-cipher): the notebook-key reading
- [music21 corpus](https://github.com/cuthbertLab/music21): Essen folksong collection (ABC encodings by Seymour Shlien)
- Jerrold Northrop Moore, *Edward Elgar: Letters of a Lifetime* (Oxford, 1990)
- [Project Gutenberg #345: Dracula](https://www.gutenberg.org/ebooks/345) and [#10: The King James Bible](https://www.gutenberg.org/ebooks/10)
- [eBible.org: Latin Clementine Vulgate](https://ebible.org/find/details.php?id=latVUC)
