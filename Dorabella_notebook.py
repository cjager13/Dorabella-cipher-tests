#!/usr/bin/env python3
"""
Dorabella NOTEBOOK-KEY reruns: every "Elgar's key" test again, with the key Elgar actually
wrote in his notebook.

WHY (audit of 16 Nov 2026)
Rounds 1-2 took "Elgar's published key" from the letters BLTACEIARWUNISNF... on the Cipher
Mysteries page. That string comes from a tidy pigpen-style key (ABC DEF GHI JKL on the upright
cross, MNO PQR STU WXY on the diagonal cross, some groups running backwards). It is NOT the key
in Elgar's notebook. The notebook writes the alphabet straight round the eight directions:

    opening right   A B C        opening left      N O P
    down-right      D E F        up-left           Q R S
    down            G H I/J      up                T U/V W
    down-left       K L M        up-right          X Y Z      (1, 2, 3 loops in each group)

and the notebook's own line "A VERY OLD CYPHER" starts with the one-loop symbol opening right,
which is A in this key and C in the other. This is the key printed on the usual picture of the
cipher (Dorabella_cipher.jpg), and it gives the well-known reading BPECAHTCKYFRQDRIRRHPPR...
So the round 2 families built on "Elgar's key" used the wrong base key. This script reruns them.
Tests that search for ANY key are not affected and are not rerun.

  python Dorabella_notebook.py three     --selftest | --stop-when-done | --report-only
  python Dorabella_notebook.py date      (same options)     the date as the key
  python Dorabella_notebook.py scales    (same options)     musical scales as the key
  python Dorabella_notebook.py numbers   (same options)     birthdays, ages, line lengths
  python Dorabella_notebook.py notes     (same options)     note names as the key
  python Dorabella_notebook.py tickshift (same options)     ticks as shifts
  python Dorabella_notebook.py spokes    (same options)     NEW: every walk round the eight spokes
  add  --reading ours  to use this project's own reading of the 5 disputed symbols (see below)

Needs the scripts it reruns (Dorabella_three.py, Dorabella_date.py, Dorabella_scales.py,
Dorabella_tickshift.py, Dorabella_circlekeys.py) and everything they need. Ctrl+C any time; the
same command resumes. Output folders: nb_three_out, nb_date_out, nb_scales_out, nb_numbers_out,
nb_notes_out, nb_tickshift_out, nb_spokes_out  (with "_ours" added for --reading ours).

PRE-REGISTERED (frozen 16 Nov 2026, before any run)
  Key: the notebook key above, in this project's codes: ABC on direction 3, DEF 4, GHI 5, KLM 6,
  NOP 7, QRS 8, TUW 1, XYZ 2; one loop = first letter of the group (J read as I, V as U).
  three      families TE, PA, PB of Dorabella_three.py, unchanged except for the key. Line to
             beat: p < 0.00625 (the original 8-family line of that script).
  date, scales, numbers, notes
             family DE of Dorabella_date.py / Dorabella_scales.py, unchanged except for the key.
             Line: p < 0.0167 (the original 3-family line), and the random-number comparison
             must also be passed, as before.
  tickshift  the investigator's rule exactly as stated (the directions never change; each
             group's ticks point out or in; move 1 or move by the loop count), now on the
             notebook's own layout: 2^8 x 2 = 512 keys. Line: p < 0.05.
  spokes     NEW family SW: the alphabet walked round the eight directions from any starting
             direction (8), clockwise or anticlockwise (2), each group running up or down
             (2^8 = 256; the notebook's right-hand rows show groups written both ways):
             4,096 keys x 5 reading orders, no search. The notebook key is one of them.
             Line: p < 0.05.
  Transcription for all of these: fixed-key tests cannot choose between two readings inside a
  search, so the reading matters. The primary run uses the CONSENSUS reading of Hauer et al.
  (2025, five independent transcriptions), which differs from Dorabella_transcription.csv at
  five symbols: 23 (14 not 15), 32 (11 not 18), 34 (12 not 11), 85 (14 not 15), 86 (12 not 13).
  Four of the five are this project's own alternates; position 32 is not (it was "settled" as
  18 on 1 Oct 2026, but every published reading has 11). "--reading ours" repeats a test with
  the CSV's main readings; that run is secondary and is reported as such.
  Everything else (copies, controls, p-values, power) is exactly as in the original scripts.
"""
import os
import sys

TESTS = ('three', 'date', 'scales', 'numbers', 'notes', 'tickshift', 'spokes')
# The choice is also kept in the environment, so worker processes (which import this file
# again, after the options below have been taken off the command line) make the same choice.
WHICH = next((a for a in sys.argv[1:] if a in TESTS), None) or os.environ.get('DORABELLA_NB_TEST')
READING = None
for _i, _a in enumerate(sys.argv):
    if _a == '--reading' and _i + 1 < len(sys.argv):
        READING = sys.argv[_i + 1]
    elif _a.startswith('--reading='):
        READING = _a.split('=', 1)[1]
READING = READING or os.environ.get('DORABELLA_NB_READING') or 'consensus'
if READING not in ('consensus', 'ours'):
    raise SystemExit("--reading must be 'consensus' or 'ours'")
if WHICH:
    os.environ['DORABELLA_NB_TEST'] = WHICH
os.environ['DORABELLA_NB_READING'] = READING

import Dorabella_three as Z                                     # noqa: E402

NOTEBOOK_SETS = [('ABC', 3, '^'), ('DEF', 4, '^'), ('GHI', 5, '^'), ('KLM', 6, '^'),
                 ('NOP', 7, '^'), ('QRS', 8, '^'), ('TUW', 1, '^'), ('XYZ', 2, '^')]
NOTEBOOK_GROUPS = [s[0] for s in NOTEBOOK_SETS]
# the widely printed notebook-key reading of the cipher (benzedrine.ch, dcode.fr)
NOTEBOOK_READING = ("BPECAHTCKYFRQDRIRRHPPRDXYXGFS" "TRTHTCKLCERREHGQTRFRHUSQDXKKXFS"
                    "ESHUSEDUWGSERHUQSDCPGSHCDXC")
CONSENSUS = {23: 14, 32: 11, 34: 12, 85: 14, 86: 12}            # 1-based position: code
KEY_NOTE = ("KEY: Elgar's NOTEBOOK key (ABC DEF GHI KLM NOP QRS TUW XYZ round directions 3,4,5,6,7,8,1,2; "
            "1 loop = first letter)")
READ_NOTE = ("TRANSCRIPTION: consensus reading (Hauer et al. 2025): symbols 23, 32, 34, 85, 86 read 14, 11, 12, 14, 12"
             if READING == 'consensus' else
             "TRANSCRIPTION: this project's own main readings (secondary run)")

# ---------------------------------------------------------------- applied at import (workers too)
Z.ELGAR_SETS = NOTEBOOK_SETS
Z.ELGAR = Z.build_elgar()
Z.ELGAR_INV = [[s for s in range(24) if Z.ELGAR[s] == L] for L in range(24)]
Z.STANDARD_DECRYPT = NOTEBOOK_READING
assert sorted(Z.ELGAR) == list(range(24)), "the notebook key gives every symbol its own letter"
_orig_load_seq = Z.load_seq


def load_seq(path='Dorabella_transcription.csv'):
    seq, source = _orig_load_seq(path)
    if READING == 'consensus':
        seq = list(seq)
        for pos, code in CONSENSUS.items():
            seq[pos - 1] = code
        source += ' + consensus reading at 23, 32, 34, 85, 86'
    return seq, source


Z.load_seq = load_seq


def _rename(text):
    for a, b in (("Elgar's own published key", "Elgar's notebook key"),
                 ("Elgar's published key", "Elgar's notebook key"),
                 ("Elgar's own key", "Elgar's notebook key"),
                 ("Elgar's key", "Elgar's notebook key"),
                 ("circle-and-cross keys", "spoke-walk keys")):
        text = text.replace(a, b)
    return text


def _wrap_print(mod):
    """Screen messages of the original script, with the key named correctly."""
    def _print(*a, **k):
        print(*[_rename(x) if isinstance(x, str) else x for x in a], **k)
    mod.print = _print


def _retitle(text):
    """Report wording: these runs use the notebook key."""
    head, _, rest = _rename(text).partition('\n')
    return f"{head}   [NOTEBOOK-KEY RERUN]\n{KEY_NOTE}\n{READ_NOTE}\n{rest}"


def _wrap_report(mod):
    _wrap_print(mod)
    old = mod.build_report

    def build_report(*a, **k):
        return _retitle(old(*a, **k))
    mod.build_report = build_report


def _strip_own_args():
    out, skip = [], False
    for a in sys.argv:
        if skip:
            skip = False
            continue
        if a in TESTS:
            continue
        if a == '--reading':
            skip = True
            continue
        if a.startswith('--reading='):
            continue
        out.append(a)
    sys.argv[:] = out


def _default_out(name):
    if not any(a == '--out' or a.startswith('--out=') for a in sys.argv):
        sys.argv += ['--out', f"nb_{name}{'_ours' if READING == 'ours' else ''}_out"]


MOD = None
if WHICH == 'three':
    MOD = Z
    _wrap_report(Z)
elif WHICH in ('date', 'scales', 'numbers', 'notes'):
    import Dorabella_date as D                                  # noqa: E402
    MOD = D
    _wrap_report(D)
elif WHICH == 'tickshift':
    import Dorabella_tickshift as T                             # noqa: E402
    MOD = T
    T.GROUPS = [('notebook groups', NOTEBOOK_GROUPS)]
    T.NKEYS = len(T.GROUPS) * 256 * len(T.STEPS)
    T.FAMILY_DESC = {'TS': f"every tick-shift key on Elgar's notebook key ({T.NKEYS:,} keys x 5 reading orders), "
                           "no search"}
    assert T.PLACES == [(d, ud) for _, d, ud in NOTEBOOK_SETS]
    _wrap_report(T)
elif WHICH == 'spokes':
    import Dorabella_circlekeys as C                            # noqa: E402
    MOD = C
    C.NKEYS = 8 * 2 * 256
    C.DEFAULT_TARGETS = {'CK': dict(shuffle=2000, doublet=5000, controls=100)}
    C.FAMILY_DESC = {'CK': f"SW: every walk round the eight spokes ({C.NKEYS:,} keys x 5 reading orders), no search"}

    def _decode_index(ix):
        ix, bits = divmod(ix, 256)
        start, sense = divmod(ix, 2)
        return start, sense, bits

    def _key_parts(ix):
        start, sense, bits = _decode_index(ix)
        step = 1 if sense == 0 else -1
        return [(letters, (start + step * k) % 8 + 1, '^' if bits >> k & 1 else 'v')
                for k, letters in enumerate(NOTEBOOK_GROUPS)]

    def _label(ix):
        start, sense, bits = _decode_index(ix)
        lay = ' '.join(f"{l}@{d}{u}" for l, d, u in _key_parts(ix))
        tag = "  <- ELGAR'S NOTEBOOK KEY" if C.key_of(ix) == Z.ELGAR else ''
        return f"from direction {start + 1}, {'clockwise' if sense == 0 else 'anticlockwise'} | {lay}{tag}"

    C.decode_index, C.key_parts, C.label = _decode_index, _key_parts, _label
    _wrap_print(C)
    _old_ck_report = C.build_report

    def _ck_report(*a, **k):
        t = _old_ck_report(*a, **k)
        for x, y in (("CIRCLE-AND-CROSS KEYS REPORT  (every 'eight sets of three' picture as the key)",
                      "SPOKE-WALK KEYS REPORT  (family SW: the alphabet walked round the eight directions)"),
                     ("circle-and-cross keys", "spoke-walk keys"), ("circle keys", "spoke-walk keys")):
            t = t.replace(x, y)
        return _retitle(t)
    C.build_report = _ck_report


def main():
    if WHICH is None:
        raise SystemExit("choose a test:  python Dorabella_notebook.py " + ' | '.join(TESTS) +
                         "   [--selftest | --stop-when-done | --report-only] [--reading ours]")
    _strip_own_args()
    print(KEY_NOTE)
    print(READ_NOTE)
    if WHICH == 'three':
        _default_out('three')
        if '--stop-when-done' in sys.argv:              # Dorabella_three.py has no such option: added here
            sys.argv.remove('--stop-when-done')
            _old_stage_list = Z.stage_list

            def stage_list(state, only):
                out = _old_stage_list(state, only)
                if not out:
                    raise SystemExit("all targets met; stopping (--stop-when-done). Progress and report are saved.")
                return out
            Z.stage_list = stage_list
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'TE', 'PA', 'PB']
        Z.main()
    elif WHICH == 'date':
        _default_out('date')
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'DE']
        MOD.main()
    elif WHICH in ('scales', 'numbers', 'notes'):
        import Dorabella_scales as SC
        _default_out(WHICH)
        if '--only' not in sys.argv:
            sys.argv += ['--only', 'DE']
        SC.configure(WHICH)
        print(f"key set '{WHICH}': {len(MOD.DIGIT_KEYS)} digit keys, {len(MOD.NUMBER_KEYS)} number keys")
        MOD.main()
    elif WHICH == 'tickshift':
        _default_out('tickshift')
        MOD.main()
    elif WHICH == 'spokes':
        _default_out('spokes')
        MOD.main()


if __name__ == '__main__':
    Z.mp.freeze_support()
    main()
