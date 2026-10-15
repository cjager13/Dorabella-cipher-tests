#!/usr/bin/env python3
"""
Dorabella MORE CRIBS: religious words, names of God and Jesus, and Elgar's and Dora's
streets and addresses, through exactly the same crib test as Dorabella_cribs.py.

Needs Dorabella_cribs.py and everything it needs, in the same folder.

  python Dorabella_morecribs.py --selftest          # eligible cribs, one planted control, timing
  python Dorabella_morecribs.py --stop-when-done    # results in morecribs_out\\
  python Dorabella_morecribs.py --report-only

Ctrl+C any time; the same command resumes.

THE IDEA
Elgar was a practising Catholic and church organist, Dora was a rector's daughter, and
the 1920 notebook page carries a circle with a cross. A note between them might name God
or Jesus, quote a prayer, or name a church, a street or a house. A crib is placed wherever
its pattern of repeated letters fits the cipher, its letters are fixed, the rest of the
key is searched, and the result is scored; scrambled copies get the identical search.

PRE-REGISTERED (frozen 12 Oct 2026, before any run)
  The crib list below, all of it. Tests C1 (any key), C2 (Elgar-style sets-of-three keys)
  and C3 (loop-group keys) exactly as in Dorabella_cribs.py; three tests, so the line to
  beat is p < 0.05/3 = 0.0167. The "MY DEAR DORA" side report of that script is left out.
  Eligibility rule, unchanged and decided on scrambled copies only: a crib is used if it
  has at least 5 letters and fits, on average, no more than 20 places in a pair-preserving
  copy. Short words (GOD, AMEN, LORD...) fit almost anywhere and are dropped by the rule.
  Monograms are tested spelled out (IHS as IESUS HOMINUM SALVATOR, INRI, AMDG...), since
  three or four letters cannot narrow anything down. Places come from the Elgar Society's
  list of Elgar's houses up to 1897 and from Dora's home at the Wolverhampton rectory;
  later homes (Birchwood 1898, Craeg Lea 1899) are left out.
  Single words: every word of 5 or more letters inside a multi-word crib is also a crib of
  its own (IESUS HOMINUM SALVATOR adds HOMINUM and SALVATOR), unless already listed; the
  eligibility rule then applies to them as to every other crib. The same rule is applied to
  the multi-word cribs of Dorabella_cribs.py (YOURS SINCERELY adds SINCERELY, STEUART POWELL
  adds STEUART and POWELL...), since that test only used them whole.
"""
import sys

import Dorabella_cribs as C

CRIBS = {
    'names of God': """GOD | THE LORD | OUR LORD | LORD GOD | THE LORD GOD | GOD ALMIGHTY |
        ALMIGHTY GOD | THE ALMIGHTY | ALMIGHTY | JEHOVAH | HEAVENLY FATHER | OUR FATHER |
        GOD THE FATHER | THE CREATOR | PROVIDENCE | DEUS | DOMINUS | DOMINE | DEO""",
    'names of Jesus': """JESUS | JESU | CHRIST | JESUS CHRIST | CHRIST JESUS | OUR LORD JESUS CHRIST |
        OUR SAVIOUR | THE SAVIOUR | SAVIOUR | REDEEMER | THE REDEEMER | EMMANUEL | IMMANUEL |
        MESSIAH | LAMB OF GOD | SON OF GOD | SON OF MAN | GOOD SHEPHERD | LIGHT OF THE WORLD |
        PRINCE OF PEACE | KING OF KINGS | JESUS OF NAZARETH | CHRISTUS""",
    'the Spirit, Mary and the saints': """HOLY GHOST | HOLY SPIRIT | TRINITY | HOLY TRINITY |
        BLESSED TRINITY | VIRGIN MARY | BLESSED VIRGIN | OUR LADY | MOTHER OF GOD | HOLY MOTHER |
        HOLY FAMILY | SAINT GEORGE | SAINT PETER | SAINT CECILIA | SAINT JOSEPH | ANGELS""",
    'monograms spelled out': """IESUS HOMINUM SALVATOR | JESUS SAVIOUR OF MEN | IN HOC SIGNO |
        IN HOC SIGNO VINCES | IESUS NAZARENUS REX IUDAEORUM | JESUS OF NAZARETH KING OF THE JEWS |
        AD MAIOREM DEI GLORIAM | TO THE GREATER GLORY OF GOD | DEO GRATIAS | LAUS DEO |
        SOLI DEO GLORIA | DEO VOLENTE""",
    'prayers and liturgy': """AMEN | ALLELUIA | HALLELUJAH | HOSANNA | HOSANNA IN THE HIGHEST |
        GLORIA | GLORIA IN EXCELSIS | GLORIA IN EXCELSIS DEO | GLORY BE TO GOD | GLORY TO GOD |
        GLORY TO GOD IN THE HIGHEST | GLORY BE TO THE FATHER | PRAISE THE LORD | PRAISE BE TO GOD |
        THANKS BE TO GOD | THANK GOD | GOD BLESS YOU | GOD BLESS | GODSPEED | GOD WILLING |
        GOD SAVE THE QUEEN | OUR FATHER WHICH ART IN HEAVEN | HALLOWED BE THY NAME |
        THY KINGDOM COME | THY WILL BE DONE | GIVE US THIS DAY | DELIVER US FROM EVIL |
        FOR EVER AND EVER | WORLD WITHOUT END | THE LORD BE WITH YOU | AND WITH THY SPIRIT |
        LIFT UP YOUR HEARTS | LORD HAVE MERCY | CHRIST HAVE MERCY | KYRIE ELEISON |
        CHRISTE ELEISON | AGNUS DEI | DONA NOBIS PACEM | SANCTUS | BENEDICTUS | PAX VOBISCUM |
        DOMINUS VOBISCUM | ET CUM SPIRITU TUO | REQUIESCAT IN PACE | REST IN PEACE | TE DEUM |
        MAGNIFICAT | NUNC DIMITTIS | MISERERE | DE PROFUNDIS | AVE MARIA | HAIL MARY |
        FULL OF GRACE | SALVE REGINA | TANTUM ERGO | O SALUTARIS | PRAISE TO THE HOLIEST""",
    "Elgar's sacred music and texts to 1897": """AVE VERUM CORPUS | AVE MARIS STELLA |
        ECCE SACERDOS MAGNUS | O SALUTARIS HOSTIA | TE DEUM AND BENEDICTUS | THE LIGHT OF LIFE |
        LUX CHRISTI | THE BANNER OF SAINT GEORGE | GERONTIUS | THE DREAM OF GERONTIUS""",
    'churches, clergy and services': """SAINT GEORGES | SAINT GEORGES CHURCH | SANSOME PLACE |
        WORCESTER CATHEDRAL | THREE CHOIRS | THREE CHOIRS FESTIVAL | HEREFORD | GLOUCESTER |
        SAINT PETERS | SAINT PETERS CHURCH | COLLEGIATE CHURCH | THE RECTOR | REVEREND |
        CHURCH | CHAPEL | CATHEDRAL | SUNDAY | SERVICE | EVENSONG | MATINS | VESPERS | ORGAN |
        ORGANIST | CHOIR | ANTHEM | HYMN | PSALM | PRAYER | BLESSING | SERMON""",
    "Elgar's homes to 1897": """THE FIRS | BROADHEATH | LOWER BROADHEATH | EDGAR STREET |
        COLLEGE PRECINCTS | HIGH STREET | TEN HIGH STREET | HIGH STREET WORCESTER |
        ELGAR BROTHERS | LORETTA VILLA | CHESTNUT WALK | FIELD TERRACE | BATH ROAD |
        MARLOES ROAD | WEST KENSINGTON | KENSINGTON | SAETERMO | THE LEES | GREAT MALVERN |
        OAKLANDS | FOUNTAIN ROAD | UPPER NORWOOD | AVONMORE ROAD | FORLI | ALEXANDRA ROAD |
        MALVERN LINK | FORLI MALVERN LINK | FORLI ALEXANDRA ROAD""",
    "Dora's home and Wolverhampton": """WATERLOO ROAD | THE RECTORY | RECTORY WATERLOO ROAD |
        MOLINEUX | PAUL STREET | WOLVERHAMPTON | STAFFORDSHIRE | WANDERERS |
        WOLVERHAMPTON WANDERERS | FOOTBALL""",
    'the Malvern hills': """BRITISH CAMP | HEREFORDSHIRE BEACON | WORCESTERSHIRE BEACON |
        MALVERN HILLS | MALVERN WELLS | WELLS ROAD | NORTH HILL | PRIORY | MALVERN PRIORY""",
}

def _single_words(cribs):
    seen = {' '.join(x.split()) for txt in cribs.values() for x in txt.split('|')}
    out = []
    for txt in cribs.values():
        for item in txt.split('|'):
            words = item.split()
            if len(words) < 2:
                continue
            for w in words:
                if len(w) >= 5 and w not in seen:
                    seen.add(w)
                    out.append(w)
    return ' | '.join(out)


CRIBS['single words from the phrases above'] = _single_words(CRIBS)
CRIBS['single words from the first crib list'] = _single_words(C.CRIBS)

# applied at import, so the worker processes (which import this module) get them too
C.CRIBS = CRIBS
C.ALL_CRIBS = C.crib_list()
C.H1_CRIBS = ()
_orig_report = C.build_report
_orig_stages = C.stages
STOP_WHEN_DONE = '--stop-when-done' in sys.argv
if STOP_WHEN_DONE:
    sys.argv.remove('--stop-when-done')


def build_report(state, sd):
    t = _orig_report(state, sd)
    t = t.replace('DORABELLA CRIB TEST REPORT',
                  'DORABELLA MORE-CRIBS REPORT  (religious words, names of God and Jesus, streets and addresses)')
    t = t.replace('; H1 is exploratory whatever its p.', '.')
    return '\n'.join(line for line in t.split('\n')
                     if not line.startswith('H1 ') and not line.startswith('  H1 ')) + '\n'


def stages(state):
    out = _orig_stages(state)
    if not out and STOP_WHEN_DONE:
        raise SystemExit("all targets met; stopping (--stop-when-done). Progress and report are saved.")
    return out


C.build_report = build_report
C.stages = stages

if __name__ == '__main__':
    C.mp.freeze_support()
    if not any(a.startswith('--out') for a in sys.argv):
        sys.argv += ['--out', 'morecribs_out']
    C.main()
