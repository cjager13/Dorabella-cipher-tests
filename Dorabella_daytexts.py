#!/usr/bin/env python3
"""
Build the "texts of the day" for Wednesday 14 July 1897 (the date Elgar wrote on the cipher)
as a source folder for Dorabella_quotes.py.

  python Dorabella_daytexts.py --kjv Bible_KJV\\Bible_KJV.txt
  python Dorabella_quotes.py --source Day_14July1897 --out day_out --selftest
  python Dorabella_quotes.py --source Day_14July1897 --out day_out --stop-when-done

PRE-REGISTERED (chosen from the date alone, before any result on these texts):
  Psalms of the day (Book of Common Prayer, 14th day of the month), Coverdale:
      Morning 71, 72   Evening 73, 74
  Psalms from the date's numbers 7 14 97, plus 98, Coverdale: 7, 14, 97, 98
  Lessons for 14 July, 1871 table (in force in 1897), King James Version:
      Morning  Job 38:39-39:30   Acts 18:1-24
      Evening  Job 40            Matthew 6:19-7:7
  Lessons for 14 July, 1662 table (in case the church still followed it):
      Morning  Ecclesiastes 7    John 2
      Evening  Ecclesiastes 8    2 Thessalonians 1
  The week's Collect, Epistle and Gospel (4th Sunday after Trinity, 11 July 1897):
      Collect "O God, the protector of all that trust in thee"   Romans 8:18-23   Luke 6:36-42

Each passage is written to its own file, so no phrase runs from one passage into another.
The Psalms are built in (from the cleaned Coverdale Psalter); the Bible lessons are cut
from your own KJV file (Project Gutenberg #10 layout: verses start "chapter:verse").
"""
import argparse
import os
import re

PSALMS = {
    7: ('Domine, Deus meus',
         'O Lord my God, in thee have I put my trust; save me from all them that persecute me, and deliver me, lest he devour my soul like a lion, and tear it in pieces, while there is none to help. O Lord my God, if I have done any such thing, or if there be any wickedness in my hands, if I have rewarded evil unto him that dealt friendly with me (yea, I have delivered him that without any cause is mine enemy), then let mine enemy persecute my soul, and take me; yea, let him tread my life down upon the earth, and lay mine honour in the dust. Stand up, O Lord, in thy wrath, and lift up thyself, because of the indignation of mine enemies; arise up for me in the judgement that thou hast commanded. And so shall the congregation of the people come about thee; for their sakes therefore lift up thyself again. The Lord shall judge the people; give sentence with me, O Lord, according to my righteousness, and according to the innocency that is in me. O let the wickedness of the ungodly come to an end; but guide thou the just. For the righteous God trieth the very hearts and reins. My help cometh of God, who preserveth them that are true of heart. God is a righteous Judge, strong, and patient; and God is provoked every day. If a man will not turn, he will whet his sword; he hath bent his bow, and made it ready. He hath prepared for him the instruments of death; he ordaineth his arrows against the persecutors. Behold, he travaileth with iniquity; he hath conceived sorrow, and brought forth ungodliness. He hath graven and dug up a pit, and is fallen himself into the destruction that he made for others. For his travail shall come upon his own head, and his wickedness shall fall on his own pate. I will give thanks unto the Lord, according to his righteousness; and I will praise the Name of the Lord Most High.'),
    14: ('Dixit insipiens',
         'The fool hath said in his heart, There is no God. They are corrupt, and become abominable in their doings; there is none that doeth good, no not one. The Lord looked down from heaven upon the children of men, to see if there were any that would understand, and seek after God. But they are all gone out of the way, they are altogether become abominable; there is none that doeth good, no not one. Their throat is an open sepulchre. With their tongues have they deceived; the poison of asps is under their lips. Their mouth is full of cursing and bitterness; their feet are swift to shed blood. Destruction and unhappiness is in their ways, and the way of peace have they not known; there is no fear of God before their eyes. Have they no knowledge, that they are all such workers of mischief, eating up my people as it were bread, and call not upon the Lord? There were they brought in great fear, even where no fear was; for God is in the generation of the righteous. As for you, ye have made a mock at the counsel of the poor, because he putteth his trust in the Lord. Who shall give salvation unto Israel out of Sion? When the Lord turneth the captivity of his people, then shall Jacob rejoice, and Israel shall be glad.'),
    71: ('In te, Domine, speravi',
         'In THEE, O Lord, have I put my trust; let me never be put to confusion, but rid me, and deliver me, in thy righteousness; incline thine ear unto me, and save me. Be thou my stronghold, whereunto I may always resort; thou hast promised to help me, for thou art my house of defence, and my castle. Deliver me, O my God, out of the hand of the ungodly, out of the hand of the unrighteous and cruel man. For thou, O Lord God, art the thing that I long for; thou art my hope, even from my youth. Through thee have I been holden up ever since I was born; thou art he that took me out of my mother’s womb. My praise shall be always of thee. I am become as it were a monster unto many, but my sure trust is in thee. O let my mouth be filled with thy praise, that I may sing of thy glory and honour all the day long. Cast me not away in the time of age; forsake me not when my strength faileth me. For mine enemies speak against me; and they that lay wait for my soul take their counsel together, saying, God hath forsaken him. Persecute him, and take him, for there is none to deliver him. Go not far from me, O God; my God, haste thee to help me. Let them be confounded and perish that are against my soul; let them be covered with shame and dishonour that seek to do me evil. As for me, I will patiently abide always, and will praise thee more and more. My mouth shall daily speak of thy righteousness and salvation; for I know no end thereof. I will go forth in the strength of the Lord God, and will make mention of thy righteousness only. Thou, O God, hast taught me from my youth up until now; therefore will I tell of thy wondrous works. Forsake me not, O God, in mine old age, when I am gray-headed, until I have showed thy strength unto this generation, and thy power to all them that are yet for to come. Thy righteousness, O God, is very high, and great things are they that thou hast done; O God, who is like unto thee? O what great troubles and adversities hast thou showed me! And yet didst thou turn and refresh me, yea, and broughtest me from the deep of the earth again. Thou hast brought me to great honour, and comforted me on every side. Therefore will I praise thee, and thy faithfulness, O God, playing upon an instrument of music; unto thee will I sing upon the harp, O thou Holy One of Israel. My lips will be fain when I sing unto thee; and so will my soul whom thou hast delivered. My tongue also shall talk of thy righteousness all the day long; for they are confounded and brought unto shame that seek to do me evil.'),
    72: ('Deus, judicium',
         'GIVE THE KING thy judgements, O God, and thy righteousness unto the King’s son. Then shall he judge thy people according unto right, and defend the poor. The mountains also shall bring peace, and the little hills righteousness unto the people. He shall keep the simple folk by their right, defend the children of the poor, and punish the wrong doer. They shall fear thee, as long as the sun and moon endureth, from one generation to another. He shall come down like the rain into a fleece of wool, even as the drops that water the earth. In his time shall the righteous flourish, yea, and abundance of peace, so long as the moon endureth. His dominion shall be also from the one sea to the other, and from the flood unto the world’s end. They that dwell in the wilderness shall kneel before him; his enemies shall lick the dust. The kings of Tarshish and of the isles shall give presents; the kings of Arabia and Saba shall bring gifts. All kings shall fall down before him; all nations shall do him service. For he shall deliver the poor when he crieth, the needy also, and him that hath no helper. He shall be favourable to the simple and needy, and shall preserve the souls of the poor. He shall deliver their souls from falsehood and wrong; and dear shall their blood be in his sight. He shall live, and unto him shall be given of the gold of Arabia; prayer shall be made ever unto him, and daily shall he be praised. There shall be an heap of corn in the earth, high upon the hills; his fruit shall shake like Lebanon, and shall be green in the city like grass upon the earth. His Name shall endure for ever; his Name shall remain under the sun among the posterities, which shall be blessed in him. And all the heathen shall praise him. Blessed be the Lord God, even the God of Israel, which only doeth wondrous things; and blessed be the Name of his Majesty for ever. And all the earth shall be filled with his Majesty. Amen, Amen.'),
    73: ('Quam bonus Israel!',
         'Truly God is loving unto Israel, even unto such as are of a clean heart. Nevertheless, my feet were almost gone; my treadings had well-nigh slipped. And why? I was grieved at the wicked; I do also see the ungodly in such prosperity. For they are in no peril of death, but are lusty and strong. They come in no misfortune like other folk; neither are they plagued like other men. And this is the cause that they are so holden with pride, and overwhelmed with cruelty. Their eyes swell with fatness, and they do even what they lust. They corrupt others, and speak of wicked blasphemy; their talking is against the Most High. For they stretch forth their mouth unto the heaven, and their tongue goeth through the world. Therefore fall the people unto them, and thereout suck they no small advantage. Tush, say they, how should God perceive it? Is there knowledge in the Most High? Lo, these are the ungodly, these prosper in the world, and these have riches in possession; and I said, Then have I cleansed my heart in vain, and washed my hands in innocency. All the day long have I been punished, and chastened every morning. Yea, and I had almost said even as they; but lo, then I should have condemned the generation of thy children. Then thought I to understand this; but it was too hard for me, until I went into the sanctuary of God. Then understood I the end of these men, namely, how thou dost set them in slippery places, and castest them down, and destroyest them. O how suddenly do they consume, perish, and come to a fearful end! Yea, even like as a dream when one awaketh, so shalt thou make their image to vanish out of the city. Thus my heart was grieved, and it went even through my reins. So foolish was I, and ignorant, even as it were a beast before thee. Nevertheless, [am always by thee; for thou hast holden me by my right hand. Thou shalt guide me with thy counsel, and after that receive me with glory. Whom have I in heaven but thee? And there is none upon earth that I desire in comparison of thee. My flesh and my heart faileth; but God is the strength of my heart, and my portion for ever. For lo, they that forsake thee shall perish; thou hast destroyed all them that commit fornication against thee. But it is good for me to hold me fast by God, to put my trust in the Lord God, and to speak of all thy works in the gates of the daughter of Sion.'),
    74: ('Ut quid, Deus?',
         'O God, wherefore art thou absent from us so long? Why is thy wrath so hot against the sheep of thy pasture? O think upon thy congregation, whom thou hast purchased, and redeemed of old. Think upon the tribe of thine inheritance, and Mount Sion, wherein thou hast dwelt. Lift up thy feet, that thou mayest utterly destroy every enemy which hath done evil in thy sanctuary. Thine adversaries roar in the midst of thy congregations, and set up their banners for tokens. He that hewed timber afore out of the thick trees was known to bring it to an excellent work. But now they break down all the carved work thereof with axes and hammers. They have set fire upon thy holy places, and have defiled the dwelling-place of thy Name, even unto the ground. Yea, they said in their hearts, let us make havoc of them altogether. Thus have they burnt up all the houses of God in the land. We see not our tokens. There is not one prophet more; no, not one is there among us, that understandeth any more. O God, how long shall the adversary do this dishonour? Shall the enemy blaspheme thy Name for ever? Why withdrawest thou thy hand? Why pluckest thou not thy right hand out of thy bosom to consume the enemy? For God is my King of old; the help that is done upon earth, he doeth it himself. Thou didst divide the sea through thy power; thou brokest the heads of the dragons in the waters. Thou smotest the heads of Leviathan in pieces, and gavest him to be meat for the people of the wilderness. Thou broughtest out fountains and waters out of the hard rocks; thou driedst up mighty waters. The day is thine, and the night is thine; thou hast prepared the light and the sun. Thou hast set all the borders of the earth; thou hast made summer and winter. Remember this, O Lord, how the enemy hath rebuked, and how the foolish people hath blasphemed thy Name. O deliver not the soul of thy turtle-dove unto the multitude of the enemies; and forget not the congregation of the poor for ever. Look upon the covenant; for all the earth is full of darkness and cruel habitations. O let not the simple go away ashamed; but let the poor and needy give praise unto thy Name. Arise, O God, maintain thine own cause; remember how the foolish man blasphemeth thee daily. Forget not the voice of thine enemies; the presumption of them that hate thee increaseth ever more and more.'),
    97: ('Dominus regnavit',
         'The Lord is king. The earth may be glad thereof; yea, the multitude of the isles may be glad thereof. Clouds and darkness are round about him; righteousness and judgement are the habitation of his seat. There shall go a fire before him, and burn up his enemies on every side. His lightnings gave shine unto the world; the earth saw it, and was afraid. The hills melted like wax at the presence of the Lord, at the presence of the Lord of the whole earth. The heavens have declared his righteousness, and all the people have seen his glory. Confounded be all they that worship carved images, and that delight in vain gods; worship him, all ye gods. Sion heard of it, and rejoiced; and the daughters of Judah were glad, because of thy judgements, O Lord. For thou, Lord, art higher than all that are in the earth; thou art exalted far above all gods. O ye that love the Lord, see that ye hate the thing which is evil. The Lord preserveth the souls of his saints; he shall deliver them from the hand of the ungodly. There is sprung up a light for the righteous, and joyful gladness for such as are true-hearted. Rejoice in the Lord, ye righteous, and give thanks for a remembrance of his holiness.'),
    98: ('Cantate Domino',
         'O sing unto the Lord a new song; for he hath done marvellous things. With his own right hand, and with his holy arm, hath he gotten himself the victory. The Lord declared his salvation; his righteousness hath he openly showed in the sight of the heathen. He hath remembered his mercy and truth toward the house of Israel; and all the ends of the world have seen the salvation of our God. Show yourselves joyful unto the Lord, all ye lands; sing, rejoice, and give thanks. Praise the Lord upon the harp; sing to the harp with a psalm of thanksgiving. With trumpets also and shawms, O show yourselves joyful before the Lord the King. Let the sea make a noise, and all that therein is, the round world, and they that dwell therein. Let the floods clap their hands, and let the hills be joyful together before the Lord; for he is come to judge the earth. With righteousness shall he judge the world, and the people with equity.'),
}

COLLECT_TRINITY4 = ("O God, the protector of all that trust in thee, without whom nothing is strong, "
                    "nothing is holy: Increase and multiply upon us thy mercy; that, thou being our "
                    "ruler and guide, we may so pass through things temporal, that we finally lose not "
                    "the things eternal: Grant this, O heavenly Father, for Jesus Christ's sake our Lord. Amen.")

# canonical book number (1-based), with a phrase from the start of that book to check the count
BOOKS = {'Job': (18, 1, 1, 'man in the land of Uz'),
         'Ecclesiastes': (21, 1, 1, 'words of the Preacher'),
         'Matthew': (40, 1, 1, 'book of the generation'),
         'Luke': (42, 1, 1, 'Forasmuch as many'),
         'John': (43, 1, 1, 'In the beginning was the Word'),
         'Acts': (44, 1, 1, 'former treatise'),
         'Romans': (45, 1, 1, 'Paul, a servant of Jesus Christ'),
         '2 Thessalonians': (53, 1, 3, 'bound to thank God always for you')}

# (file name, book, [(chapter, first verse, last verse or None = to the end of the chapter)])
LESSONS = [
    ('1871_morning_1_Job_38v39-39', 'Job', [(38, 39, None), (39, 1, None)]),
    ('1871_morning_2_Acts_18v1-24', 'Acts', [(18, 1, 24)]),
    ('1871_evening_1_Job_40', 'Job', [(40, 1, None)]),
    ('1871_evening_2_Matthew_6v19-7v7', 'Matthew', [(6, 19, None), (7, 1, 7)]),
    ('1662_morning_1_Ecclesiastes_7', 'Ecclesiastes', [(7, 1, None)]),
    ('1662_morning_2_John_2', 'John', [(2, 1, None)]),
    ('1662_evening_1_Ecclesiastes_8', 'Ecclesiastes', [(8, 1, None)]),
    ('1662_evening_2_2Thessalonians_1', '2 Thessalonians', [(1, 1, None)]),
    ('trinity4_epistle_Romans_8v18-23', 'Romans', [(8, 18, 23)]),
    ('trinity4_gospel_Luke_6v36-42', 'Luke', [(6, 36, 42)]),
]


def parse_kjv(path):
    """Return a list of books, each a dict {(chapter, verse): text}, in file order."""
    with open(path, encoding='utf-8', errors='ignore') as fh:
        raw = fh.read()
    s = raw.find('*** START OF')
    if s >= 0:
        raw = raw[raw.find('\n', s) + 1:]
    e = raw.find('*** END OF')
    if e >= 0:
        raw = raw[:e]
    marks = list(re.finditer(r'(?<![\w:])(\d+):(\d+)\s', raw))
    books, cur = [], None
    for k, m in enumerate(marks):
        ch, vs = int(m.group(1)), int(m.group(2))
        end = marks[k + 1].start() if k + 1 < len(marks) else len(raw)
        text = ' '.join(raw[m.end():end].split())
        if (ch, vs) == (1, 1):
            cur = {}
            books.append(cur)
        if cur is None:
            continue
        cur[(ch, vs)] = text
    return books


def strip_book_title(text):
    """The last verse of a book is followed by the next book's title in the file; the next
    book starts at 1:1, so a trailing title can only follow a chapter's last verse. Cut a
    trailing run of title-like words (e.g. 'The Book of Psalms')."""
    return re.sub(r'\s+(The (First|Second|Third)?\s*(Book|Epistle|Gospel|Acts|Lamentations|Revelation|'
                  r'Proverbs|Song|Ecclesiastes)\b.*)$', '', text)


def passage(book, ranges):
    out = []
    for ch, a, b in ranges:
        v = a
        while (ch, v) in book and (b is None or v <= b):
            out.append(strip_book_title(book[(ch, v)]))
            v += 1
        if v == a:
            raise SystemExit(f"verse {ch}:{a} not found")
    return ' '.join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kjv', default=os.path.join('Bible_KJV', 'Bible_KJV.txt'))
    ap.add_argument('--out', default='Day_14July1897')
    args = ap.parse_args()
    books = parse_kjv(args.kjv)
    print(f"{args.kjv}: {len(books)} books found (expected 66)")
    ok = True
    for name, (num, ch, vs, sig) in BOOKS.items():
        got = books[num - 1].get((ch, vs), '') if num <= len(books) else ''
        good = sig.lower() in got.lower()
        ok &= good
        print(f"  check {name:<16} book {num:>2}, {ch}:{vs} {'OK' if good else 'MISMATCH'}: {got[:60]}")
    if not ok:
        raise SystemExit("The KJV file is not laid out as expected. Send me its first 40 lines and "
                         "a few lines from the start of Job.")
    os.makedirs(args.out, exist_ok=True)
    files = []
    for n, (title, text) in PSALMS.items():
        grp = 'day14' if n in (71, 72, 73, 74) else 'date_numbers'
        files.append((f"psalm_{n:03d}_{grp}_coverdale", f"{title}\n{text}"))
    for fn, book, rng in LESSONS:
        files.append((fn, passage(books[BOOKS[book][0] - 1], rng)))
    files.append(('trinity4_collect', COLLECT_TRINITY4))
    for fn, text in files:
        with open(os.path.join(args.out, fn + '.txt'), 'w', encoding='utf-8') as fh:
            fh.write(text + '\n')
        letters = sum(ch.isalpha() for ch in text)
        print(f"  {fn:<38} {letters:>6} letters   {text[:50]}... {text[-40:]}")
    print(f"{len(files)} files written to {args.out}")


if __name__ == '__main__':
    main()
