import itertools
import math

# --- Your transcribed cipher, as loop-count+direction codes ---
line1 = [23,37,24,33,13,25,11,33,16,22,34,28,18,14,28,35,28,28,25,37,37,28,15,12,22,12,15,34,38]
line2 = [11,28,11,25,11,33,16,26,33,24,28,28,24,25,15,18,11,28,34,28,25,21,38,18,14,12,16,16,12,34,38]
line3 = [24,38,25,21,38,24,14,21,31,15,38,24,28,25,21,18,38,15,33,37,15,38,25,33,14,13,33]

full_sequence = line1 + line2 + line3
assert len(full_sequence) == 87

# --- The three letter groups (each has exactly 8 letters/tokens) ---
GROUP_AH = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
GROUP_IQ = ['IJ', 'K', 'L', 'M', 'N', 'O', 'P', 'Q']
GROUP_RZ = ['R', 'S', 'T', 'UV', 'W', 'X', 'Y', 'Z']
GROUPS = [GROUP_AH, GROUP_IQ, GROUP_RZ]
GROUP_NAMES = {id(GROUP_AH): 'A-H', id(GROUP_IQ): 'I/J-Q', id(GROUP_RZ): 'R-Z'}

# --- Simplified English bigram frequency table (%), rest fall back to FLOOR ---
BIGRAM_FREQ = {
    'TH':3.56,'HE':3.07,'IN':2.43,'ER':2.05,'AN':1.99,'RE':1.85,'ND':1.35,
    'ON':1.31,'EN':1.26,'AT':1.12,'OU':1.12,'ED':1.08,'HA':1.07,'TO':1.06,
    'OR':1.06,'IT':1.05,'IS':1.02,'HI':0.93,'ES':0.93,'NG':0.89,'VE':0.83,
    'CO':0.79,'DE':0.76,'RA':0.75,'AL':0.75,'RO':0.73,'LI':0.70,'IC':0.70,
    'LE':0.69,'AR':0.68,'NT':0.66,'IO':0.66,'ST':0.66,'EA':0.65,'ME':0.63,
    'SE':0.63,'FO':0.62,'TE':0.62,'NE':0.62,'WA':0.60,'CA':0.55,'AS':0.55,
    'LA':0.55,'RI':0.55,'WI':0.55,'BE':0.53,'SI':0.51,'DI':0.51,'EL':0.51,
    'UR':0.50,'CH':0.48,'AC':0.48,'TA':0.47,'EC':0.47,'PE':0.46,'SO':0.45,
    'NO':0.45,'LO':0.44,'FA':0.42,'MA':0.42,'PR':0.42,'IL':0.41,'MO':0.40,
    'GE':0.40,'AM':0.40,'PA':0.39,'SS':0.39,'NS':0.38,'US':0.38,'GA':0.37,
    'CT':0.37,'RS':0.36,'WE':0.36,'UT':0.36,'AD':0.36,'OM':0.36,'ID':0.35,
    'ET':0.35,'IE':0.35,'EE':0.33,'IV':0.33,'NC':0.33,'UN':0.33,'EM':0.33,
    'SP':0.32,'GH':0.31,'TI':0.31,'OW':0.30,'AI':0.30,'PL':0.30,'IG':0.29,
    'OS':0.29,'UL':0.28,'OP':0.28,
}
FLOOR = 0.02  # assumed frequency for any bigram not in the table above


def bigram_score(text):
    """Average log-frequency per adjacent letter pair. Higher = more English-like."""
    text = text.replace('IJ', 'I').replace('UV', 'U')
    total, count = 0.0, 0
    for a, b in zip(text, text[1:]):
        total += math.log(BIGRAM_FREQ.get(a + b, FLOOR))
        count += 1
    return total / count if count else float('-inf')


def build_key(rotation, reverse, group):
    """Return {direction 1-8: letter} for one wheel at a given rotation/direction."""
    g = list(reversed(group)) if reverse else list(group)
    g = g[rotation:] + g[:rotation]
    return {direction: g[direction - 1] for direction in range(1, 9)}


def decode(sequence, wheel_map):
    out = []
    for code in sequence:
        loop, direction = divmod(code, 10)
        out.append(wheel_map[loop][direction])
    return ''.join(out)


def run_bruteforce(independent_reverse=False, top_n=25):
    """
    independent_reverse=False -> your 3,072 x 2 = 6,144 model
        (all three wheels flip direction together, as one global choice)
    independent_reverse=True  -> the 3,072 x 8 = 24,576 model
        (each wheel can independently be forward or reversed)
    """
    results = []
    group_perms = list(itertools.permutations(GROUPS))  # 6 wheel<->group assignments

    if independent_reverse:
        rev_combos = list(itertools.product([False, True], repeat=3))  # 8 combos
    else:
        rev_combos = [(r, r, r) for r in (False, True)]  # 2 combos, all wheels matching

    for group_perm in group_perms:
        for rev_combo in rev_combos:
            for rot1 in range(8):
                for rot2 in range(8):
                    for rot3 in range(8):
                        wheel_map = {
                            1: build_key(rot1, rev_combo[0], group_perm[0]),
                            2: build_key(rot2, rev_combo[1], group_perm[1]),
                            3: build_key(rot3, rev_combo[2], group_perm[2]),
                        }
                        decoded = decode(full_sequence, wheel_map)
                        score = bigram_score(decoded)
                        results.append((score, decoded, group_perm, rev_combo, (rot1, rot2, rot3)))

    results.sort(key=lambda r: r[0], reverse=True)
    print(f"Total combinations evaluated: {len(results)}\n")
    for score, decoded, group_perm, rev_combo, rots in results[:top_n]:
        assign = [GROUP_NAMES[id(g)] for g in group_perm]
        print(f"score={score:.4f}  wheels(1,2,3)={assign}  reverse={rev_combo}  rot={rots}")
        print(f"  {decoded[:29]}")
        print(f"  {decoded[29:60]}")
        print(f"  {decoded[60:]}")
        print()


if __name__ == '__main__':
    # run_bruteforce(independent_reverse=False, top_n=25)  # 6,144 lines
    run_bruteforce(independent_reverse=True, top_n=25)  # 24,576 lines