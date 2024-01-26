import common
from common import rotl, parity

def correlation_sign_and(a, b, c, u, v, w):

    def complement(t):
        return 0xffff ^ t

    p = parity(((complement(a) & u) ^ (c & v)) & ((complement(b) & v) ^ (c & u))) +\
        parity((u & v) & (c ^ (a & b & complement(c))))
    return (-1) ** p

def correlation_sign(masks, differences):
    s = 1
    for i in range(len(masks)):
        s *= correlation_sign_and(*differences[i], *masks[i])
    return s

diffs = [
    (0x8000, 0x0003),
    (0x0000, 0x8000),
    (0x8000, 0x0000),
    (0x0082, 0x8000),
    (0x820c, 0x0082),
    (0x08a8, 0x820c),
    (0xa1ac, 0x08a8)
]

word_size = 16
max_weight = 7

and_differences = []
for i in range(len(diffs) - 1):
    (x, y) = diffs[i]
    z = diffs[i + 1][0]
    and_differences.append(
        (rotl(x, 1, word_size), rotl(x, 8, word_size), rotl(x, 2, word_size) ^ y ^ z)
    )

for (a, b) in diffs:
    print(" " * 4, "{:04x} {:04x}".format(a, b))

groups = dict()
model, weight = common.simon_quasidifferential_trails(diffs, word_size)
for target in range(max_weight):
    print("Solving for weight ", target)
    for solution in common.solve_all_simon(model, weight, target, len(diffs) - 1, word_size):
        trail = []
        for i in range(len(solution) - 1):
            (u, v, w)  = solution[i]
            (_, v_, _) = solution[i + 1]

            x = rotl(u ^ w, 1, word_size)
            y = rotl(rotl(v, word_size - 2, word_size) ^ w ^ v_, 8, word_size)
            trail.append((x, y, v))

        trail_key = tuple([x for (_, _, x) in trail])
        if trail_key not in groups:
            groups[trail_key] = []
        groups[trail_key].append((target, trail))

clusters = dict()
for (trail_key, group) in groups.items():
    print("Group", " ".join(["{:04x}".format(k) for k in trail_key]))
    c = 0
    for (loss, trail) in group:
        c += correlation_sign(trail, and_differences) / 2. ** loss
    if c != 0:
        for (loss, trail) in group:
            print("Loss: {}".format(loss))                           
            for (v, w, x) in trail:
                print(" " * 4, "{:04x} {:04x} {:04x}".format(v, w, x))
        print("Group correlation (relative): ", c)
        print()
        abs_c = abs(c)
        if abs_c not in clusters:
            clusters[abs_c] = []
        clusters[abs_c].append((c / abs_c, trail_key))

for (abs_weight, groups) in clusters.items():
    print("Group correlation (relative): {} [{} groups]".format(abs_weight, len(groups)))
    for (sign, trail_key) in groups:
        s_trail = " ".join("{:04x}".format(t) for t in trail_key)
        print(" " * 4, "{:+}".format(sign), s_trail)
