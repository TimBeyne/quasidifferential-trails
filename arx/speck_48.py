import common

# Weight 46 characteristic
diffs = [
    (0x504200, 0x004240),
    (0x001202, 0x020002),
    (0x000010, 0x100000),
    (0x000000, 0x800000),
    (0x800000, 0x800004),
    (0x808004, 0x808020),
    (0x8400A0, 0x8001A4),
    (0x608DA4, 0x608080),
    (0x042003, 0x002400),
    (0x012020, 0x000020),
    (0x200100, 0x200000),
    (0x202001, 0x202000)
]

# Weight 47 characteristic
#diffs = [
#    (0x504200, 0x004240),
#    (0x001202, 0x020002),
#    (0x000010, 0x100000),
#    (0x000000, 0x800000),
#    (0x800000, 0x800004),
#    (0x808004, 0x808020),
#    (0x8400A0, 0x8001A4),
#    (0xE08DA4, 0xE08080),
#    (0x042007, 0x002400),
#    (0x012020, 0x000020),
#    (0x200100, 0x200000),
#    (0x202001, 0x202000)
#]

word_size = 24
max_weight_loss = 4

sols = []
btor, weight = common.speck_quasidifferential_trails(diffs, word_size)
for w in range(max_weight_loss):
    for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size):
        sols.append((w, solution))

print("Characteristic:")
for (a, b) in diffs:
    print(" " * 4, "{:06x} {:06x}".format(a, b))

for (weight, solution) in sols:
    s = common.compute_sign_speck(diffs, solution, word_size)
    print("Weight: {} [{:+}]".format(weight, s))
    for (u, v) in solution:
        print(" " * 4, "{:06x} {:06x}".format(u, v))
