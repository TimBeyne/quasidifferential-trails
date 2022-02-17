import common

# 15 round characteristic
diffs = [
    (0x04092400, 0x20040104),
    (0x20000820, 0x20200001),
    (0x00000009, 0x01000000),
    (0x08000000, 0x00000000),
    (0x00080000, 0x00080000),
    (0x00080800, 0x00480800),
    (0x00480008, 0x02084008),
    (0x06080808, 0x164A0848),
    (0xF2400040, 0x40104200),
    (0x00820200, 0x00001202),
    (0x00009000, 0x00000010),
    (0x00000080, 0x00000000),
    (0x80000000, 0x80000000),
    (0x80800000, 0x80800004),
    (0x80008004, 0x84008020),
    (0x808080A0, 0xA08481A4)
]

word_size = 32
max_weight_loss = 3

sols = []
btor, weight = common.speck_quasidifferential_trails(diffs, word_size)
for w in range(max_weight_loss):
    for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size):
        sols.append((w, solution))

print("Characteristic:")
for (a, b) in diffs:
    print(" " * 4, "{:08x} {:08x}".format(a, b))

for (weight, solution) in sols:
    s = common.compute_sign_speck(diffs, solution, word_size)
    print("Weight: {} [{:+}]".format(weight, s))
    for (u, v) in solution:
        print(" " * 4, "{:08x} {:08x}".format(u, v))
