import common

# 15 rounds characteristic
diffs = [
    (0x082020000000, 0x000120200000),
    (0x000900000000, 0x000001000000),
    (0x000008000000, 0x000000000000),
    (0x000000080000, 0x000000080000),
    (0x000000080800, 0x000000480800),
    (0x000000480008, 0x000002084008),
    (0x0800FE080808, 0x0800EE4A0848),
    (0x000772400040, 0x400000104200),
    (0x000000820200, 0x000000001202),
    (0x000000009000, 0x000000000010),
    (0x000000000080, 0x000000000000),
    (0x800000000000, 0x800000000000),
    (0x808000000000, 0x808000000004),
    (0x800080000004, 0x840080000020),
    (0x808080800020, 0xA08480800124),
    (0x800400008124, 0x842004008801)
]

word_size = 48
max_weight_loss = 1

sols = []
btor, weight = common.speck_quasidifferential_trails(diffs, word_size)
for w in range(max_weight_loss):
    for solution in common.solve_all(btor, weight, w, len(diffs) - 1, word_size):
        sols.append((w, solution))

print("Characteristic:")
for (a, b) in diffs:
    print(" " * 4, "{:012x} {:012x}".format(a, b))

for (weight, solution) in sols:
    s = common.compute_sign_speck(diffs, solution, word_size)
    print("Weight: {} [{:+}]".format(weight, s))
    for (u, v) in solution:
        print(" " * 4, "{:012x} {:012x}".format(u, v))
