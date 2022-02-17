import pickle, pyboolector

# 14 rounds (key-recovery, weight 63)
diffs = [
    0x0020000600000000,
    0x0060000200000000,
    0x0200006000000000,
    0x0600002000000000,
    0x2000060000000000,
    0x6000020000000000,
    0x0000600000000002,
    0x0000200000000006,
    0x0006000000000020,
    0x0002000000000060,
    0x0060000000000200,
    0x0020000000000600,
    0x0600000000002000,
    0x0200000000006000,
    0x6000000000020000,
    0x2000000000060000,
    0x0000000000200006,
    0x0000000000600002,
    0x0000000002000060,
    0x000000000c000020,
    0x0000000000008600,
    0x0000000000001200,
    0x0000000000003000,
    0x0000000000008000,
    0x0000000000000008,
    0x0000000000000001,
    0x0000000000000001,
    0x0000000000000006,
    0x0004000000000020
]

# 14 rounds (key-recovery, second characteristic, weight 66)
#diffs = [
#    0x0020000600000000,
#    0x0060000200000000,
#    0x0200006000000000,
#    0x0600002000000000,
#    0x2000060000000000,
#    0x6000020000000000,
#    0x0000600000000002,
#    0x0000200000000006,
#    0x0006000000000020,
#    0x0002000000000060,
#    0x0060000000000200,
#    0x0020000000000600,
#    0x0600000000002000,
#    0x0200000000006000,
#    0x6000000000020000,
#    0x2000000000060000,
#    0x0000000000200006,
#    0x0000000000600002,
#    0x0000000002000060,
#    0x000000000c000020,
#    0x0000000000008600,
#    0x0000000000009200,
#    0x0000000000003008,
#    0x0000000000008001,
#    0x0000000000000009,
#    0x0000000000000001,
#    0x0000000000000001,
#    0x0000000000000006,
#    0x0004000000000020
#]

nb_bits_state = 64
nb_bits_cost = 16
min_weight = 0
max_weight = 5

nb_rounds = (len(diffs) - 1) // 2

weight_tables = pickle.load(open("weights_rectangle.p", "rb"))
D = pickle.load(open("qdt_rectangle.p", "rb"))

btor = pyboolector.Boolector()
btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)

us = [btor.Var(btor.BitVecSort(nb_bits_state), "u%d" % i) for i in range(nb_rounds + 1)]
vs = [btor.Var(btor.BitVecSort(nb_bits_state), "v%d" % i) for i in range(nb_rounds)]

def compute_sign(diffs, trail):
    correlation = 1
    for i in range(nb_rounds):
        a = [(diffs[2*i] >> j) & 0xf for j in range(0, 64, 4)]
        b = [(diffs[2*i + 1] >> j) & 0xf for j in range(0, 64, 4)]
        u = [(trail[2*i] >> j) & 0xf for j in range(0, 64, 4)]
        v = [(trail[2*i + 1] >> j) & 0xf for j in range(0, 64, 4)]
        for j in range(16):
            c = D[16*b[j]+v[j], 16*a[j]+u[j]]
            if c == 0:
                print(a[j], b[j], u[j], v[j])
                correlation = 0
            elif c < 0:
                correlation *= -1
    return correlation

def permute_bits(x, y):
    for i in range(16):
        btor.Assert(y[4*i    ] == x[4*i                ])
        btor.Assert(y[4*i + 1] == x[4*((i-1)  % 16) + 1])
        btor.Assert(y[4*i + 2] == x[4*((i-12) % 16) + 2])
        btor.Assert(y[4*i + 3] == x[4*((i-13) % 16) + 3])

def qdt_constraints(u, v, i, j):
    a = (diffs[2*i  ] >> j) & 0xf
    b = (diffs[2*i+1] >> j) & 0xf

    # Distinguish between two cases:
    # - If a = b = 0, then we essentially model the correlation matrix
    # - Otherwise, the weight is always zero when the transition is allowed
    if a == 0 and b == 0:
        weight0 = (u == 0) & (v == 0)
        weight1 = btor.Const(0)
        for (x, y) in weight_tables[(0, 0)][1]:
            weight1 |= (u == y) & (v == x)
        weight2 = btor.Const(0)
        for (x, y) in weight_tables[(0, 0)][2]:
            weight2 |= (u == y) & (v == x)
        btor.Assert(weight0 | weight1 | weight2)
        return btor.Cond(
            weight1,
            btor.Const(1, nb_bits_cost),
            btor.Cond(weight2, btor.Const(2, nb_bits_cost), btor.Const(0, nb_bits_cost))
        )
    else:
        allowed = btor.Const(0)
        for (x, y) in weight_tables[(b, a)][0]:
            allowed |= (u == y) & (v == x)
        btor.Assert(allowed)
        return btor.Const(0, nb_bits_cost)

cost = btor.Const(0, nb_bits_cost)
for i in range(nb_rounds):
    permute_bits(vs[i], us[i + 1])
    for j in range(0, nb_bits_state, 4):
        sbox_in  = us[i][j+3:j]
        sbox_out = vs[i][j+3:j]
        w = qdt_constraints(sbox_in, sbox_out, i, j)
        cost += w

btor.Assert(us[0] == 0)
btor.Assert(us[nb_rounds] == 0)
btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)
btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

for target in range(min_weight, max_weight):
    mask_strings = []
    # Find all solutions
    previous = []
    while True:
        btor.Assume(cost == target)
        distinct = btor.Const(1)
        for _, ws in previous:
            temp = btor.Const(0)
            for i in range(1, nb_rounds):
                temp |= (us[i] != btor.Const(ws[i - 1], nb_bits_state))
            distinct &= temp
        btor.Assume(distinct)
        
        r = btor.Sat()
        if r == btor.SAT:
            print("Solution: [#{} of weight {}]".format(len(previous) + 1, target))
            trail = []
            for i in range(nb_rounds):
                u = int(us[i].assignment, base = 2)
                v = int(vs[i].assignment, base = 2)
                print("u{:2} = {:016x}".format(i, u))
                print("v{:2} = {:016x}".format(i, v))
                trail.append(u)
                trail.append(v)
            s = compute_sign(diffs, trail)
            trail.append(int(us[nb_rounds].assignment, base = 2))
            print("u{:2} = {:016x}".format(nb_rounds, int(us[nb_rounds].assignment, base = 2)))
            print("Sign: {}".format(s))
            previous.append((s, [us[i].assignment for i in range(1, nb_rounds)]))
        else:
            with open("trails-{}.p".format(target), "wb") as f:
                pickle.dump(previous, f)
            print("No trails with weight equal to {}.".format(target))
            break

