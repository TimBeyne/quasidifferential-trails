import pickle, pyboolector
import math

# 10 rounds
# First dominant characteristic
diffs = [
    0x0000000000000000000000000000000010000000000000000100000000000000,
    0x0000000000000000000000000000000040000000000000000c00000000000000,
    0x000000000000000000000000c000000000000000040000000000000000000000,
    0x00000000000000000000000040000000000000000c0000000000000000000000,
    0x0000000000000000c00000000000000004000000000000000000000000000000,
    0x000000000000000040000000000000000c000000000000000000000000000000,
    0x00000000c0000000000000000400000000000000000000000000000000000000,
    0x0000000040000000000000000c00000000000000000000000000000000000000,
    0xc000000000000000040000000000000000000000000000000000000000000000,
    0x40000000000000000c0000000000000000000000000000000000000000000000,
    0x00000000040000000000000000000000000000000000000000000000c0000000,
    0x000000000a000000000000000000000000000000000000000000000040000000,
    0x000000002000000000000000000000000000000000000000c000000000000000,
    0x0000000010000000000000000000000000000000000000002000000000000000,
    0x0000000010000000000000000000000000000000000000020000000000000000,
    0x0000000040000000000000000000000000000000000000010000000000000000,
    0x4000000000000000000000000000000000000000000000010000000000000000,
    0xa000000000000000000000000000000000000000000000040000000000000000,
    0x000000000000000000000000000000000000000c000000000000000000000002,
    0x0000000000000000000000000000000000000001000000000000000000000001,
    0x0000000000000000000000000000000000000001000000000000000000000001
]

## Second dominant characteristic
#diffs = [
#    0x0000000000000000000000000000000010000000000000000100000000000000,
#    0x0000000000000000000000000000000040000000000000000c00000000000000,
#    0x000000000000000000000000c000000000000000040000000000000000000000,
#    0x00000000000000000000000040000000000000000c0000000000000000000000,
#    0x0000000000000000c00000000000000004000000000000000000000000000000,
#    0x000000000000000040000000000000000a000000000000000000000000000000,
#    0x00000000c0000000000000000000000020000000000000000000000000000000,
#    0x0000000020000000000000000000000010000000000000000000000000000000,
#    0x0000000200000000000000000000000010000000000000000000000000000000,
#    0x0000000100000000000000000000000040000000000000000000000000000000,
#    0x0000000100000000000000004000000000000000000000000000000000000000,
#    0x000000040000000000000000c000000000000000000000000000000000000000,
#    0x000000000000000040000000000000000000000000000000000000000000000c,
#    0x0000000000000000c00000000000000000000000000000000000000000000004,
#    0x0000000040000000000000000000000000000000000000000000000c00000000,
#    0x00000000c0000000000000000000000000000000000000000000000400000000,
#    0x40000000000000000000000000000000000000000000000c0000000000000000,
#    0xa000000000000000000000000000000000000000000000040000000000000000,
#    0x000000000000000000000000000000000000000c000000000000000000000002,
#    0x0000000000000000000000000000000000000001000000000000000000000001,
#    0x0000000000000000000000000000000000000001000000000000000000000001
#]
## Third dominant characteristic
#diffs = [
#    0x0000000000000000000000000000000010000000000000000100000000000000,
#    0x0000000000000000000000000000000040000000000000000c00000000000000,
#    0x000000000000000000000000c000000000000000040000000000000000000000,
#    0x00000000000000000000000040000000000000000a0000000000000000000000,
#    0x0000000000000000c00000000000000000000000200000000000000000000000,
#    0x0000000000000000200000000000000000000000100000000000000000000000,
#    0x0000000000000002000000000000000000000000100000000000000000000000,
#    0x0000000000000001000000000000000000000000400000000000000000000000,
#    0x0000000000000001000000000000000040000000000000000000000000000000,
#    0x00000000000000040000000000000000c0000000000000000000000000000000,
#    0x0000000c00000000000000004000000000000000000000000000000000000000,
#    0x000000040000000000000000c000000000000000000000000000000000000000,
#    0x000000000000000040000000000000000000000000000000000000000000000c,
#    0x0000000000000000c00000000000000000000000000000000000000000000004,
#    0x0000000040000000000000000000000000000000000000000000000c00000000,
#    0x00000000c0000000000000000000000000000000000000000000000400000000,
#    0x40000000000000000000000000000000000000000000000c0000000000000000,
#    0xa000000000000000000000000000000000000000000000040000000000000000,
#    0x000000000000000000000000000000000000000c000000000000000000000002,
#    0x0000000000000000000000000000000000000001000000000000000000000001,
#    0x0000000000000000000000000000000000000001000000000000000000000001
#]
## Fourth dominant characteristic
#diffs = [
#    0x0000000000000000000000000000000010000000000000000100000000000000,
#    0x0000000000000000000000000000000040000000000000000c00000000000000,
#    0x000000000000000000000000c000000000000000040000000000000000000000,
#    0x00000000000000000000000040000000000000000c0000000000000000000000,
#    0x0000000000000000c00000000000000004000000000000000000000000000000,
#    0x000000000000000040000000000000000c000000000000000000000000000000,
#    0x00000000c0000000000000000400000000000000000000000000000000000000,
#    0x0000000040000000000000000c00000000000000000000000000000000000000,
#    0xc000000000000000040000000000000000000000000000000000000000000000,
#    0x40000000000000000a0000000000000000000000000000000000000000000000,
#    0x00000000000000002000000000000000000000000000000000000000c0000000,
#    0x0000000000000000100000000000000000000000000000000000000020000000,
#    0x0000000000000000100000000000000000000000000000000000000200000000,
#    0x0000000000000000400000000000000000000000000000000000000100000000,
#    0x0000000040000000000000000000000000000000000000000000000100000000,
#    0x00000000c0000000000000000000000000000000000000000000000400000000,
#    0x40000000000000000000000000000000000000000000000c0000000000000000,
#    0xa000000000000000000000000000000000000000000000040000000000000000,
#    0x000000000000000000000000000000000000000c000000000000000000000002,
#    0x0000000000000000000000000000000000000001000000000000000000000001,
#    0x0000000000000000000000000000000000000001000000000000000000000001
#]
## Fifth dominant characteristic
#diffs = [
#    0x0000000000000000000000000000000010000000000000000100000000000000,
#    0x0000000000000000000000000000000040000000000000000c00000000000000,
#    0x000000000000000000000000c000000000000000040000000000000000000000,
#    0x00000000000000000000000040000000000000000c0000000000000000000000,
#    0x0000000000000000c00000000000000004000000000000000000000000000000,
#    0x000000000000000040000000000000000c000000000000000000000000000000,
#    0x00000000c0000000000000000400000000000000000000000000000000000000,
#    0x0000000040000000000000000a00000000000000000000000000000000000000,
#    0xc000000000000000000000002000000000000000000000000000000000000000,
#    0x2000000000000000000000001000000000000000000000000000000000000000,
#    0x0000000000000000000000001000000000000000000000000000000000000002,
#    0x0000000000000000000000004000000000000000000000000000000000000001,
#    0x0000000000000000400000000000000000000000000000000000000000000001,
#    0x0000000000000000c00000000000000000000000000000000000000000000004,
#    0x0000000040000000000000000000000000000000000000000000000c00000000,
#    0x00000000c0000000000000000000000000000000000000000000000400000000,
#    0x40000000000000000000000000000000000000000000000c0000000000000000,
#    0xa000000000000000000000000000000000000000000000040000000000000000,
#    0x000000000000000000000000000000000000000c000000000000000000000002,
#    0x0000000000000000000000000000000000000001000000000000000000000001,
#    0x0000000000000000000000000000000000000001000000000000000000000001
#]

nb_bits_state = 256
nb_bits_cost = 16
min_weight = 0
max_weight = 1

constants = [
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x41, 0x03, 0x06,
    0x0c, 0x18, 0x30, 0x61, 0x42, 0x05, 0x0a, 0x14, 0x28, 0x51, 0x23, 0x47,
    0x0f, 0x1e, 0x3c, 0x79, 0x72, 0x64, 0x48, 0x11, 0x22, 0x45, 0x0b, 0x16,
    0x2c, 0x59, 0x33, 0x67, 0x4e, 0x1d, 0x3a, 0x75, 0x6a, 0x54, 0x29, 0x53,
    0x27, 0x4f, 0x1f, 0x3e, 0x7d, 0x7a, 0x74, 0x68, 0x50, 0x21, 0x43, 0x07,
    0x0e, 0x1c, 0x38, 0x71, 0x62, 0x44, 0x09, 0x12, 0x24, 0x49, 0x13, 0x26,
    0x4d, 0x1b, 0x36, 0x6d, 0x5a, 0x35, 0x6b, 0x56, 0x2d, 0x5b, 0x37, 0x6f,
    0x5e, 0x3d, 0x7b, 0x76, 0x6c, 0x58, 0x31, 0x63, 0x46, 0x0d, 0x1a, 0x34,
    0x69, 0x52, 0x25, 0x4b, 0x17, 0x2e, 0x5d, 0x3b, 0x77, 0x6e, 0x5c, 0x39,
    0x73, 0x66, 0x4c, 0x19, 0x32, 0x65, 0x4a, 0x15, 0x2a, 0x55, 0x2b, 0x57,
    0x2f, 0x5f, 0x3f, 0x7f, 0x7e, 0x7c, 0x78, 0x70, 0x60, 0x40
]


nb_rounds = (len(diffs) - 1) // 2

weight_tables = pickle.load(open("weights_knot.p", "rb"))
D = pickle.load(open("qdt_knot.p", "rb"))

btor = pyboolector.Boolector()
btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)

us = [btor.Var(btor.BitVecSort(nb_bits_state), "u%d" % i) for i in range(nb_rounds + 1)]
vs = [btor.Var(btor.BitVecSort(nb_bits_state), "v%d" % i) for i in range(nb_rounds)]

def compute_sign(diffs, trail):
    correlation = 1
    for i in range(nb_rounds):
        a = [(diffs[2*i] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        b = [(diffs[2*i + 1] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        u = [(trail[2*i] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        v = [(trail[2*i + 1] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        for j in range(64):
            c = D[16*b[j]+v[j], 16*a[j]+u[j]]
            if c == 0:
                print(a[j], b[j], u[j], v[j])
                correlation = 0
            elif c < 0:
                correlation *= -1
    return correlation

def parity(constant, u):
    return bin(constant & u).count('1') % 2

def compute_correlation(diffs, trail):
    correlation = 1
    for i in range(nb_rounds):
        a = [(diffs[2*i] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        b = [(diffs[2*i + 1] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        u = [(trail[2*i] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        v = [(trail[2*i + 1] >> j) & 0xf for j in range(0, nb_bits_state, 4)]
        for j in range(64):
            assert D[16*b[j], 16*a[j]] != 0
            correlation *= D[16*b[j]+v[j], 16*a[j]+u[j]] / D[16*b[j], 16*a[j]] # Relative
        # Sign due to constant addition
        correlation *= (-1) ** parity(
            constants[i],
            ((u[7] & 1) << 7) | ((u[6] & 1) << 6) | ((u[5] & 1) << 5) | ((u[4] & 1) << 4) |
            ((u[3] & 1) << 3) | ((u[2] & 1) << 2) | ((u[1] & 1) << 1) | (u[0] & 1)
        )
    return correlation

def permute_bits(x, y):
    for i in range(64):
        btor.Assert(y[4*i    ] == x[4*i                ])
        btor.Assert(y[4*i + 1] == x[4*((i-1)  % 64) + 1])
        btor.Assert(y[4*i + 2] == x[4*((i-8)  % 64) + 2])
        btor.Assert(y[4*i + 3] == x[4*((i-25) % 64) + 3])

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
correlation = 0
for target in range(min_weight, max_weight):
    mask_strings = []
    # Find all solutions
    previous = []
    while True:
        btor.Assume(cost == target)
        distinct = btor.Const(1)
        for (_, ws) in previous:
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
                print("u{:2} = {:064x}".format(i, u))
                print("v{:2} = {:064x}".format(i, v))
                trail.append(u)
                trail.append(v)
            s = compute_sign(diffs, trail)
            trail.append(int(us[nb_rounds].assignment, base = 2))
            print("u{:2} = {:064x}".format(nb_rounds, int(us[nb_rounds].assignment, base = 2)))
            print("Sign: {}".format(s))
            correlation += compute_correlation(diffs, trail)
            if correlation != 0:
                print("Correlation: {} ({})".format(math.log2(abs(correlation)), correlation/abs(correlation)))
            else:
                print("Correlation: zero.")
            previous.append((s, [us[i].assignment for i in range(1, nb_rounds)]))
        else:
            with open("trails-{}.p".format(target), "wb") as f:
                pickle.dump(previous, f)
            print("No trails with weight equal to {}.".format(target))
            break

