import pyboolector
from pyboolector import Boolector

def is_power_of_two(x):
    return x & (x - 1) == 0

def next_power_of_two(x):
    """ Suppose x < 2^{32} """
    x -= 1
    x |= x >> 1
    x |= x >> 2
    x |= x >> 4
    x |= x >> 8
    x |= x >> 16
    x += 1
    return x

def shift_left(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return x << offset
    else:
        n = next_power_of_two(nb_bits)
        y = btor.Concat(btor.Const(0, n - nb_bits), x) << offset
        return y[nb_bits - 1:]

def shift_right(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return x >> offset
    else:
        n = next_power_of_two(nb_bits)
        y = btor.Concat(btor.Const(0, n - nb_bits), x) >> offset
        return y[nb_bits - 1:]

def rotate_right(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return btor.Ror(x, offset)
    else:
        return shift_left(btor, x, nb_bits - offset, nb_bits) | \
               shift_right(btor, x, offset, nb_bits)

def rotate_left(btor, x, offset, nb_bits):
    if is_power_of_two(nb_bits):
        return btor.Rol(x, offset)
    else:
        return shift_left(btor, x, offset, nb_bits) | \
               shift_right(btor, x, nb_bits - offset, nb_bits)

def M_pseudoinverse(btor, t, nb_bits):
    t = t ^ shift_left(btor, t, 1, nb_bits)
    return shift_right(btor, t, 1, nb_bits)

def M_transpose(btor, t, nb_bits):
    t = shift_right(btor, t, 1, nb_bits)
    i = 1
    while i < nb_bits:
        t = t ^ shift_right(btor, t, i, nb_bits)
        i *= 2
    return t

def hamming_weight_16(btor, x, nb_bits):
    x -= shift_right(btor, x, 1, nb_bits) & 0x5555
    x = (x & 0x3333) + (shift_right(btor, x, 2, nb_bits) & 0x3333)
    x = (x + shift_right(btor, x, 4, nb_bits)) & 0x0F0F
    x += shift_right(btor, x, 8, nb_bits)
    return x & 0x003F

def hamming_weight_32(btor, x, nb_bits):
    x -= shift_right(btor, x, 1, nb_bits) & 0x55555555
    x = (x & 0x33333333) + (shift_right(btor, x, 2, nb_bits) & 0x33333333)
    x = (x + shift_right(btor, x, 4, nb_bits)) & 0x0F0F0F0F
    x += shift_right(btor, x, 8, nb_bits)
    x += shift_right(btor, x, 16, nb_bits)
    return x & 0x0000003F

def hamming_weight(btor, x, nb_bits):
    if nb_bits == 16:
        return hamming_weight_16(btor, x, nb_bits)
    elif nb_bits == 24:
        return hamming_weight(btor, btor.Concat(btor.Const(0, 8), x), 32)[23:]
    elif nb_bits == 32:
        return hamming_weight_32(btor, x, nb_bits)
    elif nb_bits == 48:
        return btor.Concat(btor.Const(0, 16), hamming_weight(btor, x[31:], 32))\
             + btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[:32], 16))
    elif nb_bits == 64:
        return btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[31:], 32))\
             + btor.Concat(btor.Const(0, 32), hamming_weight(btor, x[:32], 32))

        
def bitwise_and(btor, a, b, c, u, v, w, nb_bits):
    a = btor.Const(a, nb_bits)
    b = btor.Const(b, nb_bits)
    c = btor.Const(c, nb_bits)

    btor.Assert((u | v) & ~(a | b | w) == 0)
    btor.Assert((a & u) ^ (b & v) == (c & w))

    return hamming_weight(btor, w & ~a & ~b, nb_bits)

def modular_addition(btor, a, b, c, u, v, w, nb_bits, i):
    a_ = b ^ c
    b_ = a ^ c
    c_ = M_pseudoinverse(btor, a ^ b ^ c, nb_bits)

    u_ = u ^ w
    v_ = v ^ w
    w_ = M_transpose(btor, u ^ v ^ w, nb_bits)

    n = nb_bits - 1

    btor.Assert((u_ | v_) & ~(a_ | b_ | w_) == 0)
    btor.Assert((a_ & u_) ^ (b_ & v_) == (c_ & w_))
    btor.Assert(((a_[n] == 0) & (b_[n] == 0)) | (a_[n] & u_[n] == u_[n] ^ v_[n]))

    weight = hamming_weight(btor, w_ & ~a_ & ~b_, nb_bits) 
    extra = btor.Cond(
        (a_[n] | b_[n]) & ((u_[n] ^ v_[n]) == (a_[n] & u_[n])) & (u_[n] != v_[n]),
        btor.Const(1, nb_bits),
        btor.Const(0, nb_bits)
    )
    return weight - extra

def speck_quasidifferential_trails(diffs, nb_bits):
    btor = Boolector()
    btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)
    btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

    nb_rounds = len(diffs) - 1

    u = [btor.Var(btor.BitVecSort(nb_bits), "u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Var(btor.BitVecSort(nb_bits), "v%d" % i) for i in range(nb_rounds + 1)]
    
    weight = btor.Const(0, nb_bits)
    for i in range(nb_rounds):
        a = btor.Const(diffs[i][0], nb_bits)
        b = btor.Const(diffs[i][1], nb_bits)
        c = btor.Const(diffs[i + 1][0], nb_bits)
        if nb_bits == 16:
            u_ = btor.Ror(u[i], 7)
            v_ = btor.Ror(v[i + 1], 2) ^ v[i] 
            a  = btor.Ror(a, 7)
        else:
            u_ = rotate_right(btor, u[i], 8, nb_bits)
            v_ = rotate_right(btor, v[i + 1], 3, nb_bits) ^ v[i] 
            a  = rotate_right(btor, a, 8, nb_bits)
        w_ = u[i + 1] ^ v[i + 1]

        weight += modular_addition(
            btor, a, b, c, u_, v_, w_, nb_bits, i
        )
    
    btor.Assert(u[0] == 0)
    btor.Assert(v[0] == 0)
    btor.Assert(u[nb_rounds] == 0)
    btor.Assert(v[nb_rounds] == 0)

    return btor, weight

def simon_quasidifferential_trails(diffs, nb_bits):
    btor = Boolector()
    btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)
    btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

    nb_rounds = len(diffs) - 1

    u = [btor.Var(btor.BitVecSort(nb_bits), "u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Var(btor.BitVecSort(nb_bits), "v%d" % i) for i in range(nb_rounds + 1)]
    w = [btor.Var(btor.BitVecSort(nb_bits), "w%d" % i) for i in range(nb_rounds)]
    
    weight = btor.Const(0, nb_bits)
    for i in range(nb_rounds):
        u_ = btor.Rol(u[i] ^ w[i], 1)
        v_ = btor.Rol(btor.Ror(v[i], 2) ^ w[i] ^ v[i + 1], 8)
        t = btor.Const(diffs[i][0], nb_bits)
        a = btor.Rol(t, 1)
        b = btor.Rol(t, 8)
        c = btor.Rol(t, 2) ^ diffs[i + 1][0] ^ diffs[i][1]
        weight += bitwise_and(
            btor, a, b, c, u_, v_, v[i], nb_bits
        )
        btor.Assert(u[i + 1] == v[i])
    
    btor.Assert(u[0] == 0)
    btor.Assert(v[0] == 0)
    btor.Assert(u[nb_rounds] == 0)
    btor.Assert(v[nb_rounds] == 0)
    
    return btor, weight

def distinctness_constraint(btor, masks, solutions, word_size):
    distinctness_condition = btor.Const(1)
    for solution in solutions:
        condition = btor.Const(0)
        for i in range(len(masks)):
            for j in range(len(masks[i])):
                condition |= (masks[i][j] != btor.Const(solution[j][i], word_size))
        distinctness_condition &= condition
    btor.Assume(distinctness_condition)

def solve_all(btor, weight, w, nb_rounds, word_size):
    # Get variables
    u = [btor.Match_by_symbol("u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Match_by_symbol("v%d" % i) for i in range(nb_rounds + 1)]

    solutions = []
    while True:
        btor.Assume(weight == w)
        distinctness_constraint(btor, [u, v], solutions, word_size)

        r = btor.Sat()
        if r != btor.SAT:
            return solutions

        solutions.append([
            (int(u[i].assignment, base = 2), int(v[i].assignment, base = 2))\
            for i in range(nb_rounds + 1)
        ])

def solve_all_simon(btor, weight, target_weight, nb_rounds, word_size):
    # Get variables
    u = [btor.Match_by_symbol("u%d" % i) for i in range(nb_rounds + 1)]
    v = [btor.Match_by_symbol("v%d" % i) for i in range(nb_rounds + 1)]
    w = [btor.Match_by_symbol("w%d" % i) for i in range(nb_rounds)]

    solutions = []
    while True:
        btor.Assume(weight == target_weight)
        distinctness_constraint(btor, [u, v, w], solutions, word_size)

        r = btor.Sat()
        if r != btor.SAT:
            return solutions

        solutions.append([
            (int(u[i].assignment, base = 2),
             int(v[i].assignment, base = 2),
             int(w[i].assignment, base = 2) if i < nb_rounds else None)\
            for i in range(nb_rounds + 1)
        ])

def parity(x):
    return bin(x).count('1') % 2

def rotl(x, r, word_size):
   mask = (1 << word_size) - 1
   return ((x << r) & mask) | ((x & mask) >> (word_size - r))
    

def compute_sign_speck(differences, masks, word_size):
    """ Compute the sign of a quasidifferential trail for Speck. """

    def pseudoinverseM(t):
        t = t ^ ((t << 1) % 2 ** word_size)
        return t >> 1

    def complement(t):
        return (2 ** word_size - 1) ^ t

    s = 1
    for i in range(len(masks) - 1):
        if word_size == 16:
            u = rotl(masks[i][0], word_size - 7, word_size)
            v = rotl(masks[i + 1][1], word_size - 2, word_size) ^ masks[i][1]
        else:
            u = rotl(masks[i][0], word_size - 8, word_size)
            v = rotl(masks[i + 1][1], word_size - 3, word_size) ^ masks[i][1]
        w = masks[i + 1][0] ^ masks[i + 1][1]
        (u_, v_, _) = (u ^ w, v ^ w, u ^ v ^ w)

        (l , b) = differences[i]
        c       = differences[i + 1][0]
        if word_size == 16:
            a = rotl(l, word_size - 7, word_size) 
        else:
            a = rotl(l, word_size - 8, word_size) 
        (a_, b_, c_) = (b ^ c, a ^ c, pseudoinverseM(a ^ b ^ c))
        p1 = parity(((complement(a_) & u_) ^ (c_ & v_)) & ((complement(b_) & v_) ^ (c_ & u_)))
        p2 = parity((u_ & v_) & (c_ ^ (a_ & b_ & complement(c_))))
        s *= (-1) ** (p1 ^ p2)

    return s
