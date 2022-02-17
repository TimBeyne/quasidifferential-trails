import itertools

def interleave_bits(x, y, n):
    z = 0
    for i in range(n):
        z |= (x & (1 << i)) << i | (y & (1 << i)) << (i + 1);
    return z

def to_quasidifferential_basis(x):
    if len(x) == 1:
        return x

    assert len(x) % 4 == 0

    l = int(len(x) / 4) 

    x_00 = to_quasidifferential_basis(x[   :  l])
    x_01 = to_quasidifferential_basis(x[  l:2*l])
    x_10 = to_quasidifferential_basis(x[2*l:3*l])
    x_11 = to_quasidifferential_basis(x[3*l:   ])

    return vector(x.base_ring(),
        (x_00 + x_11).list() +\
        (x_01 + x_10).list() +\
        (x_00 - x_11).list() +\
        (x_01 - x_10).list())

def interleaved_transition_matrix(F, n, m):
    T = matrix(QQ, 2 ** (2*m), 2 ** (2*n))
    for x in range(2 ** n):
        for y in range(2 ** n):
            i = interleave_bits(x, y, n)
            j = interleave_bits(F(x), F(y), m)
            T[j, i] = 1
    return T

def quasidifferential_transition_matrix(F, n, m):
    D = interleaved_transition_matrix(F, n, m)

    # Transform columns
    for i in range(2 ** (2*n)):
        D.set_column(i, to_quasidifferential_basis(D.column(i)))
    # Transform rows
    for i in range(2 ** (2*m)):
        D.set_row(i, to_quasidifferential_basis(D.row(i)))

    return D / 2**n

def deinterleave_quasidifferential_transition_matrix(D, n, m, primary = 'diff'):
    R = matrix(QQ, 2 ** (2*m), 2 ** (2*n))
    for u, v in itertools.product(range(2 ** n), range(2 ** m)):
        for a, b in itertools.product(range(2 ** n), range(2 ** m)):
            if primary == 'mask':
                R[2 ** m * v + b, 2 ** n * u + a] = D[
                    interleave_bits(b, v, m), interleave_bits(a, u, n)
                ] 
            elif primary == 'diff':
                R[2 ** m * b + v, 2 ** n * a + u] = D[
                    interleave_bits(b, v, m), interleave_bits(a, u, n)
                ] 
    return R

# Example:
from sage.crypto.sboxes import Rectangle
D = quasidifferential_transition_matrix(Rectangle, 4, 4)
D = deinterleave_quasidifferential_transition_matrix(D, 4, 4)
print(D.str())
