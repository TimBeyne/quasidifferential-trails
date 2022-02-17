import pickle

F = GF(2)

with open("trails-0.p", "rb") as f:
    trails = pickle.load(f)

trails_by_weight = []

all_masks = []
c = []
for (s, trail) in trails:
    masks = [vector(F, u) for u in trail]
    total_weight = sum(u.list().count(1) for u in masks)
    trails_by_weight.append((total_weight, masks))

    all_masks.append(vector(F, sum((list(map(int, u)) for u in trail), [])))
    c.append(0 if s > 0 else 1)

trails_by_weight = sorted(trails_by_weight)

for weight, trail in trails_by_weight:
    print(weight)
    for u in trail:
        print("".join([str(x) for x in u]))
    print()

B = matrix(all_masks)
c = vector(F, c)
x = B.solve_right(c)

V = span(all_masks)
M = V.echelonized_basis_matrix()

for row in M.rows():
    print([(i // 64, 63 - i % 64) for i in row.nonzero_positions()])
print()
print(M * x)
