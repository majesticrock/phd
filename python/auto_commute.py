from methods_for_commute import *

N_U  = np.array([Operator(Momentum("k", True, False), True,  True),   Operator(Momentum("k", True, False), True,  False)])
N_D  = np.array([Operator(Momentum("k", False, False), False, True),   Operator(Momentum("k", False, False), False, False)])
SC   = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True,  False)])
SC_D = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True,  False)])
dagger_it(SC_D)
CDW_U   = np.array([Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)])
CDW_U_D = np.array([Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)])
CDW_D   = np.array([Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, True), False, False)])
CDW_D_D = np.array([Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, True), False, False)])
dagger_it(CDW_U_D)
dagger_it(CDW_D_D)
ETA   = np.array([Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)])
ETA_D = np.array([Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)])
dagger_it(ETA_D)

commuted_with_H = Expression(1, np.array([], dtype=Term))
ex_l = Expression(1, np.array([ Term(1, "", N_U), Term(1, "", N_D) ]))
ex_r = Expression(1, np.array([ Term(1, "", N_U), Term(1, "", N_D) ]))

for i in range(0, ex_r.terms.size):
    left  = ex_r.terms[i].operators[0]
    right = ex_r.terms[i].operators[1]
    commute_bilinear_with_H(left, right, commuted_with_H)

commuted_with_H.normalOrder()
print(f"\\left[H, {ex_r}\\right] &= {commuted_with_H}\n")

c = anti_commmute(ex_l, commuted_with_H)
c.normalOrder()
print(f"\\left\\{{ {ex_l}, \\left[ H, {ex_r} \\right] \\right\\}} &= {c.as_expectation_values()}")