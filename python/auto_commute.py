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
ETA_Q   = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, True), True, False)])
ETA_Q_D = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, True), True, False)])
dagger_it(ETA_D)
dagger_it(ETA_Q_D)

H = Expression(1, np.array([], dtype=Term))
for k in [Momentum("k", True, False), Momentum("k", True, True), Momentum("k", False, False), Momentum("k", False, True)]:
    H.append(Term(*sync_eps(k), np.array([Operator(k.copy(), True, True), Operator(k.copy(), True, False)])))
    H.append(Term(*sync_eps(k), np.array([Operator(k.copy(), False, True), Operator(k.copy(), False, False)])))

    l = k.copy()
    l.positive = not l.positive
    H.append(Term(1, r"\Delta_\text{SC}", np.array([Operator(l.copy(), False, False), Operator(k.copy(), True, False)])))
    H.append(Term(1, r"\Delta_\text{SC}", np.array([Operator(k.copy(), True, True), Operator(l.copy(), False, True)])))

    l.additional_Q()
    H.append(Term(1, r"\Delta_\eta^*"  , np.array([Operator(l.copy(), False, False), Operator(k.copy(), True, False)])))
    H.append(Term(1, r"\Delta_\eta", np.array([Operator(k.copy(), True, True), Operator(l.copy(), False, True)])))

    l.positive = not l.positive
    H.append(Term(1, r"\Delta_\text{CDW}", np.array([Operator(k.copy(), True, True), Operator(l.copy(), True, False)])))
    H.append(Term(1, r"\Delta_\text{CDW}", np.array([Operator(k.copy(), False, True), Operator(l.copy(), False, False)])))

H.sortByCoefficient()

ex_l = Expression(1, np.array([ Term(1, "", SC) ]))
ex_r = Expression(1, np.array([ Term(1, "", N_D), Term(1, "", N_U) ]))
commuted_with_H = commute(H, ex_r)
commuted_with_H.normalOrder()
commuted_with_H.sortByCoefficient()
print(f"\\begin{{align*}}\n\\left[H, {ex_r}\\right] &= {commuted_with_H}\n\\end{{align*}}\n")

c = anti_commute(ex_l, commuted_with_H)
c.normalOrder()
print(f"\\begin{{align*}}\n\\left\\{{ {ex_l}, \\left[ H, {ex_r} \\right] \\right\\}} &= {c}\n\\end{{align*}}\n")