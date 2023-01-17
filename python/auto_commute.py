from methods_for_commute import *
from copy import deepcopy

N_U    = np.array([Operator(Momentum("k", True, False), True,  True),   Operator(Momentum("k", True, False), True,  False)])
N_D    = np.array([Operator(Momentum("k", False, False), False, True),   Operator(Momentum("k", False, False), False, False)])
N_Q_U  = np.array([Operator(Momentum("k", True,  True), True,  True),   Operator(Momentum("k", True,  True), True,  False)])
N_Q_D  = np.array([Operator(Momentum("k", False, True), False, True),   Operator(Momentum("k", False, True), False, False)])
N_K   = Expression(1, np.array([Term(0.5, "", N_U),   Term(0.5, "", N_D)]))
N_K_Q = Expression(1, np.array([Term(0.5, "", N_Q_U), Term(0.5, "", N_Q_D)]))

SC     = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True,  False)])
SC_D   = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True,  False)])
SC_Q   = np.array([Operator(Momentum("k", False, True),  False, False), Operator(Momentum("k", True, True),  True,  False)])
SC_Q_D = np.array([Operator(Momentum("k", False, True),  False, False), Operator(Momentum("k", True, True),  True,  False)])
dagger_it(SC_D)
dagger_it(SC_Q_D)
F_K     = Expression(1, np.array([Term(1, "", SC    )]))
F_K_D   = Expression(1, np.array([Term(1, "", SC_D  )]))
F_K_Q   = Expression(1, np.array([Term(1, "", SC_Q  )]))
F_K_Q_D = Expression(1, np.array([Term(1, "", SC_Q_D)]))

CDW_U   = np.array([Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)])
CDW_U_D = np.array([Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)])
CDW_D   = np.array([Operator(Momentum("k", False, False), False, True), Operator(Momentum("k", False, True), False, False)])
CDW_D_D = np.array([Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, True), False, False)])
dagger_it(CDW_U_D)
dagger_it(CDW_D_D)
G_K   = Expression(1, np.array([Term(0.5, "", CDW_U),   Term(0.5, "", CDW_D)]))
G_K_D = Expression(1, np.array([Term(0.5, "", CDW_U_D), Term(0.5, "", CDW_D_D)]))

ETA     = np.array([Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)])
ETA_D   = np.array([Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)])
ETA_Q   = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, True), True, False)])
ETA_Q_D = np.array([Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, True), True, False)])
dagger_it(ETA_D)
dagger_it(ETA_Q_D)
ETA_K   = Expression(1, np.array([Term(0.5, "", ETA),   Term(0.5, "", ETA_Q)]))
ETA_K_D = Expression(1, np.array([Term(0.5, "", ETA_D), Term(0.5, "", ETA_Q_D)]))

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

basis = np.array([F_K, F_K_Q, F_K_D, F_K_Q_D, ETA_K, ETA_K_D, G_K, G_K_D, N_K, N_K_Q])
basis = np.array([
    Expression(1, np.array([Term(1, "", SC    )])),
    Expression(1, np.array([Term(1, "", SC_D  )])),
    Expression(1, np.array([Term(1, "", SC_Q  )])),
    Expression(1, np.array([Term(1, "", SC_Q_D)])),
    Expression(1, np.array([Term(1, "", N_U  )])),
    Expression(1, np.array([Term(1, "", N_D  )])),
    Expression(1, np.array([Term(1, "", N_Q_U)])),
    Expression(1, np.array([Term(1, "", N_Q_D)])),
    Expression(1, np.array([Term(1, "", CDW_U  )])),
    Expression(1, np.array([Term(1, "", CDW_U_D)])),
    Expression(1, np.array([Term(1, "", CDW_D  )])),
    Expression(1, np.array([Term(1, "", CDW_D_D)])),
    Expression(1, np.array([Term(1, "", ETA    )])),
    Expression(1, np.array([Term(1, "", ETA_D  )])),
    Expression(1, np.array([Term(1, "", ETA_Q  )])),
    Expression(1, np.array([Term(1, "", ETA_Q_D)]))
])
#basis = np.array([
#    Expression(1, np.array([Term(1, "", np.array([Operator(Momentum("k", True, False), True, True)]))])), #c_k^+
#    Expression(1, np.array([Term(1, "", np.array([Operator(Momentum("k", True, True), True, True)]))])),
#    Expression(1, np.array([Term(1, "", np.array([Operator(Momentum("k", False, False), False, False)]))])), #c_-k
#    Expression(1, np.array([Term(1, "", np.array([Operator(Momentum("k", False, True), False, False)]))]))
#])
N = ""
M = ""

for i in range(0, basis.size):
    ex_r = deepcopy(basis[i])
    commuted_with_H = commute(H, ex_r)
    commuted_with_H.normalOrder()
    commuted_with_H.sortByCoefficient()
    for j in range(i, basis.size):
        ex_r = deepcopy(basis[i])
        ex_l = deepcopy(basis[j])
        ex_l.hermitianConjugate()
        c = commute(ex_l, ex_r)
        d = commute(ex_l, commuted_with_H)
        c.normalOrder()
        buffer = c.as_data()
        #if buffer != "":
        N += f"[ # {{ {ex_l}, {ex_r} }}\n{buffer}\n]\n"
        #print(f"\\begin{{align*}}\n\\left\\{{ {ex_l}, {ex_r} \\right\\}} &= {c.as_expectation_values()}\n\\end{{align*}}\n")

        d.normalOrder()
        d.sortByCoefficient()
        buffer = d.as_data()
        M += f"[ # {{ {ex_l}, [ H, {ex_r} ]}}\n{buffer}\n]\n"

with open("data/commuting_N.txt", "w") as f:
    f.write(N[:-1])

with open("data/commuting_M.txt", "w") as f:
    f.write(M[:-1])

#ex_l =  deepcopy(F_K_D)
#ex_r =  deepcopy(F_K)
#
#x = Expression(1, np.array([Term(1, "", np.array([Operator(Momentum("k", False, False), False, False)]))]))
#commuted_with_H = commute(H, x)
#commuted_with_H.normalOrder()
#commuted_with_H.sortByCoefficient()
#print(f"\\begin{{align*}}\n {ex_l}, \\left[ H, {x} \\right]  &= {commuted_with_H}\n\\end{{align*}}\n")
#
#c = product(ex_l, ex_r)
#c.normalOrder()
#print(f"\\begin{{align*}}\n\\left\\{{ {ex_l}, \\left[ H, {ex_r} \\right] \\right\\}} &= {c.as_expectation_values()}\n\\end{{align*}}\n")