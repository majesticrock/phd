from copy import deepcopy
from methods_for_commute import *

N_U  = np.array([Operator(Momentum("k", True, False), True,  True),   Operator(Momentum("k", True, False), True,  False)])
N_D  = np.array([Operator(Momentum("k", True, False), False, True),   Operator(Momentum("k", True, False), False, False)])
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
left  = SC[0]
right = SC[1]

commuted_with_H = Expression(1, np.array([], dtype=Term))
if(left.daggered == right.daggered):
    if(left.daggered):
        if(not left.spin):
            b = right
            right = deepcopy(left)
            left = deepcopy(b)
        else:
            commuted_with_H.global_prefactor = -1
    elif(left.spin):
        b = right
        right = deepcopy(left)
        left = deepcopy(b)
        commuted_with_H.global_prefactor = -1

    commuted_with_H.append([ Term(*sync_eps(left.momentum), np.array([right, left])), Term(*sync_eps(left.momentum, -1), np.array([left, right])) ])
    ##########
    buffer_l = deepcopy(left)
    buffer_l.additional_Q()
    buffer_r = deepcopy(right)
    buffer_r.additional_Q()
    commuted_with_H.append([ Term(1, "\\Delta_\\text{CDW}", np.array([right, buffer_l])), Term(-1, "\\Delta_\\text{CDW}", np.array([left, buffer_r])) ])
    ##########
    buffer_l = deepcopy(left)
    buffer_l.make_sc()
    buffer_r = deepcopy(right)
    buffer_r.make_sc()
    commuted_with_H.append([ Term(1, "\\Delta_\\text{SC}", np.array([buffer_l, right])), Term(1, "\\Delta_\\text{SC}", np.array([buffer_r, left])) ])
    if(left.momentum == buffer_r.momentum):
        commuted_with_H.append(Term(-1, "\\Delta_\\text{SC}", np.array([], dtype=Operator)))
    ##########
    buffer_l = deepcopy(buffer_l)
    buffer_r = deepcopy(buffer_r)
    buffer_l.additional_Q()
    buffer_r.additional_Q()
    if(left.daggered):
        commuted_with_H.append([ Term(1, "\\Delta_\\eta^*", np.array([buffer_l, right])), Term(1, "\\Delta_\\eta^*", np.array([buffer_r, left])) ])
        if(left.momentum == buffer_r.momentum):
            commuted_with_H.append(Term(-1, "\\Delta_eta^*", np.array([], dtype=Operator)))
    else:
        commuted_with_H.append([ Term(1, "\\Delta_\\eta", np.array([buffer_l, right])), Term(1, "\\Delta_\\eta", np.array([buffer_r, left])) ])
        if(left.momentum == buffer_r.momentum):
            commuted_with_H.append(Term(-1, "\\Delta_eta", np.array([], dtype=Operator)))
else:
    #ensure that we have a normal ordered term
    if(right.daggered):
        b = right
        right = deepcopy(left)
        left = deepcopy(b)
        commuted_with_H.global_prefactor = -1
        if(left.momentum == right.momentum and left.spin == right.spin):
            commuted_with_H.append([Term(-1, "", np.array([], dtype=Operator))])

    commuted_with_H.append([ Term(*sync_eps(left.momentum), np.array([left, right])), Term(*sync_eps(left.momentum, -1), np.array([left, right])) ])
    ##########
    buffer_l = deepcopy(left)
    buffer_l.additional_Q()
    buffer_r = deepcopy(right)
    buffer_r.additional_Q()
    commuted_with_H.append([ Term(1, "\\Delta_\\text{CDW}", np.array([buffer_l, right])), Term(-1, "\\Delta_\\text{CDW}", np.array([left, buffer_r])) ])
    ##########
    buffer_l = deepcopy(left)
    buffer_r = deepcopy(right)
    buffer_l.make_sc()
    buffer_r.make_sc()

    if(left.spin):
        commuted_with_H.append([ Term(1, "\\Delta_\\text{SC}", np.array([buffer_l, right])) ])
    else:
        commuted_with_H.append([ Term(1, "\\Delta_\\text{SC}", np.array([right, buffer_l])) ])

    if(right.spin):
        commuted_with_H.append([ Term(-1, "\\Delta_\\text{SC}", np.array([left, buffer_r])) ])
    else:
        commuted_with_H.append([ Term(-1, "\\Delta_\\text{SC}", np.array([buffer_r, left])) ])

    ##########
    buffer_l = deepcopy(buffer_l)
    buffer_r = deepcopy(buffer_r)
    buffer_l.additional_Q()
    buffer_r.additional_Q()
    if(left.spin):
        commuted_with_H.append([ Term(1, "\\Delta_\\eta^*", np.array([buffer_l, right])) ])
    else:
        commuted_with_H.append([ Term(1, "\\Delta_\\eta^*", np.array([right, buffer_l])) ])

    if(right.spin):
        commuted_with_H.append([ Term(1, "\\Delta_\\eta", np.array([left, buffer_r])) ])
    else:
        commuted_with_H.append([ Term(1, "\\Delta_\\eta", np.array([buffer_r, left])) ])

commuted_with_H.normalOrder()
print(f"\\left[H, {str_duo(left, right)}\\right] &= {commuted_with_H}\n")

ex_l = Expression(1, np.array([ Term(1, "", SC_D) ]))
ex_r = Expression(1, np.array([ Term(1, "", SC) ]))
c = anti_commmute(ex_l, commuted_with_H)
c.normalOrder()
print(f"\\left\\{{ {ex_l}, \\left[ H, {ex_r} \\right] \\right\\}} &= {c}")