from dataclasses import dataclass
from copy import deepcopy
import numpy as np

@dataclass
class Momentum:
    value: str
    positve: bool
    add_Q: bool

    def __str__(self):
        ret = ""
        if(not self.positve):
            ret += "-"
        ret += self.value
        if(self.add_Q):
            if(self.positve):
                ret += "+"
            else:
                ret += "-"
            ret += "Q"
        return ret

    def additional_Q(self):
        self.add_Q = not self.add_Q

@dataclass
class Operator:
    momentum: Momentum
    spin: bool
    daggered: bool

    def __str__(self):
        ret = "c_{" + str(self.momentum)
        if(self.spin):
            ret += "\\uparrow"
        else:
            ret += "\\downarrow"
        ret += "}"
        if(self.daggered):
            ret += "^\\dagger"
        return ret

    def additional_Q(self):
        self.momentum.additional_Q()
    
    def make_sc(self):
        self.spin = not self.spin
        self.momentum.positve = not self.momentum.positve
        self.daggered = not self.daggered

    def spin_as_string(self):
        if(self.spin):
            return "\\uparrow"
        return "\\downarrow"

@dataclass
class Term:
    prefactor: int
    coefficient: str
    operators: np.ndarray

@dataclass
class Expression:
    global_prefactor: int
    terms: np.ndarray
    def append(self, values):
        np.append(self.terms, values)

def sync_eps(momentum: Momentum, base=1):
    if(momentum.add_Q):
        return [-1*base, f"\\epsilon_{momentum.value}"]
    return [1*base, f"\\epsilon_{momentum.value}"]

def str_duo(l, r):
    if(l.daggered and not r.daggered):
        if(l.spin == r.spin):
            if(l.momentum.value == r.momentum.value and l.momentum.positve == r.momentum.positve):
                if(l.momentum.add_Q == r.momentum.add_Q):
                    ret = f"n_{{{l.momentum}{l.spin_as_string()}}}"
                else:
                    if(l.momentum.add_Q):
                        ret = f"g_{{{r.momentum}{l.spin_as_string()}}}^\\dagger"
                    else:
                        ret = f"g_{{{l.momentum}{l.spin_as_string()}}}"
                    
                return ret
    
    else:
        if(l.spin != r.spin):
            if(l.momentum.value == r.momentum.value and l.momentum.positve != r.momentum.positve):
                if(not r.daggered):
                    if(r.spin):
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"f_{{{r.momentum}}}"
                        else:
                            ret = f"\\eta_{{{r.momentum}}}"
                    else:
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"(- f_{{{l.momentum}}})"
                        else:
                            ret = f"(- \\eta_{{{l.momentum}}})"
                else:
                    if(not r.spin):
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"f_{{{l.momentum}}}"
                        else:
                            ret = f"\\eta_{{{l.momentum}}}"
                    else:
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"(- f_{{{r.momentum}}})"
                        else:
                            ret = f"(- \\eta_{{{r.momentum}}})"
                    ret += "^\\dagger"
                return ret
    
    return f"{l} {r}"

def dagger_it(src):
    src[0], src[1] = src[1], src[0]
    src[0].daggered = not src[0].daggered
    src[1].daggered = not src[1].daggered

N_U = [Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, False), True,  False)]
N_D = [Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, False), False, False)]

SC   = [Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True, False)]
SC_D = [Operator(Momentum("k", False, False), False, False), Operator(Momentum("k", True, False), True, False)]
dagger_it(SC_D)

CDW_U   = [Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)]
CDW_U_D = [Operator(Momentum("k", True, False), True,  True), Operator(Momentum("k", True, True), True,  False)]
CDW_D   = [Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, True), False, False)]
CDW_D_D = [Operator(Momentum("k", True, False), False, True), Operator(Momentum("k", True, True), False, False)]
dagger_it(CDW_U_D)
dagger_it(CDW_D_D)

ETA   = [Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)]
ETA_D = [Operator(Momentum("k", False, True), False, False), Operator(Momentum("k", True, False), True, False)]
dagger_it(ETA_D)

left  = ETA[0]
right = ETA[1]

comm_ex = Expression(1, np.array([], dtype=Term))

commutator = f"\\left[H, {str_duo(left, right)}\\right] &= "
swapped = False
sign = ["+", "-"]
if(left.daggered == right.daggered):
    if(left.daggered):
        sign = ["-", "+"]
        commutator += "- "
        comm_ex.global_prefactor = -1
    
    comm_ex.append([Term(*sync_eps(left.momentum), np.array([right, left])), Term(*sync_eps(left.momentum, -1), np.array([left, right]))])
    
    commutator += f"\\epsilon_{{{left.momentum}}} {str_duo(right, left)} {sign[1]} \\epsilon_{{{right.momentum}}} {str_duo(left, right)} \\\\"
    ##########
    buffer_l = deepcopy(left)
    buffer_l.additional_Q()
    buffer_r = deepcopy(right)
    buffer_r.additional_Q()
    comm_ex.append([Term(1, "\\Delta_\\text{{CDW}}", np.array([right, buffer_l])), Term(-1, "\\Delta_\\text{{CDW}}", np.array([right, buffer_l]))])

    commutator += f"\n &{sign[0]} \\Delta_\\text{{CDW}} \\left( {str_duo(right, buffer_l)} - {str_duo(left, buffer_r)} \\right) \\\\"
    ##########
    # ensure that spin down is upfront
    if(left.spin):
        b = right
        right = deepcopy(left)
        left = deepcopy(b)
        swapped = True
        commutator += f"\n &{sign[1]} "
    else:
        commutator += f"\n &{sign[0]} "
    buffer_l = deepcopy(left)
    buffer_l.make_sc()
    buffer_r = deepcopy(right)
    buffer_r.make_sc()
    commutator += f"\\Delta_\\text{{SC}} \\left( {str_duo(buffer_l, right)} + {str_duo(buffer_r, left)}"
    if(left.momentum == buffer_r.momentum):
        commutator += " - 1"
    commutator += " \\right) \\\\"
    ##########
    buffer_l.additional_Q()
    buffer_r.additional_Q()
    if(swapped):
        commutator += f"\n &{sign[1]} "
    else:
        commutator += f"\n &{sign[0]} "
    commutator += "\\Delta_\\eta"
    if(left.daggered):
        commutator += "^*"
    commutator += f" \\left( {str_duo(buffer_l, right)} + {str_duo(buffer_r, left)}"
    if(left.momentum == buffer_r.momentum):
        commutator += " - 1"
    commutator += " \\right)"
else:
    #ensure that we have a normal ordered term
    if(right.daggered):
        b = right
        right = deepcopy(left)
        left = deepcopy(b)
        sign = ["-", "+"]
        if(left.momentum == right.momentum and left.spin == right.spin):
            commutator += "1 "
        commutator += "- "
    commutator += f"\\epsilon_{{{left.momentum}}} {str_duo(left, right)} {sign[1]} \\epsilon_{{{right.momentum}}} {str_duo(left, right)} \\\\"
    ##########
    buffer_l = deepcopy(left)
    buffer_l.additional_Q()
    buffer_r = deepcopy(right)
    buffer_r.additional_Q()
    commutator += f"\n &{sign[0]} \\Delta_\\text{{CDW}} \\left( {str_duo(buffer_l, right)} - {str_duo(left, buffer_r)} \\right) \\\\"
    ##########
    buffer_l = deepcopy(left)
    buffer_r = deepcopy(right)
    buffer_l.make_sc()
    buffer_r.make_sc()
    commutator += f"\n &{sign[0]} \\Delta_\\text{{SC}} \\left( "
    if(left.spin):
        commutator += f"{str_duo(buffer_l, right)}"
    else:
        commutator += f"{str_duo(right, buffer_l)}"
    commutator += f" {sign[1]} "
    if(right.spin):
        commutator += f"{str_duo(left, buffer_r)}"
    else:
        commutator += f"{str_duo(buffer_r, left)}"
    commutator += f" \\right) \\\\"
    ##########
    buffer_l.additional_Q()
    buffer_r.additional_Q()
    commutator += f"\n &{sign[0]} \\Delta_\\eta^* "
    if(left.spin):
        commutator += f"{str_duo(buffer_l, right)}"
    else:
        commutator += f"{str_duo(right, buffer_l)}"
    commutator += f" {sign[1]} \\Delta_\\eta "
    if(right.spin):
        commutator += f"{str_duo(left, buffer_r)}"
    else:
        commutator += f"{str_duo(buffer_r, left)}"

print(commutator)
