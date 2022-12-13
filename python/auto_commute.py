from dataclasses import dataclass

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

    def __eq__(self, other):
        if isinstance(other, Momentum):
            if(self.value != other.value):
                return False
            if(self.positve != other.positve):
                return False
            if(self.add_Q != other.add_Q):
                return False
            return True
        return False

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

left = Operator(Momentum("k", False, False), False, False)
right = Operator(Momentum("k", True, False), True, False)

from copy import deepcopy
commutator = f"\\left[H, {left} {right}\\right] &= "
swapped = False
if(left.daggered == right.daggered): 
    commutator += f"\\epsilon_{{{left.momentum}}} {right}{left} - \\epsilon_{{{right.momentum}}} {left}{right} \\\\"
    ##########
    buffer_l = deepcopy(left)
    buffer_l.additional_Q()
    buffer_r = deepcopy(right)
    buffer_r.additional_Q()
    commutator += f"\n &+ \\Delta_\\text{{CDW}} \\left( {right} {buffer_l} - {left} {buffer_r} \\right) \\\\"
    ##########
    # ensure that spin down is upfront
    if(left.spin):
        b = right
        right = deepcopy(left)
        left = deepcopy(b)
        swapped = True
        commutator += "\n &- "
    else:
        commutator += "\n &+ "
    buffer_l = deepcopy(left)
    buffer_l.make_sc()
    buffer_r = deepcopy(right)
    buffer_r.make_sc()
    commutator += f"\\Delta_\\text{{SC}} \\left( {buffer_l} {right} + {buffer_r} {left}"
    if(left.momentum == buffer_r.momentum):
        commutator += " - 1"
    commutator += " \\right) \\\\"
    ##########
    buffer_l.additional_Q()
    buffer_r.additional_Q()
    if(swapped):
        commutator += "\n &- "
    else:
        commutator += "\n &+ "
    commutator += "\\Delta_\\eta"
    if(left.daggered):
        commutator += "^*"
    commutator += f" \\left( {buffer_l} {right} + {buffer_r} {left}"
    if(left.momentum == buffer_r.momentum):
        commutator += " - 1"
    commutator += " \\right)"

    if(left.daggered):
        commutator = "-\\left[ " + commutator + " \\right]"
else:
    commutator = f"\\left[H, {left} {right}\\right] &= "


print(commutator)
