from dataclasses import dataclass
from copy import deepcopy
import numpy as np

def str_duo(l, r):
    if(l.daggered and not r.daggered):
        if(l.spin == r.spin):
            if(l.momentum.value == r.momentum.value and l.momentum.positive == r.momentum.positive):
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
            if(l.momentum.value == r.momentum.value and l.momentum.positive != r.momentum.positive):
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

@dataclass
class Momentum:
    value: str
    positive: bool
    add_Q: bool

    def __str__(self):
        ret = ""
        if(not self.positive):
            ret += "-"
        ret += self.value
        if(self.add_Q):
            if(self.positive):
                ret += "+"
            else:
                ret += "-"
            ret += "Q"
        return ret

    def additional_Q(self):
        self.add_Q = not self.add_Q

    def copy(self):
        return Momentum(self.value, self.positive, self.add_Q)

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
        self.momentum.positive = not self.momentum.positive
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

    def isIdentity(self):
        return self.operators.size == 0
    
    def append(self, values):
        self.operators = np.append(self.operators, values)

    def __eq__(self, other):
        if(self.coefficient != other.coefficient):
            return False
        if(self.operators.size != other.operators.size):
            return False
        for i in range(0, self.operators.size):
            if(self.operators[i] != other.operators[i]):
                return False
        return True

    def __str__(self):
        ret = ""
        if(self.prefactor < 0):
            ret = "- "
        if(abs(self.prefactor) != 1):
            ret += f"{abs(self.prefactor)} "
        if(self.coefficient != ""):
            ret += f"{self.coefficient} \\cdot "
        if(self.isIdentity()):
            return ret + "\\mathbb{1}"
        if(self.operators.size == 2):
            ret += str_duo(self.operators[0], self.operators[1])
            return ret
        for o in self.operators:
            ret += f"{o} "
        return ret

    def to_string_no_coefficient(self):
        ret = ""
        if(self.prefactor < 0):
            ret = "- "
        if(abs(self.prefactor) != 1):
            ret += f"{abs(self.prefactor)} "
        if(self.isIdentity()):
            return ret + "\\mathbb{1}"
        if(self.operators.size == 2):
            ret += str_duo(self.operators[0], self.operators[1])
            return ret
        for o in self.operators:
            ret += f"{o} "
        return ret

@dataclass
class Expression:
    global_prefactor: int
    terms: np.ndarray

    def __str__(self):
        first = True
        ret = ""
        if(self.global_prefactor != 1):
            ret = f"{self.global_prefactor} \\cdot \\Big["
        for i, t in enumerate(self.terms):
            if(i > 0):
                if(self.terms[i].coefficient != self.terms[i - 1].coefficient):
                    ret += "\\right) \\\\\n&"
                    first = True
            if(first):
                first = False
                ret += f"+{t.coefficient} \\left("
            if(t.prefactor >= 0):
                ret += "+ "
            ret += f"{t.to_string_no_coefficient()}"
        ret += "\\right) "
        if(self.global_prefactor != 1):
            ret += "\\Big]"
        return ret

    def append(self, values):
        self.terms = np.append(self.terms, values)

    def sortByCoefficient(self):
        i = 0
        while i < self.terms.size:
            start = i + 1
            while(start < self.terms.size and self.terms[i].coefficient == self.terms[start].coefficient):
                start += 1
            for j in range(start + 1, self.terms.size):
                if(self.terms[i].coefficient == self.terms[j].coefficient):
                    self.terms[start], self.terms[j] = self.terms[j], self.terms[start]
                    i += 1
                    start += 1
                j += 1
            i += 1


    def sortTerms(self):
        for t in self.terms:
            for i in range(0, t.operators.size):
                for j in range(i + 1, t.operators.size):
                    if(t.operators[i].daggered == t.operators[j].daggered):
                        if(not t.operators[i].daggered):
                            if(t.operators[i].spin and not t.operators[j].spin):
                                t.operators[i], t.operators[j] = t.operators[j], t.operators[i]
                                if((j - i) % 2 == 1):
                                    t.prefactor *= -1
                        else:
                            if(not t.operators[i].spin and t.operators[j].spin):
                                t.operators[i], t.operators[j] = t.operators[j], t.operators[i]
                                if((j - i) % 2 == 1):
                                    t.prefactor *= -1
            for i in range(0, t.operators.size):
                for j in range(i + 1, t.operators.size):
                    if(t.operators[i].daggered == t.operators[j].daggered and t.operators[i].spin == t.operators[j].spin):
                        if(t.operators[i].momentum.add_Q and not t.operators[j].momentum.add_Q):
                            t.operators[i], t.operators[j] = t.operators[j], t.operators[i]
                            if((j - i) % 2 == 1):
                                t.prefactor *= -1
            
            for i in range(0, t.operators.size):
                for j in range(i + 1, t.operators.size):
                    if(t.operators[i].daggered == t.operators[j].daggered and t.operators[i].spin == t.operators[j].spin and t.operators[i].momentum.add_Q == t.operators[j].momentum.add_Q):
                        if(not t.operators[i].momentum.positive and t.operators[j].momentum.positive):
                            t.operators[i], t.operators[j] = t.operators[j], t.operators[i]
                            if((j - i) % 2 == 1):
                                t.prefactor *= -1

    def normalOrder(self):
        t = 0
        while t < self.terms.size:
            n = self.terms[t].operators.size
            while(n > 1):
                new_n = 0
                for i in range(1, self.terms[t].operators.size):
                    if(self.terms[t].operators[i - 1] == self.terms[t].operators[i]):
                        self.terms = np.delete(self.terms, t)
                        t -= 1
                        n = 0
                        break
                    if(not self.terms[t].operators[i - 1].daggered and self.terms[t].operators[i].daggered):
                        new_n = i
                        self.terms[t].operators[i], self.terms[t].operators[i - 1] = self.terms[t].operators[i - 1], self.terms[t].operators[i]
                        self.terms[t].prefactor *= -1
                        if(self.terms[t].operators[i].momentum == self.terms[t].operators[i - 1].momentum and self.terms[t].operators[i].spin == self.terms[t].operators[i - 1].spin):
                            self.append([ Term(-self.terms[t].prefactor, self.terms[t].coefficient, np.delete(self.terms[t].operators, [i - 1, i])) ])
                n = new_n
            t += 1

        self.sortTerms()

        i = 0
        while i < self.terms.size:
            skip = False
            for k in range(1, self.terms[i].operators.size):
                if(self.terms[i].operators[k] == self.terms[i].operators[k - 1]):
                    self.terms = np.delete(self.terms, i)
                    skip = True
                    break
                    
            if(not skip):
                j = i + 1
                while j < self.terms.size:   
                    if(self.terms[i] == self.terms[j]):
                        self.terms[i].prefactor += self.terms[j].prefactor
                        self.terms = np.delete(self.terms, j)
                    else:
                        j += 1
                if(self.terms[i].prefactor == 0):
                    self.terms = np.delete(self.terms, i)
                else:
                    i += 1

    def as_expectation_values(self) -> str:
        ret = ""
        for i, t in enumerate(self.terms):
            if(i > 0):
                ret += "\\\\\n&"
            if(t.prefactor >= 0):
                ret += "+ "
            else:
                ret += "- "
            if(abs(t.prefactor) != 1):
                ret += f"{abs(t.prefactor)} "
            if(t.coefficient != ""):
                ret += f"{t.coefficient} \\cdot "
            ret += t.wick()

        return ret
            

def sync_eps(momentum: Momentum, base=1):
    if(momentum.add_Q):
        return [-1*base, f"\\epsilon_{momentum.value}"]
    return [1*base, f"\\epsilon_{momentum.value}"]

def dagger_it(src):
    src[0], src[1] = src[1], src[0]
    src[0].daggered = not src[0].daggered
    src[1].daggered = not src[1].daggered

def anti_commute(l: Expression, r: Expression):
    commuted = Expression(l.global_prefactor * r.global_prefactor, np.array([], dtype=Term))
    for lt in l.terms:
        for rt in r.terms:
            commuted.append([ Term(lt.prefactor * rt.prefactor, lt.coefficient + rt.coefficient, np.append(lt.operators, rt.operators)) ])
            commuted.append([ Term(lt.prefactor * rt.prefactor, lt.coefficient + rt.coefficient, np.append(rt.operators, lt.operators)) ])

    return commuted

def commute(l: Expression, r: Expression):
    commuted = Expression(l.global_prefactor * r.global_prefactor, np.array([], dtype=Term))
    for lt in l.terms:
        for rt in r.terms:
            commuted.append([ Term(lt.prefactor * rt.prefactor, lt.coefficient + rt.coefficient, np.append(lt.operators, rt.operators)) ])
            commuted.append([ Term(-lt.prefactor * rt.prefactor, lt.coefficient + rt.coefficient, np.append(rt.operators, lt.operators)) ])

    return commuted