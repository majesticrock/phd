from dataclasses import dataclass
import numpy as np

def str_duo(l, r, as_data=False):
    if(l.daggered and not r.daggered):
        if(l.spin == r.spin):
            if(l.momentum.value == r.momentum.value and l.momentum.positive == r.momentum.positive):
                if(l.momentum.add_Q == r.momentum.add_Q):
                    ret = f"n_{{{l.momentum}{l.spin_as_string()}}}"
                    if(as_data):
                        ret = f"{{n{l.spin_as_string(True)},0,{l.momentum.as_data()}}}"
                else:
                    if(l.momentum.add_Q):
                        ret = f"g_{{{r.momentum}{l.spin_as_string()}}}^\\dagger"
                        if(as_data):
                            ret = f"{{g{l.spin_as_string(True)},1,{r.momentum.as_data()}}}"
                    else:
                        ret = f"g_{{{l.momentum}{l.spin_as_string()}}}"
                        if(as_data):
                            ret = f"{{g{l.spin_as_string(True)},0,{l.momentum.as_data()}}}"         
                return ret
    else:
        if(l.spin != r.spin):
            if(l.momentum.value == r.momentum.value and l.momentum.positive != r.momentum.positive):
                if(not r.daggered):
                    if(r.spin):
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"f_{{{r.momentum}}}"
                            if(as_data):
                                ret = f"{{f,0,{r.momentum.as_data()}}}"
                        else:
                            ret = f"\\eta_{{{r.momentum}}}"
                            if(as_data):
                                ret = f"{{eta,0,{r.momentum.as_data()}}}"
                    else: # Should not come up anymore
                        print("MEH")
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"(- f_{{{l.momentum}}})"
                        else:
                            ret = f"(- \\eta_{{{l.momentum}}})"
                else:
                    if(not r.spin):
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"f_{{{l.momentum}}}"
                            if(as_data):
                                ret = f"{{f,1,{l.momentum.as_data()}}}"
                        else:
                            ret = f"\\eta_{{{l.momentum}}}"
                            if(as_data):
                                ret = f"{{eta,1,{l.momentum.as_data()}}}"
                    else: # Should not come up anymore
                        print("MEH")
                        if(l.momentum.add_Q == r.momentum.add_Q):
                            ret = f"(- f_{{{r.momentum}}})"
                        else:
                            ret = f"(- \\eta_{{{r.momentum}}})"
                    if(not as_data):
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

    def as_data(self) -> str:
        return f"{1 if self.positive else 0},{1 if self.add_Q else 0}"

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

    def spin_as_string(self, as_data=False):
        if(as_data):
            return "1" if self.spin else "0"
        return "\\uparrow" if self.spin else "\\downarrow"

@dataclass
class Term:
    prefactor: float
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

    def to_string_no_coefficient(self) -> str:
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

    def to_string_no_prefactor(self) -> str:
        ret = ""
        if(self.isIdentity()):
            return ret + "\\mathbb{1}"
        if(self.operators.size == 2):
            ret += str_duo(self.operators[0], self.operators[1])
            return ret
        for o in self.operators:
            ret += f"{o} "
        return 

    def wick(self) -> str:
        ret = ""
        if(self.operators.size == 2):
            ret += rf"\langle {self.to_string_no_prefactor()} \rangle"
        else:
            did_something = False
            ret = r"\Big( "
            for i in range(1, self.operators.size):
                # conservation of spin and momentum
                y = Term(1, "", np.array([self.operators[0], self.operators[i]]))
                if(self.operators[0].daggered == self.operators[i].daggered and self.operators[0].spin != self.operators[i].spin and self.operators[0].momentum.positive != self.operators[i].momentum.positive):
                    x = Term(1, "", np.array([], dtype=Operator))
                    for j in range(1, self.operators.size):
                        if(j != i):
                            x.append(self.operators[j])
                    if(i % 2 == 0):
                        ret += " - "
                    else:
                        ret += " + "
                    ret += rf"\langle {y} \rangle"
                    ret += x.wick()
                    did_something = True
                elif(self.operators[0].daggered != self.operators[i].daggered and self.operators[0].spin == self.operators[i].spin and self.operators[0].momentum.positive == self.operators[i].momentum.positive):
                    x = Term(1, "", np.array([], dtype=Operator))
                    for j in range(1, self.operators.size):
                        if(j != i):
                            x.append(self.operators[j])
                    if(i % 2 == 0):
                        ret += " - "
                    else:
                        ret += " + "
                    ret += rf"\langle {y} \rangle"
                    ret += x.wick()
                    did_something = True
            
            if(not did_something):
                ret = ""
            else: 
                ret += r" \Big)"
        return ret

    def coeff_as_data(self) -> str:
        if(self.coefficient == r"\epsilon_k"):
            return "0"
        if(self.coefficient == r"\Delta_\text{CDW}"):
            return "1"
        if(self.coefficient == r"\Delta_\text{SC}"):
            return "2"
        if(self.coefficient == r"\Delta_\eta^*"):
            return "3"
        if(self.coefficient == r"\Delta_\eta"):
            return "4"
        if(self.coefficient == ""):
            return "-1"

    def wick_as_data(self) -> str:
        ret = ""
        if(self.operators.size == 2):
            return str_duo(self.operators[0], self.operators[1], True)
        else:
            did_something = False
            for i in range(1, self.operators.size):
                if(i > 1 and ret != ""):
                    if(ret[-1] == "}"):
                        ret += ","
                y = Term(1, "", np.array([self.operators[0], self.operators[i]]))
                # conservation of spin and momentum
                if(self.operators[0].daggered == self.operators[i].daggered and self.operators[0].spin != self.operators[i].spin and self.operators[0].momentum.positive != self.operators[i].momentum.positive):
                    x = Term(1, "", np.array([], dtype=Operator))
                    for j in range(1, self.operators.size):
                        if(j != i):
                            x.append(self.operators[j])
                    if(i % 2 == 0):
                        ret += "{-,"
                    else:
                        ret += "{+,"
                    ret += str_duo(y.operators[0], y.operators[1], True)
                    ret += ","
                    ret += x.wick_as_data()
                    ret += "}"
                    did_something = True
                elif(self.operators[0].daggered != self.operators[i].daggered and self.operators[0].spin == self.operators[i].spin and self.operators[0].momentum.positive == self.operators[i].momentum.positive):
                    x = Term(1, "", np.array([], dtype=Operator))
                    for j in range(1, self.operators.size):
                        if(j != i):
                            x.append(self.operators[j])
                    if(i % 2 == 0):
                        ret += "{-,"
                    else:
                        ret += "{+,"
                    ret += str_duo(y.operators[0], y.operators[1], True)
                    ret += ","
                    ret += x.wick_as_data()
                    ret += "}"
                    did_something = True
            if(did_something and ret[-1] == ","):
                ret = ret[:len(ret) - 1]
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
        increment_t = True
        while t < self.terms.size:
            n = self.terms[t].operators.size
            while(n > 1):
                new_n = 0
                for i in range(1, self.terms[t].operators.size):
                    if(self.terms[t].operators[i - 1] == self.terms[t].operators[i]):
                        self.terms = np.delete(self.terms, t)
                        increment_t = False
                        new_n = 0
                        break
                    if(not self.terms[t].operators[i - 1].daggered and self.terms[t].operators[i].daggered):
                        new_n = i
                        self.terms[t].operators[i], self.terms[t].operators[i - 1] = self.terms[t].operators[i - 1], self.terms[t].operators[i]
                        self.terms[t].prefactor *= -1
                        if(self.terms[t].operators[i].momentum == self.terms[t].operators[i - 1].momentum and self.terms[t].operators[i].spin == self.terms[t].operators[i - 1].spin):
                            self.append([ Term(-self.terms[t].prefactor, self.terms[t].coefficient, np.delete(self.terms[t].operators, [i - 1, i])) ])
                n = new_n
            if(increment_t):
                t += 1
            else:
                increment_t = True

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
            to_add = ""
            if(i > 0):
                to_add += "\\\\\n&"
            if(t.prefactor >= 0):
                to_add += "+ "
            else:
                to_add += "- "
            if(abs(t.prefactor) != 1):
                to_add += f"{abs(t.prefactor)} "
            if(t.coefficient != ""):
                to_add += f"{t.coefficient} \\cdot "
            buffer = t.wick()
            if(buffer != ""):
                ret += to_add + buffer

        return ret
    
    def as_data(self):
        ret = ""
        for i, t in enumerate(self.terms):
            to_add = ""
            if (i==0):
                to_add += f"{{\n\t{t.coeff_as_data()}"
            elif(t.coefficient != self.terms[i - 1].coefficient):
                to_add += f"\n}}\n{{\n\t{t.coeff_as_data()}"
            buffer = t.wick_as_data()
            if(buffer != ""):
                to_add += f"\n\t{{\n\t\t{t.prefactor}\n\t\t{buffer}\n\t}}"
                ret += to_add
        ret += "\n}"
        if(ret == "\n}"):
            return ""
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