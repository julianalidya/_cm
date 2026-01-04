from __future__ import annotations

class FiniteField:
    def __init__(self, p: int):
        if not isinstance(p, int) or p <= 1:
            raise ValueError("p must be an integer > 1")
        self.p = p

    def __call__(self, value: int):
        return FFElem(value, self)

    def elements(self):
        return [self(i) for i in range(self.p)]

    def nonzero_elements(self):
        return [self(i) for i in range(1, self.p)]

    def __repr__(self):
        return f"GF({self.p})"


class FFElem:
    __slots__ = ("v", "F")

    def __init__(self, value: int, field: FiniteField):
        self.F = field
        self.v = value % field.p

    def _coerce(self, other):
        if isinstance(other, FFElem):
            if other.F.p != self.F.p:
                raise TypeError("Different fields")
            return other
        if isinstance(other, int):
            return FFElem(other, self.F)
        return NotImplemented

    def __add__(self, other):
        other = self._coerce(other)
        return FFElem(self.v + other.v, self.F)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        other = self._coerce(other)
        return FFElem(self.v - other.v, self.F)

    def __rsub__(self, other):
        other = self._coerce(other)
        return other.__sub__(self)

    def __neg__(self):
        return FFElem(-self.v, self.F)

    def __mul__(self, other):
        other = self._coerce(other)
        return FFElem(self.v * other.v, self.F)

    def __rmul__(self, other):
        return self.__mul__(other)

    def inv(self):
        if self.v == 0:
            raise ZeroDivisionError
        p = self.F.p
        t, new_t = 0, 1
        r, new_r = p, self.v
        while new_r:
            q = r // new_r
            t, new_t = new_t, t - q * new_t
            r, new_r = new_r, r - q * new_r
        return FFElem(t, self.F)

    def __truediv__(self, other):
        other = self._coerce(other)
        return self * other.inv()

    def __rtruediv__(self, other):
        other = self._coerce(other)
        return other.__truediv__(self)

    def __pow__(self, n: int):
        if n < 0:
            return (self.inv()) ** (-n)
        r = FFElem(1, self.F)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def __eq__(self, other):
        if isinstance(other, FFElem):
            return self.F.p == other.F.p and self.v == other.v
        if isinstance(other, int):
            return self.v == other % self.F.p
        return False

    def __hash__(self):
        return hash((self.F.p, self.v))

    def __int__(self):
        return self.v

    def __repr__(self):
        return f"{self.v}"


class AddGroup:
    def __init__(self, field):
        self.F = field

    def elements(self):
        return self.F.elements()

    def op(self, a, b):
        return a + b

    def identity(self):
        return self.F(0)

    def inv(self, a):
        return -a


class MulGroup:
    def __init__(self, field):
        self.F = field

    def elements(self):
        return self.F.nonzero_elements()

    def op(self, a, b):
        return a * b

    def identity(self):
        return self.F(1)

    def inv(self, a):
        return a.inv()


def check_group(G):
    elems = G.elements()
    e = G.identity()
    for a in elems:
        for b in elems:
            if G.op(a, b) not in elems:
                return False
    for a in elems:
        if G.op(a, e) != a or G.op(e, a) != a:
            return False
    for a in elems:
        if G.op(a, G.inv(a)) != e or G.op(G.inv(a), a) != e:
            return False
    for a in elems:
        for b in elems:
            for c in elems:
                if G.op(G.op(a, b), c) != G.op(a, G.op(b, c)):
                    return False
    return True


def check_distributivity(field):
    elems = field.elements()
    for a in elems:
        for b in elems:
            for c in elems:
                if a * (b + c) != (a * b) + (a * c):
                    return False
    return True


if __name__ == "__main__":
    F = FiniteField(7)
    a = F(3)
    b = F(5)
    print(a + b, a - b, a * b, a / b)
    print(check_group(AddGroup(F)))
    print(check_group(MulGroup(F)))
    print(check_distributivity(F))
