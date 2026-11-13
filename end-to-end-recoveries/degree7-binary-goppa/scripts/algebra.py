"""Small exact arithmetic used only for public-input preparation and tests."""
from math import comb


def rref(rows, width):
    rows = list(rows)
    pivots = []
    for col in range(width):
        at = len(pivots)
        pivot = next((i for i in range(at, len(rows)) if rows[i] >> col & 1), None)
        if pivot is None:
            continue
        rows[at], rows[pivot] = rows[pivot], rows[at]
        for i in range(len(rows)):
            if i != at and rows[i] >> col & 1:
                rows[i] ^= rows[at]
        pivots.append(col)
    return rows[:len(pivots)], pivots


def kernel(rows, width):
    reduced, pivots = rref(rows, width)
    free = [j for j in range(width) if j not in pivots]
    return [(1 << j) | sum(((reduced[i] >> j) & 1) << p
                          for i, p in enumerate(pivots)) for j in free]


def restrict(row, indices):
    return sum(((row >> j) & 1) << i for i, j in enumerate(indices))


def columns(rows, width):
    return [sum(((row >> j) & 1) << i for i, row in enumerate(rows))
            for j in range(width)]


def retained_orders(degree, multiplicity):
    keep = []
    for u in reversed(range(multiplicity)):
        if not any(comb(degree-u, v-u) & 1 for v in keep):
            keep.append(u)
    return sorted(keep)


def choose(n, r):
    return comb(n, r) if 0 <= r <= n else 0


class Field:
    def __init__(self, m, modulus):
        self.m, self.modulus, self.q = m, modulus, 1 << m

    def mul(self, a, b):
        result = 0
        while b:
            if b & 1:
                result ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return result

    def power(self, a, exponent):
        value = 1
        while exponent:
            if exponent & 1:
                value = self.mul(value, a)
            a = self.mul(a, a)
            exponent >>= 1
        return value

    def inv(self, a):
        if not a:
            raise ValueError('zero denominator')
        result = self.power(a, self.q-2)
        if self.mul(a, result) != 1:
            raise ValueError('invalid field arithmetic')
        return result

    @staticmethod
    def trim(a):
        a = list(a)
        while a and not a[-1]:
            a.pop()
        return a

    def add(self, a, b):
        result = list(a) + [0]*max(0, len(b)-len(a))
        for i, value in enumerate(b):
            result[i] ^= value
        return self.trim(result)

    def rem(self, a, b):
        a, b = self.trim(a), self.trim(b)
        inverse = self.inv(b[-1])
        while len(a) >= len(b):
            off, scale = len(a)-len(b), self.mul(a[-1], inverse)
            for i, value in enumerate(b):
                a[off+i] ^= self.mul(scale, value)
            a = self.trim(a)
        return a

    def gcd(self, a, b):
        while b:
            a, b = b, self.rem(a, b)
        return a

    def square_mod(self, a, g):
        a = self.trim(a)
        if not a:
            return []
        result = [0]*(2*len(a)-1)
        for i, value in enumerate(a):
            result[2*i] = self.mul(value, value)
        return self.rem(result, g)

    def irreducible(self, g):
        degree = len(g)-1
        divisors = [p for p in range(2, degree+1)
                    if degree % p == 0 and all(p % d for d in range(2, p))]
        h = [0, 1]
        for r in range(1, degree+1):
            for _ in range(self.m):
                h = self.square_mod(h, g)
            if r in [degree//p for p in divisors]:
                if len(self.gcd(g, self.add(h, [0, 1]))) != 1:
                    return False
        return h == [0, 1]

    def evaluate(self, poly, x):
        result = 0
        for value in reversed(poly):
            result = self.mul(result, x) ^ value
        return result

    def binary_parity(self, support, multipliers, degree):
        rows = [0]*(self.m*degree)
        for j, (a, value) in enumerate(zip(support, multipliers)):
            for e in range(degree):
                for bit in range(self.m):
                    rows[e*self.m+bit] |= ((value >> bit) & 1) << j
                value = self.mul(value, a)
        return rref(rows, len(support))[0]
