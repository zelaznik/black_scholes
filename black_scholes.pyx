from math import e, pi
inv_sqrt_2pi = (2*pi) ** (-0.5)

cdef extern from "math.h":
    double log(double x)
    
cdef extern from "math.h":
    double sqrt(double x)
    
cdef extern from "math.h":
    double exp(double x)

cdef double pdf(double z):
    return inv_sqrt_2pi * exp(-0.5*z*z)
    
cdef extern from "stats.h":
    double CDF(double z)

def std_normal(z):
    return CDF(z)

cdef class Option:
    ''' Option(call_flag, S, K, vol, q, r, t)'''
    cdef double S, K, vol, q, r, t
    cdef double d1, d2
    cdef int call_flag, y

    def __cinit__(object self, int call_flag, double S, \
        double K, double vol, double q, double r, double t):
        self.call_flag = call_flag
        self.S, self.K, self.vol, self.q, self.r, self.t = \
        S, K, vol, q, r, t
        self.y = (1 if call_flag else -1)

    def __repr__(self):
        cls = ('Call' if self.call_flag else 'Put')
        return '%s(S=%r, K=%r, vol=%r, q=%r, r=%r, t=%r)' % \
        (cls, self.S, self.K, self.vol, self.q, self.r, self.t)
        
    def __iter__(self):
        yield bool(self.call_flag)
        yield self.S
        yield self.K
        yield self.vol
        yield self.q
        yield self.r
        yield self.t

    def __hash__(self):
        cdef int i
        cdef long value
        cdef double item
        value = 0x345678
        for item in self:
            i += 1
            value = (1000003 * value) ^ hash(item)
        value ^= i
        if value == -1:
            value = -2
        return value

    cdef double d(self, int i):
        cdef int a = 3 - 2*i
        cdef double numerator, denominator

        numerator = log(self.S / self.K) + \
        (self.r - self.q + a * 0.5 * self.vol * self.vol) * self.t

        denominator = self.vol * sqrt(self.t)
        if denominator != 0:
            return numerator / denominator

        #If the denominator == 0, it means that
        #the total variance == 0 as well.
        #This means that Nd1 and Nd2 are either
        #both equal to 1 or both equal to 0
        if numerator >= 0:
            return 999.0
        elif numerator < 0:
            return -999.0
        
    cdef double eqt(self):
        return exp(- self.q * self.t)

    cdef double ert(self):
        return exp(- self.r * self.t)

    cpdef double price(self):
        cdef int y = self.y
        return y * (self.S * self.eqt() * CDF(y * self.d(1)) - self.K * self.ert() * CDF(y * self.d(2)))

    cpdef double delta(self):
        cdef int y = self.y
        return y * self.eqt() * CDF(y * self.d(1))

    cpdef double gamma(self):
        return self.eqt() * pdf(self.d(1)) / (self.S * self.vol * sqrt(self.t))

    cpdef double vega(self):
        return self.S * self.eqt() * pdf(self.d(1)) * sqrt(self.t)

    cpdef double rho(self):
        cdef int y = self.y
        return y * self.K * self.t * self.ert() * CDF(y * self.d(2))

    cpdef double theta(self):
        cdef int y = 2 * self.call_flag - 1
        cdef double d1, d2
        cdef double part_a, part_b, part_c

        d1 = self.d(1)
        d2 = self.d(2)
        part_a = - self.eqt() * self.S * pdf(d1) * self.vol / (2 * sqrt(self.t))
        part_b = - y * self.r * self.K * self.ert() * CDF(y * d2)
        part_c = y * self.q * self.S * self.eqt() * CDF(y * d1)
        return part_a + part_b + part_c