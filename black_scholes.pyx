from math import e, pi
inv_sqrt_2pi = (2*pi) ** (-0.5)

cdef extern from "math.h":
    double log(double x)
    
cdef extern from "math.h":
    double sqrt(double x)

cdef double CDF(double z):
    ''' Sums the integral of the talor series for e ^ (-0.5 * z**2) '''
    cdef double total, item, z_sq
    cdef int m, k

    if z < -6:
        return 0
    if z > 6:
        return 1
    total = 0.0
    z_sq = z * z
    for k in range(0, 100, 2):
        item = 1
        for m in range(1, k+1):
            item *= (2*m)
        item = z ** (2*k+1) / ((2*k+1) * item)
        item = item * (1 -  ((2*k+1) * z_sq) / ((2*k+2)*(2*k+3)))
        if k > 50 and item == 0:
            break
        total += item

    return 0.5 + inv_sqrt_2pi * total

cdef double pdf(double z):
    return inv_sqrt_2pi * e ** (-0.5*z*z)
    
def std_normal(double z):
    return CDF(z)

cdef class Option:
    ''' Option(call_flag, S, K, vol, q, r, t)'''
    cdef double _S, _K, _vol, _q, _r, _t
    cdef int _call_flag
    
    property call_flag:
        def __get__(self):
            return self._call_flag
    
    property S:
        def __get__(self):
            return self._S
    
    property K:
        def __get__(self):
            return self._K
            
    property vol:
        def __get__(self):
            return self._vol
            
    property q:
        def __get__(self):
            return self._q
            
    property r:
        def __get__(self):
            return self._r
            
    property t:
        def __get__(self):
            return self._t

    def __cinit__(object self, int call_flag, double S, \
        double K, double vol, double q, double r, double t):
        self._call_flag = call_flag
        self._S, self._K, self._vol, self._q, self._r, self._t = \
        S, K, vol, q, r, t
        
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
        
    def __len__(self):
        return 7

    def __hash__(self):
        cdef int value
        cdef double item
        value = 0x345678
        for item in self:
            value = (1000003 * value) ^ hash(item)
        value ^= len(self)
        if value == -1:
            value = -2
        return value

    cdef double d(self, int i):
        cdef int y
        cdef double numerator, denominator

        y = -2 * i + 3 #y = lambda i: {1:1, 2:-1}[i]
        numerator = log(self.S / self.K) + \
        (self.r - self.q + y * 0.5 * self.vol ** 2) * self.t

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
        return e ** (- self.q * self.t)

    cdef double ert(self):
        return e ** (- self.r * self.t)

    property price:
        def __get__(self):
            if self.call_flag:
                return self.S * self.eqt() * CDF(self.d(1)) - self.K * self.ert() * CDF(self.d(2))
            else:
                return self.K * self.ert() * CDF(-self.d(2)) - self.S * self.eqt() * CDF(-self.d(1))

    property delta:
        def __get__(self):
            if self.call_flag:
                return self.eqt() * CDF(self.d(1))
            else:
                return - self.eqt() * CDF(-self.d(1))

    property gamma:
        def __get__(self):
            return self.eqt() * pdf(self.d(1)) / (self.S * self.vol * sqrt(self.t))

    property vega:
        def __get__(self):
            return self.S * self.eqt() * pdf(self.d(1)) * sqrt(self.t)
            
    property rho:
        def __get__(self):
            cdef int y = 2 * self.call_flag - 1
            return y * self.K * self.t * self.ert() * CDF(y * self.d(2))
            
    property theta:
        def __get__(self):
            cdef int y = 2 * self.call_flag - 1
            cdef double d1, d2
            cdef double part_a, part_b, part_c

            d1 = self.d(1)
            d2 = self.d(2)
            part_a = - self.eqt() * self.S * pdf(d1) * self.vol / (2 * sqrt(self.t))
            part_b = - y * self.r * self.K * self.ert() * CDF(y * d2)
            part_c = y * self.q * self.S * self.eqt() * CDF(y * d1)
            return part_a + part_b + part_c