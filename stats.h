double const inv_sq_2pi = 0.3989422804014327;

double CDF(double z) {
    /* Takes the integral of the taylor series for
       e ** (-1/2 * z**2) and sums the first 100
       terms of the series.
    f(n, z) = (((-1)**n) * z ** (2n+1)) / ((2n+1)(*2**n)(n!))
    CDF(z) = f(0,z) + f(1,z) + ... + f(99, z)*/
    int k;
    double m, total, item, z2, z4, a, b;
    
    if (z < -6) {
        return 0;
    }
    if (z >  6) {
        return 1;
    }

    m = 1;         // m(k) == (2**k)/factorial(k)
    b = z;         // b(k) == z ** (2*k + 1)
    z2 = z * z;    // cache of z squared
    z4 = z2 * z2;  // cache of z to the 4th
    total = 0;

    /* The series is conditionally convergent
    we group the terms into pairs, one positive,
    one negative, so that we avoid excessive rounding
    errors and overflow errors.*/
    for (k=0; k<100; k+=2) {
        a = 2*k + 1;
        item = b / (a*m);
        item *= (1 - (a*z2)/((a+1)*(a+2)));
        total += item;
        m *= (4*(k+1)*(k+2));
        b *= z4;
    }
    
    return 0.5 + inv_sq_2pi * total;
}