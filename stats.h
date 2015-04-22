double const inv_sq_2pi = 0.3989422804014327;

double CDF(double z) {
    int k;
    double m, total, item, z2, z4, a, b;
    
    if (z < -6) {
        return 0;
    }
    if (z >  6) {
        return 1;
    }

    m = 1;
    b = z;
    z2 = z * z;
    z4 = z2 * z2;
    total = 0;
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

double fabs(double z) {
    if (z < 0) {
        return -z;
    }
    return z;
}

double nerf(double xi)
{
    double      ret;
    double      x;
    int         i;
    double a[7] = {.0002765672,.0001520143,.0092705272,.0422820123,.0705230784,1.,.0000430638};

    x=fabs(xi);
    if( x > 10. ) {
        ret = (xi < 0) ? -1. : 1.;
    }
    else {
        ret = a[6];
        for (i = 0; i < 6; i++) {
            ret = a[i] + x * ret;
        }
        for (i = 0; i < 4; i++) {
            ret = ret * ret;
        }
        ret = 1.-1./ret;
        ret = (xi < 0.) ? -ret : ret;
    }
    return ret;
}
