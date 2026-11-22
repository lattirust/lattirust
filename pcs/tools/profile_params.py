"""Exact layout enumeration, and profile-static native f64 parameter recipe."""
from math import ceil, floor, log2, sqrt
from fractions import Fraction

P=2**32-99

def nearest_even_integer(x):
    if floor(x)%2==0: return floor(x)
    return ceil(x)+int(x==floor(x))

def integer_nth_root(x,n):
    a,b=0,x+1
    while b-a>1:
        c=(a+b)//2
        if c**n<=x: a=c
        else: b=c
    return a

def native_tail_norm(r,n,betasq,k,b):
    d=64
    t1=floor(32/log2(b)+.5)
    b1=nearest_even_integer(integer_nth_root(P,t1))
    tmp=sqrt(float(24*n*d))*(betasq/float(r*n*d))
    t2=floor(log2(tmp)/log2(b)+.5)
    b2=nearest_even_integer(pow(tmp,1./t2))
    tri=r*(r+1)//2
    e1=(b1*b1*t1)/12.*float(r*k*d)+(b2*b2*t2)/12.*float(tri*d)
    e2=(b1*b1*t1)/12.*float(tri*d)
    return (2./float(b*b))*(betasq*51.)+e1*e2

def shape(r,n,k,t1,t2):
    # Source enumeration resolves ties at first minimum. For both fixed
    # profiles the exact ratio quotient occurs with integer values.
    tri=r*(r+1)//2
    m=r*k*t1+(t1+t2)*tri
    best=None
    for mu in range(1,m):
        for nu in range(1,2*n):
            v=(abs(m/mu-n/nu),mu,nu)
            if best is None or v<best:
                best=v
            if not v[0]:
                assert Fraction(m,mu)==Fraction(n,nu)
                return 2*nu+mu, max(ceil(n/nu),ceil(m/mu))
    raise ValueError('profile does not have an exact stock ratio match')
