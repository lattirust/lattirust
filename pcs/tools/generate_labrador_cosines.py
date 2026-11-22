"""Generate the fixed-point cosine table used by the LattiRust verifier.

All operations use integer interval arithmetic, including the half-angle
construction of pi/64 and subsequent powers of the complex root. Every
printed integer is within 1 of 2^20 * cos(i*pi/64).
"""
from math import isqrt

S=1<<20
T=1<<80
def const(n):return n*T,n*T
def plus(a,b):return a[0]+b[0],a[1]+b[1]
def minus(a,b):return a[0]-b[1],a[1]-b[0]
def times(a,b):
    outs=[x*y for x in a for y in b]
    return min(outs)//T,-((-max(outs))//T)
def sqrt(a):
    lo,hi=a
    return isqrt(lo*T),isqrt(hi*T)+1
def half(a):return a[0]//2,-((-a[1])//2)
cos=const(0)
for _ in range(5):cos=sqrt(half(plus(const(1),cos)))
sin=sqrt(minus(const(1),times(cos,cos)))
re,im=const(1),const(0)
out=[]
for i in range(128):
    lo,hi=re
    w=((lo+hi)*S+T)//(2*T)
    assert (w-1)*T<lo*S and hi*S<(w+1)*T
    out.append(w)
    re,im=minus(times(re,cos),times(im,sin)),plus(times(re,sin),times(im,cos))
for i in range(0,128,8):print("    "+", ".join(map(str,out[i:i+8]))+",")
