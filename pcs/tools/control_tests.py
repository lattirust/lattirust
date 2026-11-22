from random import Random
from full_count import P
from exact_fold_sampler import Evaluate, bits, rem_const, read_roots, candidate

rng=Random(719)

for n in range(2,65):
    for _ in range(4):
        raw=rng.getrandbits(64)
        C=Evaluate(bits(raw,64))
        out=rem_const(C,C.inputs(64),n)
        assert sum(C.value(t)<<j for j,t in enumerate(out))==raw%n

for x in [0,1,P//2,(P+1)//2,P-1]+[rng.randrange(P) for _ in range(12)]:
    C=Evaluate(bits(x,32))
    out=C.centered_square(C.inputs(32))
    assert sum(C.value(t)<<j for j,t in enumerate(out))==min(x,P-x)**2

import json
roots=json.load(open('fold_cosines.json'))
assert len(roots)==128
for _ in range(12):
    idxs=[rng.getrandbits(64) for _ in range(63)]
    signs=bits(rng.getrandbits(41),41)
    C=Evaluate([key for x in idxs for key in bits(x,64)]+signs)
    ok,w,_,_=candidate(C,roots)
    perm=list(range(64)); passes=True
    for i,x in zip(range(63,0,-1),idxs):
        passes &= x>=2**64%(i+1)
        j=x%(i+1)
        perm[i],perm[j]=perm[j],perm[i]
    expected=[]
    next_sign=0
    for idx in perm:
        if idx<23: expected.append(0); continue
        expected.append((1+int(idx>=54)) * (1-2*signs[next_sign]))
        next_sign+=1
    recovered=[(C.value(nz)+C.value(two))*(1-2*C.value(sgn)) for nz,two,sgn in w]
    assert recovered==expected
    for j in range(32):
        real=sum(x*roots[((2*j+1)*i)&127] for i,x in enumerate(expected))
        imag=sum(x*roots[(160-(2*j+1)*i)&127] for i,x in enumerate(expected))
        passes &= (abs(real)+51)**2+(abs(imag)+51)**2<(15*2**20)**2
    assert C.value(ok)==passes
print('PASS: 252 modular reductions, 17 centered squares, and 12 exact independent spectral controls')
