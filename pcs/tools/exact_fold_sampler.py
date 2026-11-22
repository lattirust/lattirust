"""Boolean inventory for the exact Rust short challenge sampler.

Same random bytes, one fixed circuit that evaluates all 128 candidates.
Does not count Keccak generation: 65,280 bytes per returned R64 challenge.
"""
import re, sys
from pathlib import Path
from full_count import Count, BitCircuit

N=64
TRIALS=128

def bits(x,n):
    return [(x>>i)&1 for i in range(n)]

def upset(c, xs):
    out=0
    for x in xs: out=c.or_(out,x)
    return out

def product(c, xs):
    out=1
    for x in xs: out=c.and_(out,x)
    return out

def equal(c, x, u):
    return product(c,(t if b else c.inv(t) for t,b in zip(x,bits(u,len(x)))))

def rem_const(c, word, modulus):
    k=(modulus-1).bit_length()
    if modulus==1<<k:
        return word[:k]
    m=[0]*k
    for x in reversed(word):
        # t=2*m+x in [0,2*modulus). One subtraction suffices.
        t=[x]+m
        ge=c.ge_const(t,modulus)
        sub=c.add(t,bits(-modulus,k),k)
        m=[c.mux(ge,a,b) for a,b in zip(t,sub)]
    return m

def candidate(c, roots):
    flags=[]
    perm=[bits(i,6) for i in range(N)]
    for i in range(N-1,0,-1):
        n=i+1
        word=c.inputs(64)
        flags.append(c.ge_const(word,2**64%n))
        j=rem_const(c,word,n)
        selectors=[equal(c,j,k) for k in range(n)]
        prev=perm[i]
        selected=[0]*6
        for k in range(n):
            selected=[c.xor(a,c.and_(b,selectors[k])) for a,b in zip(selected,perm[k])]
        for k in range(i):
            perm[k]=[c.mux(selectors[k],a,b) for a,b in zip(perm[k],prev)]
        perm[i]=selected
    # FIFO is shifted once for every nonzero coefficient; indices not leaked.
    fifo=c.inputs(41)
    out=[]
    for i in range(N):
        nz=c.ge_const(perm[i],23)
        twice=c.ge_const(perm[i],54)
        sign=fifo[0]
        out.append((nz,twice,sign))
        if i+1<N:
            fifo=[c.mux(nz,x,y) for x,y in zip(fifo,fifo[1:]+[0])]
    sampler=c.row
    c.row=Count()

    def dot(factors):
        out_= [0]*27
        for (nz,twice,sign), a in zip(out, factors):
            pos=[c.mux(twice,x,y) for x,y in zip(bits(a,27),bits(2*a,27))]
            neg=[c.mux(twice,x,y) for x,y in zip(bits(-a,27),bits(-2*a,27))]
            raw=[c.and_(nz,c.mux(sign,x,y)) for x,y in zip(pos,neg)]
            out_=c.add(out_,raw,27)
        return out_

    def radius(x):
        sign=x[-1]
        neg=c.add([c.inv(a) for a in x[:26]], [1], 26)
        abs_x=[c.mux(sign,a,b) for a,b in zip(x,neg)]
        return c.square(c.add(abs_x,bits(51,26),26))

    for j in range(32):
        args=[((2*j+1)*i)&127 for i in range(N)]
        re=dot([roots[a] for a in args])
        im=dot([roots[(160-a)&127] for a in args])
        rsq,ssq=radius(re),radius(im)
        norm=c.add(rsq,ssq,53)
        flags.append(c.inv(c.ge_const(norm,(15*2**20)**2)))
    good=product(c,flags)
    spectral=c.row
    return good,out,sampler,spectral

def count(roots):
    c=BitCircuit();_,_,parser,spectral=candidate(c,roots)
    c=BitCircuit()
    active=None
    good=0
    for _ in range(TRIALS):
        ok=c.input()
        out=c.inputs(3*N)
        good=c.or_(good,ok)
        active=out if active is None else [c.mux(ok,a,b) for a,b in zip(active,out)]
    c.row.f2_eq+=1
    final=c.row+Count(lift=3*N,fp_mul=N)
    total=(parser+spectral)*TRIALS+final
    return dict(parser=parser,spectral=spectral,selection_lifting=final,per_challenge=total)

def read_roots(path):
    text=Path(path).read_text()
    coeffs=list(map(int,re.findall(r'-?\d+',text.split('const COS: [i64; 128] = [',1)[1].split(']',1)[0])))
    assert len(coeffs)==128
    return coeffs

if __name__=='__main__':
    for k,v in count(read_roots(sys.argv[1])).items():
        print(k,v)

class Evaluate(BitCircuit):
    """Evaluate the same circuit on witness values, for semantic controls."""
    def __init__(self,xs):
        super().__init__()
        self.x=iter(xs)
        self.values={}
    def inputs(self,n):
        a=super().inputs(n)
        for q in a:self.values[q]=next(self.x)
        return a
    def value(self,a):
        if a in (0,1):return a
        return (a<0)^self.values[abs(a)]
    def xor(self,a,b):
        t=self.w
        c=super().xor(a,b)
        x=self.value(a)^self.value(b)
        if t!=self.w:self.values[c]=x
        else:assert self.value(c)==x
        return c
    def and_(self,a,b):
        t=self.w
        c=super().and_(a,b)
        x=self.value(a)&self.value(b)
        if t!=self.w:self.values[c]=x
        else:assert self.value(c)==x
        return c
