"""Native constraint inventory for a *specified* fixed full verifier.

A result is only complete if the independently counted transcript/sampler
stages are provided in manifest extras. Missing stages are printed by name;
they are never silently zero. This is not a single-field R1CS generator.
"""
import argparse
import dataclasses as dc
import json
from collections import Counter
from math import ceil, floor
from pathlib import Path

P=2**32-99

@dc.dataclass
class Count:
    fp_mul: int=0
    fp_aff: int=0
    f2_and: int=0
    f2_xor: int=0
    f2_eq: int=0
    lift: int=0
    ring_mul: int=0    # explanatory subcount only; lowered into fp_mul
    ring_sq: int=0     # explanatory subcount only; lowered into fp_mul
    def __add__(self, rhs):
        return Count(*(x+y for x,y in zip(dc.astuple(self), dc.astuple(rhs))))
    def __mul__(self, n):
        return Count(*(x*n for x in dc.astuple(self)))

class BitCircuit:
    """0/1 are constants; ±wire names encode polarity with no added row."""
    def __init__(self):
        self.row=Count()
        self.w=2
    def input(self):
        out=self.w
        self.w+=1
        return out
    def inputs(self,n):
        return [self.input() for _ in range(n)]
    def inv(self,a):
        return 1-a if a in (0,1) else -a
    def xor(self,a,b):
        if a==0: return b
        if b==0: return a
        if a==1: return self.inv(b)
        if b==1: return self.inv(a)
        if a==b: return 0
        if a==-b: return 1
        self.row.f2_xor+=1
        return self.input()
    def and_(self,a,b):
        if a==0 or b==0 or a==-b: return 0
        if a==1 or a==b: return b
        if b==1: return a
        self.row.f2_and+=1
        return self.input()
    def or_(self,a,b):
        return self.inv(self.and_(self.inv(a),self.inv(b)))
    def mux(self,s,a,b):
        return self.xor(a,self.and_(s,self.xor(a,b)))
    def add(self,a,b,w):
        """Carry-save relation x+y+cin; discard the known out-of-domain carry."""
        r=[]; c=0
        for i in range(w):
            x=a[i] if i<len(a) else 0
            y=b[i] if i<len(b) else 0
            t=self.xor(x,y)
            r.append(self.xor(t,c))
            if i+1<w: c=self.xor(self.and_(x,y),self.and_(t,c))
        return r
    def ge_const(self,x,t):
        if t<0: return 1
        if t>=1<<len(x): return 0
        # Fold lower-bit comparator into successive higher bits.
        ge=1
        for i,y in enumerate(x):
            ge=self.and_(y,ge) if t>>i&1 else self.or_(y,ge)
        return ge
    def square(self,x):
        # Triangular partial products; off-diagonal coefficient doubles are
        # handled as a shift. Evaluate rows through one fixed ripple circuit.
        t=[0]*(2*len(x))
        for i,a in enumerate(x):
            row=[0]*(2*len(x))
            row[2*i]=a
            for j in range(i+1,len(x)):
                row[i+j+1]=self.and_(a,x[j])
            t=self.add(t,row,len(t))
        return t
    def centered_square(self,x):
        assert len(x)==32
        s=self.ge_const(x,(P+1)//2)
        # q - x = q + ~x + 1; only its low 31 bits are needed when s=1.
        c=[(P+1)>>i&1 for i in range(31)]
        opp=self.add([self.inv(y) for y in x[:31]],c,31)
        abs_x=[self.mux(s,a,b) for a,b in zip(x[:31],opp)]
        return self.square(abs_x)

def integer_cutoff(n, bound):
    return min(n*(P//2)**2, floor(bound))

def norm(n, bound):
    """All coefficients already have linked canonical words. No field wrap."""
    circ=BitCircuit()
    single=circ.centered_square(circ.inputs(32))
    row=circ.row*n
    # Accumulator wires are independent symbols after first square. Each
    # old/new width pair repeats the identical gate shape.
    widths=Counter((((i-1)*(P//2)**2).bit_length(), (i*(P//2)**2).bit_length())
                   for i in range(2,n+1))
    for (prev,after), mult in widths.items():
        c=BitCircuit()
        c.add(c.inputs(prev),c.inputs(len(single)),after)
        row+=c.row*mult
    cutoff=integer_cutoff(n,bound)
    if cutoff<(n*(P//2)**2):
        circ=BitCircuit()
        circ.ge_const(circ.inputs((n*(P//2)**2).bit_length()),cutoff+1)
        row+=circ.row+Count(f2_eq=1)
    return row

def canonical(n):
    """Unique p input words; range, mixed equality, and reconstruction."""
    c=BitCircuit()
    c.ge_const(c.inputs(32), P)
    return (c.row+Count(fp_aff=1, f2_eq=1, lift=32))*n

def perms(n, absorb_xor=0):
    return Count(f2_and=38400*n,f2_xor=115200*n+absorb_xor)

def r_ops(d,m,s):
    # Bilinear products ai*bj are individual outputs; their signed cyclic
    # convolution sums remain affine until explicitly materialized.
    return Count(fp_mul=d*d*m + d*(d+1)//2*s,ring_mul=m,ring_sq=s)

def e_ops(n):
    # Two-level Karatsuba, output real/imag affine. Nine base products.
    return Count(fp_mul=9*n)

def hachi(p):
    d=p["hachi_d"]; rows=p["hachi_rows"]; cols=p["hachi_cols"]
    width=p["hachi_weight"]; zeta=p["folded_digits"]; delta=p["delta"]
    assert p["hachi_base"]==16 and p["original_point_private"]
    assert p["hachi_ranks"]==[1,1,1]
    assert d//4==1<<(d.bit_length()-3)
    gh_rows,gh_cols=p["gh_rows"],p["gh_cols"]
    logd=d.bit_length()-1
    active=rows*delta*2+cols*delta*zeta+5*delta
    logw=(active-1).bit_length()
    n=logd+logw
    assert 1<<n == gh_rows*gh_cols*16
    # See count.py; apply the same partially factored endpoint method to
    # any dyadic row/column geometry. No dynamic hash bytes removed.
    sparse=rows*width
    e=(d-1) + 6 + (2**logw-2) + (2*rows+cols)*delta
    e+= 3*rows+cols + 5+6 +(2*logd-1) + 2
    e+= sparse*(logd-1) + 3 + n*(32+2)
    e+= (2*n-1) +16+1+2*n
    e+= (max(rows-2,0)+max(cols-2,0)+(d//4-2))+(d//4)
    # E weights for full separated table, including 4 within-word axes.
    e+=(max(gh_rows-2,0)+16*gh_cols-2)
    b=2*(d-1)*4+(rows+cols)*3*4 + sparse*(logd+1)*4
    if p.get("original_commitment_private"): b+=(d-1)*4
    return dict(
       Hachi_Fig7_and_trace_arithmetic=e_ops(e)+Count(fp_mul=b, fp_aff=4*(2*n+4)+8),
       Hachi_public_CRS_Fp_materializations=Count(fp_aff=4*(2*rows+cols)*delta),
    )

def lab(p):
    d=p['lab_d'];r=p['lab_r'];n=p['lab_n'];g=p['lab_aggregates']
    k=p['lab_k'];k1=p['lab_k1'];k2=p['lab_k2']
    rt=p['tail_r'];nt=p['tail_n']
    gd=p['gh_depth'];gr=p['gh_rows'];gc=p['gh_cols']
    assert d==64 and p['gh_ranks']==[1,1,1] and p['gh_input_depth']==1
    assert p['lab_full_rows']==5 and p['lab_zero_ct_rows']==1
    assert r*n>=gc+2*gr*gd
    tri=r*(r+1)//2
    # Three dynamically generated row supports from GH (others are fixed).
    dyn_gh=3*gr*gd+gc
    phase1=g*256*d*r*n
    result={
        'LAB_JL_ternary_bit_interface':Count(lift=2*256*d*r*n),
        'LAB_projection_RHS':Count(fp_mul=g*256,fp_aff=g),
        'LAB_projection_norm':norm(256,128*float(p['lab_input_norm_sq'])),
        'LAB_projection_coefficients':Count(fp_mul=phase1),
        'LAB_row_aggregation_and_fold':r_ops(d,g*r*n+dyn_gh+r*n+7+tri-r,r),
        'LAB_terminal_syndromes_and_Gram':r_ops(d,k*r+2*tri+n,n)
          +Count(fp_aff=(k+3+k1+k2)*d),
        'LAB_terminal_norm':norm(d*rt*nt,float.fromhex(p['lab_tail_norm_sq_hex'])),
    }
    return result


def lab_fp_samples(p):
    n=p['lab_aggregates']*(256+p['lab_zero_ct_rows'])+p['lab_d']*(p['lab_full_rows']+p['lab_aggregates'])
    # 6 independent words + lt-p checks; mux in reverse to prioritize first.
    c=BitCircuit()
    candidates=[c.inputs(32) for _ in range(6)]
    accepted=[c.inv(c.ge_const(a,P)) for a in candidates]
    ok=accepted[-1]
    selected=candidates[-1]
    for flag,candidate in zip(reversed(accepted[:-1]),reversed(candidates[:-1])):
        selected=[c.mux(flag,s,w) for s,w in zip(selected,candidate)]
        ok=c.or_(ok,flag)
    return (c.row+Count(fp_aff=1,f2_eq=1,lift=32))*n

def lab_transcript(p):
    k1=p['lab_k1'];k2=p['lab_k2']
    g=p['lab_aggregates'];K=p['lab_full_rows'];J=p['lab_zero_ct_rows']
    d=p['lab_d'];r=p['lab_r'];n=p['lab_n']
    A=[16+d*4*k1, 16+256*4, 8+g*d*4, 16+d*4*k2]
    S=[256*n*r*d,24*g*(256+J),1536*(K+g),65280*r]
    duplex=sum((a-1)//136 + (s+135)//136 for a,s in zip(A,S))
    op=[
        ('A',A[0],'u1'),('S',S[0],'Pi'),('A',A[1],'p'),
        ('S',24*g*J,'psi'),('S',24*g*256,'omega'),('A',A[2],'b'),
        ('S',1536*K,'alpha'),('S',1536*g,'beta'),('A',A[3],'u2'),
        ('S',S[3],'c'),
    ]
    tag_length=len('lattirust/greyhound/labrador/round/v1/')+64+sum(1+len(t+str(v)+l) for t,v,l in op)
    # Nimue uses overwrite; chain SHA3 has one block from zero state.
    tag_count=(tag_length+135)//136
    chain=(len('lattirust/greyhound/labrador/explicit-rounds/v1')+64)//136+1
    assert chain==1  # no variable absorb XOR in first block (zero state)
    c=BitCircuit()
    digit=c.inputs(4)
    numerical=c.add(digit,[0,0,0,0,1,1],8) # nibble + '0'
    alpha=c.add(digit,[1,1,1,0,1,0,1],8)  # nibble + ('a' - 10)
    s=c.ge_const(digit,10)
    for a,b in zip(numerical,alpha): c.mux(s,a,b)
    return perms(duplex+tag_count+chain)+c.row*64, (A,S,tag_length,duplex,tag_count,chain)

EXPECTED_EXTRAS={
    'initial_Hachi_Greyhound_transcripts','LaBRADOR_transcript_io',
    'Hachi_Greyhound_field_and_sparse_samplers','LaBRADOR_field_samplers',
    'LaBRADOR_fold_sampler','canonical_field_io','other_affine_materializations',
}

def run(p):
    stages=hachi(p)|lab(p)
    stages['LaBRADOR_transcript_io'], trace=lab_transcript(p)
    stages['LaBRADOR_field_samplers']=lab_fp_samples(p)
    from exact_fold_sampler import count, read_roots
    roots=(read_roots(p['fold_sampler_rust']) if p.get('fold_sampler_rust') else
           json.loads((Path(__file__).parent/'fold_cosines.json').read_text()))
    # from exact_fold_sampler imports a second copy if run as __main__.
    stage=count(roots)['per_challenge']
    stages['LaBRADOR_fold_sampler']=Count(**dc.asdict(stage))*p['lab_r']
    if p.get('include_frontend'):
        from front_end_count import bounds
        front=bounds(p)
        stages['initial_Hachi_Greyhound_transcripts']=perms(**{
            'n':front['initial_Hachi_Greyhound_transcripts']['permutations'],
            'absorb_xor':front['initial_Hachi_Greyhound_transcripts']['nonconstant_absorb_xor'],
        })
        stages['Hachi_Greyhound_field_and_sparse_samplers']=Count(**front['Hachi_Greyhound_field_and_sparse_samplers'])
        # Reserve one generic bridge slot per four bytes of the entire
        # front transcript to upper-bound the real field word occurrences.
        # Slots that are not field words are dummy; never constrain digest
        # bytes to be field canonical.
        g=p['lab_aggregates']
        d=p['lab_d']
        lab_words=d*(p['lab_k1']+p['lab_k2']+g+p['tail_r']*p['tail_n'])+256
        stages['canonical_field_io']=canonical(front['max_canonical_words']+lab_words)
        # All 16-byte matrix shape headers + g vector header +
        # terminal outer vector count. Other structural arities fixed by IR.
        stages['other_affine_materializations']=Count(f2_eq=8*(56+8+16*p['tail_r']))
    for tag, counters in p.get('transcripts',{}).items():
        stages[tag]=perms(counters['permutations'],counters['nonconstant_absorb_xor'])
    for tag, counters in p.get('extras',{}).items():
        stages[tag]=Count(**counters)
    slots=p.get('distinct_private_canonical_field_words')
    if slots is not None: stages['canonical_field_io']=canonical(slots)
    missing=EXPECTED_EXTRAS-stages.keys()
    total=sum(stages.values(),Count())
    print('profile',p['name'])
    print('LAB absorb, squeeze, IOPattern tag bytes, duplex/tag/chain permutations',trace)
    for tag,count in stages.items(): print(tag, ' '.join(f'{k}:{v}' for k,v in dc.asdict(count).items() if v))
    print('MISSING stages' if missing else 'All defined stages populated',', '.join(sorted(missing)))
    print('INCOMPLETE subtotal' if missing else 'TOTAL in this native mixed gate model',dc.asdict(total))
    print('fp_aff counts explicit equality/copy rows. Fp fixed linear forms are inlined without extra rows except where a stage counts materialization.')
    print('Ring mul/square are included in fp_mul, not extra. F2 XOR rows separate. Cross-field lift is its own gadget; do not add into a single-field R1CS total.')

if __name__=='__main__':
    a=argparse.ArgumentParser()
    a.add_argument('manifest')
    args=a.parse_args()
    with open(args.manifest) as f: run(json.load(f))
