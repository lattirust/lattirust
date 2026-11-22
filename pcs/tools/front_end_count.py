"""Count a fixed circuit for Lattirust Hachi/Greyhound SHAKE256 + sampling.

No wall-clock measurements. Gates for dynamic routing are conservative
bounds for the explicitly described sequential selectors, not compiler output.
"""
from math import ceil, log2

RATE=136

def vec(n,word):return 8+n*word

def frame(tag,n):return 16+len(tag)+n

def ceil_div(a,b):return (a+b-1)//b

class Shake:
    def __init__(self):
        self.length=0
        self.squeezes=0
    def append(self,tag,n):self.length+=frame(tag,n)
    def clone(self,tag,output):
        self.append("challenge",len(tag))
        self.squeezes+=ceil_div(output,RATE)
    def field(self,tag):
        self.clone(tag,96)
        self.append("sample",16)
    def finish(self):return self.length//RATE+self.squeezes

def bounds(p):
    n=p["original_log"]
    d,hr,hc=p["hachi_d"],p["hachi_rows"],p["hachi_cols"]
    gr,gc=p["gh_rows"],p["gh_cols"]
    omega,ghw=p["hachi_weight"],p["gh_weight"]
    ext_domain=p["sumcheck_log"]
    size_key_poly=p["hachi_key_ring_coefficients"]
    ctxt=p["external_context_bytes"]
    init=len("lattirust/hachi/base-field-statement/no-peeling/v1")+24+ctxt
    init+=vec(1,4*d)+vec(n,16)+16+4*d

    t=Shake()
    t.append("protocol",len("lattirust/hachi-shake256/v3"))
    t.append("context",32)
    t.append("phase",len("split-and-fold"))
    for tag in ("ring-degree","gadget-base","folded-digits","sparse-weight"):t.append(tag,8)
    for tag,len_ in (("A",8*hc),("B",8*hr),("D",8*hr)):
        t.append(tag,vec(1,vec(len_,vec(size_key_poly,4))))
    for tag,num in (("outer-commitment",vec(1,vec(d,4))),("ring-value",vec(d,4)),
                    ("left",vec(hc,vec(d,4))),("right",vec(hr,vec(d,4))),
                    ("projected-commitment",vec(1,vec(d,4)))):
        t.append(tag,num)
    for _ in range(hr):
        t.clone("fold-challenge",33*omega)
        t.append("sparse-sample",9*omega)
    t.append("prime-modulus",8);t.append("extension-degree",8)
    t.append("linear-source",len("derived-eq20/v1"))
    t.append("short-table-commitment",vec(1,256))
    t.field("ring-switch")
    for _ in range(ext_domain):t.field("smallness-row")
    for _ in range(3):t.field("relation-row")
    for label,deg in (("smallness",32),("relation",2)):
        t.append("sumcheck",len(label))
        for _ in range(ext_domain):
            t.append("round",vec(deg+1,16))
            t.field("sumcheck-evaluation")
        t.append("endpoint",16)
    t.append("restriction-to-line",vec(ext_domain+1,16))
    t.field("line-point")
    t.clone("continue",32)

    vals=lambda n: vec(n,8+256)
    gh=8+len("lattirust/greyhound/short-relation/v2")+(8+32)+(8+4)+6*8+4+8
    gh+= (8+vals(gc)) + (8+vals(8*gr)) + (8+vals(8*gr))
    gh+= vals(1)+8+1+vals(gr)+vals(gc)+vals(1)+8+vals(1)
    permutations=init//RATE+1+t.finish()+gh//RATE+ceil_div(32+33*gr*ghw,RATE)
    # Conservative count of absorb XORs even for public, redundant, and
    # first-in-block constant wires.
    input_bytes=init+t.length+gh

    # Four candidates per u64 index. Compute comparisons, first-valid mux,
    # long division by public n, then read/write a one-hot-selected array.
    # For use as ring coefficients allocate both sign choices per index.
    def sparse(radix,words,weight):
        w=radix.bit_length()
        lg=w-1
        core_a=4*64+4*64+12+64*2*w+w
        core_x=3*64+12+64*4*w+2*w
        route_a=4*radix+2*radix*lg+radix
        route_x=3*radix*lg+radix
        return {
            "f2_and":words*weight*(core_a+route_a),
            "f2_xor":words*weight*(core_x+route_x),
            "f2_eq":words*weight,
            "lift":words*weight*2*radix
        }
    h=sparse(d,hr,omega);g=sparse(64,gr,ghw)
    out={k:h[k]+g[k] for k in h}
    prime_draws=4*(3*ext_domain+5) # alpha, tau0, tau1, rounds, line
    for k,c in dict(f2_and=351,f2_xor=320,f2_eq=1,lift=32,fp_aff=1).items():
        out[k]=out.get(k,0)+prime_draws*c
    return {
        "initial_Hachi_Greyhound_transcripts":{
            "permutations":permutations,
            "nonconstant_absorb_xor":8*input_bytes,
        },
        "Hachi_Greyhound_field_and_sparse_samplers":out,
        "max_canonical_words":ceil_div(input_bytes,4),
        "absorbed_bytes_including_public":input_bytes,
    }
