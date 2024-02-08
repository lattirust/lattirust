from math import log2
from estimator import *
from sage.all import oo # +Infinity

def security_level_l2(n, q, length_bound, m):
    params = SIS.Parameters(n=n, q=q, length_bound=length_bound, m=m, norm=2)
    res = SIS.estimate(params)["lattice"]
    min_cost = min(res["rop"], res["red"])
    return log2(min_cost)

def security_level_linf(n, q, length_bound, m):
    print("linf")
    params = SIS.Parameters(n=n, q=q, length_bound=length_bound, m=m, norm=oo)
    print(params)
    print(SIS.estimate(params))
    res = SIS.estimate(params)["lattice"]
    min_cost = min(res["rop"], res["red"], res["sieve"])
    return log2(min_cost)