from estimator import *
from math import log2
from sage.all import oo  # +Infinity

Logging.set_level(Logging.CRITICAL)


def sis_security_level_l2(n, q, length_bound, m):
    params = SIS.Parameters(n=n, q=q, length_bound=length_bound, m=m, norm=2)
    res = SIS.estimate(params)["lattice"]
    min_cost = min(res["rop"], res["red"])
    return log2(min_cost)


def sis_security_level_linf(n, q, length_bound, m):
    params = SIS.Parameters(n=n, q=q, length_bound=length_bound, m=m, norm=oo)
    res = SIS.estimate(params)["lattice"]
    min_cost = min(res["rop"], res["red"], res["sieve"])
    return log2(min_cost)
