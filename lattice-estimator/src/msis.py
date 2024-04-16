from MSIS_security import *


def msis_security_level_l2(n, d, q, length_bound, m):
    params = MSISParameterSet(n=d, h=n, w=m, q=q, B=length_bound, norm="l2")
    (_, c_pc, c_pq, c_pp) = MSIS_summarize_attacks(params)
    min_cost = min(c_pc, c_pq, c_pp)
    return min_cost


def msis_security_level_linf(n, d, q, length_bound, m):
    params = MSISParameterSet(n=d, h=n, w=m, q=q, B=length_bound, norm="linf")
    (_, c_pc, c_pq, c_pp) = MSIS_summarize_attacks(params)
    min_cost = min(c_pc, c_pq, c_pp)
    return min_cost
