import moments
import numpy as np

"""
Utilities file that stores functions that construct and return SFSs from models
defined with moments' demographic modeling language.
"""

# Returns square matrix with dim n_pops that equals m everywhere, except zeros on diag
def get_mig_mat(n_pops: int, m: float) -> np.ndarray:
    """
    :param n_pops int: Number of populations
    :param m float: Migration rate, in coalescent units of `1 / 2N_anc`
    :return: Square `n_pops` x `n_pops` matrix representing symmetric 
        migration, equals `m` everywhere except for zeros along the diagonal
    :rtype: np.ndarray
    """
    return m * np.ones([n_pops] * 2) - np.diag([m] * n_pops)

# An ancestral population instantaneously changes in size at a specified time in
# the past. In the paper, we also refer to this model as "one-population." This 
# function is a replicate of https://moments.readthedocs.io/en/main/api/api_moments.html#moments.Demographics1D.two_epoch.
def two_epoch(params: list, ns: list, pop_ids: list=None) -> moments.Spectrum:
    """
    :param params list: List of float parameters for the model
    :param ns list: List of int haploid sample sizes for each population
    :param pop_ids list: List of string population IDs
    :return: Folded SFS from a two-epoch model
    """
    nu, T = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts, pop_ids=pop_ids)
    
    fs.integrate([nu], T)
    return fs

# Ancestral pop. splits into two pops., with symmetric migration post-split
def split(params, ns, pop_ids=None):
    nu1, nu2, T, m = params
    mig_mat = get_mig_mat(2, m)

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)

    fs = fs.split(0, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=mig_mat)

    fs.pop_ids = pop_ids
    return fs

# Ancestral pop. splits into two pops., with asymmetric migration post-split
def split_asymmig(params, ns, pop_ids=None):
    nu1, nu2, T, m_2_to_1, m_1_to_2 = params
    mig_mat = np.array([[0, m_2_to_1], [m_1_to_2, 0]])

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1])
    fs = moments.Spectrum(sts)

    fs = fs.split(0, ns[0], ns[1])
    fs.integrate([nu1, nu2], T, m=mig_mat)

    fs.pop_ids = pop_ids
    return fs

# Ancestral pop. splits into two pops., then they admix into a third pop., with
# symmetric migration at different strengths at the two-pop. and three-pop. 
# stages. The admixed population will always be the third population, i.e. idx=2.
def admixture(params, ns, pop_ids=None):
    nu1, nu2, nu_admix, T_split, T_admix, m2, m3, admix_prop = params
    mig_mat2 = get_mig_mat(2, m2)
    mig_mat3 = get_mig_mat(3, m3)

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + 2 * ns[2])
    fs = moments.Spectrum(sts)

    fs = fs.split(0, ns[0] + ns[2], ns[1] + ns[2])
    fs.integrate([nu1, nu2], T_split, m=mig_mat2)

    fs = fs.admix(0, 1, ns[2], admix_prop)
    fs.integrate([nu1, nu2, nu_admix], T_admix, m=mig_mat3)

    fs.pop_ids = pop_ids
    return fs.fold()

# Ancestral pop. splits into two pops., pop #1 and an intermediate pop., then 
# the intermediate pop. splits into pops. #2 and #3, with symmetric migration 
# at different strengths at the two-pop. and three-pop. stages
def twosplits(params, ns, pop_ids=None):
    nu1, nu_intermediate, nu2, nu3, T1, T2, m2, m3 = params

    mig_mat2 = get_mig_mat(2, m2)
    mig_mat3 = get_mig_mat(3, m3)
    
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)

    fs = fs.split(0, ns[0], ns[1] + ns[2])
    fs.integrate([nu1, nu_intermediate], T1, m=mig_mat2)

    fs = fs.split(1, ns[1], ns[2])
    fs.integrate([nu1, nu2, nu3], T2, m=mig_mat3)

    fs.pop_ids = pop_ids
    return fs
