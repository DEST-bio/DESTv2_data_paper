import msprime

def get_one_population_dem(n_anc, params):
    nu, T = params
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=n_anc * nu)
    dem.add_population_parameters_change(time=2 * n_anc * T, population="A", initial_size=n_anc)
    return dem

def get_split_dem(n_anc, params):
    nu1, nu2, T, m = params
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=n_anc*nu1)  # pop1 at present time
    dem.add_population(name="B", initial_size=n_anc*nu2)  # pop2 at present time
    dem.add_population(name="ancestral", initial_size=n_anc)  # ancestral pop
    dem.add_population_split(time=2*n_anc*T, derived=["A", "B"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["A", "B"], m/(2*n_anc))
    return dem

def get_split_asymmig_dem(n_anc, params):
    nu1, nu2, T, m_A_to_B, m_B_to_A = params
    dem = msprime.Demography()
    dem.add_population(name="A", initial_size=n_anc*nu1)  # pop1 at present time
    dem.add_population(name="B", initial_size=n_anc*nu2)  # pop2 at present time
    dem.add_population(name="C", initial_size=n_anc)  # ancestral pop
    dem.add_population_split(time=2*n_anc*T, derived=["A", "B"], ancestral="C")
    # These migration rates are backwards-in-time, so they should be swapped when
    # interpreted in natural units!
    dem.set_migration_rate(source="A", dest="B", rate=m_A_to_B/(2*n_anc))
    dem.set_migration_rate(source="B", dest="A", rate=m_B_to_A/(2*n_anc))
    return dem

def get_admixture_dem(n_anc, params):
    # It is necessary that T_split > T_admix because the split from the ancestral 
    # precedes the admixture event.
    nu1, nu2, nu_admix, T_split, T_admix, m2, m3, p_admix = params

    # Initialize populations in `Demography` object
    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=n_anc)
    dem.add_population(name="A", initial_size=n_anc*nu1)
    dem.add_population(name="B", initial_size=n_anc*nu2)
    dem.add_population(name="admixed", initial_size=n_anc*nu_admix)

    # Admixture of two derived populations
    dem.add_admixture(derived="admixed", ancestral=["A", "B"],
    				  time=2*n_anc*T_admix, proportions=[p_admix, 1 - p_admix])
    dem.set_symmetric_migration_rate(["A", "B", "admixed"], m3/(2*n_anc))

    # Split of two derived populations from ancestral population
    dem.add_population_split(time=2*n_anc*T_split, derived=["A", "B"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["A", "B"], m2/(2*n_anc))

    return dem

def get_two_splits_dem(n_anc, params):
    # `T2 < T1` must hold
    nu1, nu_intermediate, nu2, nu3, T1, T2, m2, m3 = params

    dem = msprime.Demography()
    dem.add_population(name="ancestral", initial_size=n_anc)
    dem.add_population(name="A", initial_size=n_anc*nu1)
    dem.add_population(name="intermediate", initial_size=n_anc*nu_intermediate)
    dem.add_population(name="B", initial_size=n_anc*nu2)
    dem.add_population(name="C", initial_size=n_anc*nu3)

    dem.add_population_split(time=2*n_anc*T2, derived=["B", "C"], ancestral="intermediate")
    dem.set_symmetric_migration_rate(["A", "B", "C"], m3/(2*n_anc))
    dem.add_population_split(time=2*n_anc*T1, derived=["A", "intermediate"], ancestral="ancestral")
    dem.set_symmetric_migration_rate(["A", "B"], m2/(2*n_anc))
    return dem