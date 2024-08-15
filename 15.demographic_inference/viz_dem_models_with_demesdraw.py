import msprime
import demesdraw
import matplotlib.pyplot as plt
import msprime_models as msm

def main():
    #dems = {"one-population": msm.get_one_population_dem(1e4, [2, 1]),
    #        "split": msm.get_split_dem(1e4, [1, 1, 1, 1]),
    #        "admixture": msm.get_admixture_dem(1e4, [1, 1, 1, 2, 1, 0, 0, 0.5]),
    #        "two-splits": msm.get_two_splits_dem(1e4, [1, 1, 1, 1, 2, 1, 1, 1])}
    # 0.0010982661644379046 / (2 * n_anc), 0.0686663736024 / (2 * n_anc), 
    n_anc = 1773#44
    dems = {"admixture": msm.get_admixture_dem(n_anc, 
                                               [1.6484174485803664, 
                                                0.2819202780787078, 
                                                1.2418220349493856, 
                                                0.06569726302718709, 
                                                0.00955907372305152, 
                                                0.01,
                                                0.01,
                                                0.5395910636154568])}
    output_dir = "draft_figures/"

    viz_dem_models_with_demesdraw(dems, output_dir)

def viz_dem_models_with_demesdraw(dems, output_dir):
    for name, dem in dems.items():
        print(name)
        graph = msprime.Demography.to_demes(dem)
        fig, ax = plt.subplots()
        demesdraw.tubes(graph, ax=ax, seed=1)
        plt.savefig(f"{output_dir}{name}_demesdraw.png")

if __name__ == "__main__":
    main()
