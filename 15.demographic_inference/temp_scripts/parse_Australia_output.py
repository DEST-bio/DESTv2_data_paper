import pandas as pd
import sys

def main():
    df_file = "output/moments_output_Australia.tsv"

    max_num_params = 8
    params = ["nu1", "nu2", "nu_admix", "T_split", "T_admix", "m2", "m3", "admix_prop"]
    names = ['model', 'pop_of_interest'] + \
                           ['init' + str(i) for i in range(max_num_params)] + \
                           ['est' + str(i) for i in range(max_num_params)] + \
                           ['upper_bound' + str(i) for i in range(max_num_params)] + \
                           ['ll', 'coll_pop_ll', 'func_calls', 'grad_calls',
                            'maxiter', 'hour_limit']
    names += ['dummy']
    names += ['region']

    df = pd.read_csv(df_file, sep='\s+', header=None, names=names)

    print(df.pop_of_interest.unique())

    #pop_of_interest = int(sys.argv[1])
    for pop_of_interest in range(4):
        print(pop_of_interest, end="\t")

        df_pop = df[df.pop_of_interest == pop_of_interest]

        print(df_pop.shape[0], end="\t")

        ll_filt = df_pop.ll >= df_pop.ll.max() * 1.0000000001
        df_pop = df_pop[ll_filt]

        print(df_pop.shape[0], end="\t")

        print(df_pop.ll.max())

        for i in range(len(params)):
            #print(params[i])
            qs = df_pop.loc[:, 'est' + str(i)].quantile([0.025, 0.975])
            print(f"{qs[0.025]}\t{qs[0.975]}")

if __name__ == "__main__":
    main()
