import pandas as pd
import sys

def main():
    df_file = "output/moments_output_Transatlantic_expandedAfr.tsv"

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

    print(df.shape[0])

    pop_of_interest = int(sys.argv[1])
    df = df[df.pop_of_interest == pop_of_interest]

    print(df.shape[0])

    ll_filt = df.ll >= df.ll.max() * 1.00001
#    print(df.ll.nlargest(2).iloc[1])
#    ll_filt = df.ll == df.ll.nlargest(2).iloc[1]
    df = df[ll_filt]

    print(df.shape[0])
    print()

    for i in range(max_num_params):
        print(params[i])
        foo = df.loc[:, 'est' + str(i)].quantile([0.025, 0.975])
        print(foo)
        print()

if __name__ == "__main__":
    main()
