import pandas as pd
import sys

def main():
    df_file = "output/moments_output_Europe.tsv"

    max_num_params = 8
    params = ["nu_EUW", "nu_EUE", "T_split", "m_fwd", "m_rev"]
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

    ll_filt = df.ll >= df.ll.max() * 1.000000001
    print(df.ll.max())
    df = df[ll_filt]

    print(df.shape[0])
    print()

    for i in range(len(params)):
        print(params[i])
        foo = df.loc[:, 'est' + str(i)].quantile([0.025, 0.975])
        print(foo)
        print()

if __name__ == "__main__":
    main()
