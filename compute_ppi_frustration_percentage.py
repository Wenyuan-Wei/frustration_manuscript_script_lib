import json
import pandas as pd

def get_frustration(file_path, sys, first_chain):
    df_raw = pd.read_csv(file_path, delimiter=' ')
    #Drop other cols
    df_drop_cols = df_raw[['Res1', 'Res2', 'ChainRes1', 'ChainRes2', 'AA1', 'AA2', 'FrstState']]
    # Switch values if ChainRes2 is e3 chain. Make sure e3 chain alway first
    mask = df_drop_cols['ChainRes2'] == first_chain
    df_drop_cols.loc[mask, ['Res1', 'Res2']] = df_drop_cols.loc[mask, ['Res2', 'Res1']].values
    df_drop_cols.loc[mask, ['ChainRes1', 'ChainRes2']] = df_drop_cols.loc[mask, ['ChainRes2', 'ChainRes1']].values
    df_drop_cols.loc[mask, ['AA1', 'AA2']] = df_drop_cols.loc[mask, ['AA2', 'AA1']].values

    df1 = df_drop_cols[df_drop_cols['ChainRes1'] != df_drop_cols['ChainRes2']]
    df1.insert(0, 'sys', sys)
    return df1

def extract_ppi_pair(diction_in, general_dir, first_chain_name):
    df_list = []
    for sys, val in diction_in.items():
        first_chain = val[first_chain_name]
        #ga_chain = val['G_alpha_chain']
        frust_path = f"{general_dir}/frust/frus_output/{sys.lower()}.done/FrustrationData/{sys.lower()}.pdb_mutational"
        
        df_frust = get_frustration(frust_path,sys, first_chain)
        df_list.append(df_frust)
    all_sys_df = pd.concat(df_list, axis=0, ignore_index=True)
    return all_sys_df

def frustration_level_one_sys(df_in, frststate):
    total_count = df_in['FrstState'].count()
    # Compute the frequency of the value "frststate" in the "FrstState" column
    #frststate could be "highly" or "minimally"
    state_count = df_in['FrstState'].value_counts().get(frststate, 0)
    state_fre = state_count/total_count              
    return state_fre

def main(df_raw, dict_chain, frststate):
    sys_list = list(df_raw["sys"].unique())
    #degrader is the PPI with degrader on it. Other are natural PPI without degrader
    degrader_frust_list = []
    other_frust_list= []
    degrader_name_list = []
    other_name_list= []
    for sys in sys_list:
        #Get e3,poi chain, if raw data is gpcr_ga, change it to "gpcr", "ga"
        e3_chain, poi_chain = dict_chain[sys]['E3_chain'], dict_chain[sys]['POI_chain']
        df_sys = df_raw[df_raw["sys"] == sys]
        #Get array of ppi
        ppi_array = df_sys.apply(lambda row: 
                            tuple(sorted([row["ChainRes1"], row["ChainRes2"]])), axis=1).unique()
        for ppi in ppi_array:
            df_ppi = df_sys[df_sys.apply(lambda row: 
                            set([row["ChainRes1"], row["ChainRes2"]]) == set(ppi), axis=1)]
            ppi_percentage = frustration_level_one_sys(df_ppi, frststate)
            if set(ppi) == set([e3_chain, poi_chain]):
                degrader_frust_list.append(ppi_percentage)
                degrader_name_list.append(sys)
            else:
                other_frust_list.append(ppi_percentage)
                other_name_list.append(sys)
    df_degrader = pd.DataFrame({
    'structure': degrader_name_list,
    'frust_percentage': degrader_frust_list
    })
    df_other = pd.DataFrame({
    'structure': other_name_list,
    'frust_percentage': other_frust_list
    })
    df_degrader.to_csv(f"degrader_PPI_{frststate}.csv", index=None)
    df_other.to_csv(f"other_PPI_{frststate}.csv", index=None)

if __name__ == "__main__":
    #ldd_complex.json is the json file recording chain ID for E3 and POI chain. Contents looks like:
    #"5T35": {
        #     "name": "vhl-brd4",
        #     "E3_chain": "D",
        #     "POI_chain": "A"
        # },
    with open("ldd_complex.json", "r") as file:
        dict_chain = json.load(file)
    
    #Extract inter-protein pairs
    general_dir = "../../whole_complex"
    df_ppi = extract_ppi_pair(diction_in=dict_chain, general_dir=general_dir, first_chain_name="E3_chain")
    
    main(df_ppi, dict_chain, "highly")
    main(df_ppi, dict_chain, "minimally")


