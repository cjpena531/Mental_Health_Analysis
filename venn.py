import pandas as pd
from matplotlib_venn import venn3

b = pd.read_csv("p_val_ancg_b.csv")
m = pd.read_csv("p_val_ancg_m.csv")
s = pd.read_csv("p_val_ancg_s.csv")

s = s.rename(columns={"Unnamed: 0": "join"})
m = m.rename(columns={"Unnamed: 0": "join"})
b = b.rename(columns={"Unnamed: 0": "join"})

joined = b.set_index('join').join(m.set_index('join'),how="outer", lsuffix='_b', rsuffix='_m')
joined = s.set_index('join').join(joined, how="outer")
#joined = b.join(m, on="Unnamed:0", how="outer", lsuffix='_b', rsuffix='_m')
joined = joined[['baseMean', 'baseMean_b', 'baseMean_m']]

just_s = joined[(joined['baseMean'].isnull() == False) & 
                (joined['baseMean_b'].isnull() == True) & 
                (joined['baseMean_m'].isnull() == True)]
just_b = joined[(joined['baseMean'].isnull() == True) & 
                (joined['baseMean_b'].isnull() == False) & 
                (joined['baseMean_m'].isnull() == True)]
just_m = joined[(joined['baseMean'].isnull() == True) & 
                (joined['baseMean_b'].isnull() == True) & 
                (joined['baseMean_m'].isnull() == False)]

sm = joined[(joined['baseMean'].isnull() == False) & 
                (joined['baseMean_b'].isnull() == True) & 
                (joined['baseMean_m'].isnull() == False)]
bm = joined[(joined['baseMean'].isnull() == True) & 
                (joined['baseMean_b'].isnull() == False) & 
                (joined['baseMean_m'].isnull() == False)]
sb = joined[(joined['baseMean'].isnull() == False) & 
                (joined['baseMean_b'].isnull() == False) & 
                (joined['baseMean_m'].isnull() == True)]
sbm = joined[(joined['baseMean'].isnull() == False) & 
                (joined['baseMean_b'].isnull() == False) & 
                (joined['baseMean_m'].isnull() == False)]

venn3(subsets = (len(just_s), len(just_b), len(sb), len(just_m), len(sm), len(bm), len(sbm)), 
      set_labels = ('Schizophrenia', 'Bipolar Disorder', 'Major Depression'))