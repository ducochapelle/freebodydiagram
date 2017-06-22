import pandas as pd
import numpy as np
df = pd.read_pickle('fbd')

print(df.pivot_table(index=['alfa', 'beta'],
                     columns='m_load',
                     aggfunc=lambda x: max(x.min(), x.max(), key=abs), # does this do what i think it does?
                     values=['RC3b','RC3j']))

