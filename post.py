import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import pandas as pd
import numpy as np
plt.style.use('ggplot')

absmax = lambda x: max(x.min(), x.max(), key=abs)
absmin = lambda x: min(x.min(), x.max(), key=abs)
size = lambda x: max(x)-min(x)

df = pd.read_pickle('fbd')

df['q1'] = abs(df['RC3b']) > 248
df['q2'] = abs(df['RC3j']) > 248
df['RC>248'] = df['q1'] | df['q2']

df['q3'] = abs(df['RC3b']) < 248
df['q4'] = abs(df['RC3j']) < 248
df['chart'] =  df['q3'] & df['q4']

##print(df.pivot_table(index=['alfa','beta'],
##                     columns=['m_load'],
##                     aggfunc=absmax, # does this do what i think it does?
##                     values=['RC3b','RC3j']))
##print()
##print(df.pivot_table(index=['m_load','alfa'],
##                     columns=['beta'],
##                     aggfunc=absmax, # does this do what i think it does?
##                     values=['RC3b','RC3j','RC>248']))
##print()


##df.pivot_table(index=['alfa','beta'],
##                     columns=['m_load'],
##                     aggfunc=absmin, # does this do what i think it does?
##                     values=['reachx','reachz','chart']).plot.scatter('reachx','reachz',c='chart',subplots=True)

dfi = pd.read_pickle('fbdi')
##print(dfi.pivot_table(
##    index='beta',
##    columns='alfa',
##    aggfunc=np.min,
##    values='m_load'))
# discretize data
DSCRT = 5
dfi['m_load2'] = dfi['m_load'] // DSCRT * DSCRT
dfipt = dfi[dfi.reachz>-5].pivot_table(
    index=['alfa','beta'],
    aggfunc=np.min,
    values=['reachx','reachz','m_load2','m_load'])
##dfipt.plot.scatter(
##        'reachx',
##        'reachz',
##        c='m_load2',
##        s=100,
##        cmap=plt.get_cmap('Spectral_r',size(dfipt['m_load2'])/DSCRT))

x, y, z = np.array([dfipt.reachx, dfipt.reachz, dfipt.m_load2])
n = 100
d = 0
xi = np.linspace(min(dfipt.reachx)+d,max(dfipt.reachx)-d,n)
yi = np.linspace(min(dfipt.reachz)+d,max(dfipt.reachz)-d,n)

zi = griddata(x,y,z,xi,yi,interp='linear')
##plt.contour(xi,yi,zi, 6, linewidths=.5, colors='k')
plt.contour(xi, yi, zi)
# plt.contourf(xi, yi, zi, 4,vmax=abs(zi).max(), vmin=-abs(zi).max())
plt.colorbar()
plt.show()
