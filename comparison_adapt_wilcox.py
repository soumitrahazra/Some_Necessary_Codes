import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyPLUTO as pp
import seaborn as sns

da=pd.read_csv('adapt_201811061200_2a')
db=pd.read_csv('GONG_2210A')

xas=da['Col3']
yas=db['Col3']

#pa= abs(xas)/abs(yas);
pa=xas/yas
l_max=40;
subs= np.zeros((41,41))                                                                                                                                                                            
moda=0                                                                                                                                                                                             

for l in range(1, l_max+1): 
   for m in range(0, l+1): 
        subs[l,m]=pa[moda] 
        moda=moda+1 
cmap = sns.diverging_palette(220, 10, as_cmap=True)
# Generate a mask for the upper triangle
#mask = np.triu(np.ones_like(corr, dtype=np.bool))

mask = np.zeros_like(subs) 
mask[np.triu_indices_from(mask)] = True

with sns.axes_style("white"): 
       f, ax = plt.subplots(figsize=(10, 10)) 
       #ax = sns.heatmap(subs, mask=mask, vmax=50, vmin=-50, square=True, cmap=cmap)
       ax = sns.heatmap(subs, mask=mask, square=True)
       ax.set_xlabel('m')
       ax.set_ylabel('l')
       ax.set_title('WSO to Adapt (Real Coefficients)')

plt.savefig('fig6a.png', dpi=500)
plt.show()








