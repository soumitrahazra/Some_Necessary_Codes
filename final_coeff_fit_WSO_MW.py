import numpy as np
from astropy.io import fits
import scipy.special as scisp
from scipy import interpolate
import matplotlib.pyplot as plt
import sys, os
import glob
from scipy.io import readsav
from sklearn.decomposition import PCA

data_wilcox_coef= np.loadtxt("datacoef_wilcox.txt")
data_mount_coef= np.loadtxt("datacoef_mount.txt")

wtime= data_wilcox_coef[:,0]
wrcoeff10= data_wilcox_coef[:,1]
wrcoeff11= data_wilcox_coef[:,2]
wicoeff11= data_wilcox_coef[:,3]
wrcoeff20= data_wilcox_coef[:,4]
wrcoeff21= data_wilcox_coef[:,5]
wicoeff21= data_wilcox_coef[:,6]
wrcoeff22= data_wilcox_coef[:,7]
wicoeff22= data_wilcox_coef[:,8]


mtime= data_mount_coef[:,0]
mrcoeff10= data_mount_coef[:,1]
mrcoeff11= data_mount_coef[:,2]
micoeff11= data_mount_coef[:,3]
mrcoeff20= data_mount_coef[:,4]
mrcoeff21= data_mount_coef[:,5]
micoeff21= data_mount_coef[:,6]
mrcoeff22= data_mount_coef[:,7]
micoeff22= data_mount_coef[:,8]


fit1= np.polyfit(wrcoeff10[0:450], mrcoeff10[0:450], 1)
fit2= np.polyfit(wrcoeff20[0:450], mrcoeff20[0:450], 1)
fit3= np.polyfit(wrcoeff11[0:450], mrcoeff11[0:450], 1)
fit4= np.polyfit(wicoeff11[0:450], micoeff11[0:450], 1)
print(fit1, fit2, fit3, fit4)

fit1_fun= np.poly1d(fit1)
fit2_fun= np.poly1d(fit2)
fit3_fun= np.poly1d(fit3)
fit4_fun= np.poly1d(fit4)
print(fit1_fun)
mwo_new_g10=  mrcoeff10[0:450]/abs(fit1[0])
mwo_new_g20=  mrcoeff20[0:450]/abs(fit2[0])
mwo_new_g11= mrcoeff11[0:450]/fit3[0]
mwo_new_h11= wicoeff11[0:450]/abs(fit4[0])

data1= np.array([wrcoeff20[0:450], mrcoeff20[0:450]])
datas=data1.transpose()
print(data1.shape)
print(datas.shape)
pca= PCA(n_components=2)
pca.fit_transform(datas)
print( pca.explained_variance_ratio_)
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('number of components')
plt.ylabel('cumulative explained variance')
plt.show()

#Figure 1
plt.subplot(4,1,1)
plt.scatter(wrcoeff10[0:450], mrcoeff10[0:450])
plt.plot(wrcoeff10[0:450], fit1_fun(wrcoeff10[0:450]), '-k')
#plt.text(fit1)
plt.ylabel('MWO g10')
plt.xlabel('WSO g10')
plt.subplot(4,1,2)
plt.scatter(wrcoeff20[0:450], mrcoeff20[0:450])
plt.plot(wrcoeff20[0:450], fit2_fun(wrcoeff20[0:450]), '-k')
plt.ylabel('MWO g20')
plt.xlabel('WSO g20')
plt.subplot(4,1,3)
plt.scatter(wrcoeff11[0:450], mrcoeff11[0:450])
plt.plot(wrcoeff11[0:450], fit3_fun(wrcoeff11[0:450]), '-k')
plt.ylabel('MWO g11')
plt.xlabel('WSO g11')
plt.subplot(4,1,4)
plt.scatter(wicoeff11[0:450], micoeff11[0:450])
plt.plot(wicoeff11[0:450], fit4_fun(wicoeff11[0:450]), '-k')
plt.ylabel('MWO h11')
plt.xlabel('WSO h11')
plt.tight_layout()

plt.show()


# figure 2
plt.subplot(4,1,1)
plt.plot(wtime[0:450], wrcoeff10[0:450], '-b')
plt.plot(wtime[0:450], mwo_new_g10[0:450], '-r')
plt.ylabel('g10')
plt.legend(['WSO', 'MWO'])
plt.title('MWO scaled to WSO values')
plt.subplot(4,1,2)
plt.plot(wtime[0:450], wrcoeff20[0:450], '-b')
plt.plot(wtime[0:450], mwo_new_g20[0:450], '-r')
plt.ylabel('g20')
plt.subplot(4,1,3)
plt.plot(wtime[0:450], wrcoeff11[0:450], '-b')
plt.plot(wtime[0:450], mwo_new_g11[0:450], '-r')
plt.ylabel('g11')
plt.subplot(4,1,4)
plt.plot(wtime[0:450], wicoeff11[0:450], '-b')
plt.plot(wtime[0:450], mwo_new_h11[0:450], '-r')
plt.ylabel('h11')
plt.xlabel('CR Number')
plt.show()






