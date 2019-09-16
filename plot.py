
import matplotlib
matplotlib.rc('font',size=16)
matplotlib.rc('lines',linewidth=3)
import matplotlib.pyplot as plt

from netCDF4 import Dataset

nc = Dataset('./work/example_output.nc')
data = {k:nc.variables[k][:] for k in nc.variables.keys()}
nc.close()

plt.figure(figsize=[8,8])
l1,=plt.semilogx(data['source']*1E6,data['lyZ'],'k-')
l2,=plt.semilogx(data['sink_nox']*1E6,data['lyZ'],'r--')
l3,=plt.semilogx(data['sink_ox']*1E6,data['lyZ'],'b--')
plt.xlabel('rate (ppv/day)')
plt.ylim([25,40])
plt.xlim([1E-2,10])
plt.legend([l1,l2,l3],['source','nox','ox'])
plt.savefig('siracha_balance.png')

