import matplotlib.pyplot as plt
import numpy as np

expid = '004'
profile_number = '900'

p1,dtdt1,dqdt1=np.loadtxt('../data_out/conv_tendencies_kuo65_'+expid+'.dat',unpack=True)
p2,dtdt2,dqdt2=np.loadtxt('../data_out/conv_tendencies_kuogeleyn_'+expid+'.dat',unpack=True)
p3,dtdt3,dqdt3=np.loadtxt('../data_out/conv_tendencies_kuosymmetric_'+expid+'.dat',unpack=True)
p4,dtdt4,dqdt4=np.loadtxt('../data_out/conv_tendencies_kuoanthes_'+expid+'.dat',unpack=True)


fig1 = plt.figure()
fig1.set_size_inches(10.,10.)
ax1 = fig1.add_subplot(221)
ax1.set_ylim(1000,100)
#ax1.set_xlim(200.,300)
ax1.plot(dtdt1,p1/100.0,label='Kuo65')
ax1.plot(dtdt2,p1/100.0,label='Kuo-Geleyn')
ax1.plot(dtdt3,p1/100.0,label='Kuo symmetric')
ax1.plot(dtdt4,p1/100.0,label='Kuo Anthes')
ax1.set_title('Heating rate Exp n°'+expid+' --- Profile n°'+profile_number,fontsize='small',loc='left')
ax1.set_xlabel('K/day',fontsize='small')
ax1.set_ylabel('Pressure (hPa)')
ax1.legend(loc='best',fontsize='small')
#
ax2 = fig1.add_subplot(222)
ax2.set_ylim(1000,100)
#ax2.set_xlim(0.,10)
ax2.plot(dqdt1,p1/100.0,label='Kuo65')
ax2.plot(dqdt2,p1/100.0,label='Kuo-Geleyn')
ax2.plot(dqdt3,p1/100.0,label='Kuo symmetric')
ax2.plot(dqdt4,p1/100.0,label='Kuo Anthes')
ax2.set_title('Moistening rate n°'+expid,fontsize='small',loc='left')
ax2.set_xlabel('(g/kg)/day',fontsize='small')
ax2.legend(loc='best',fontsize='small')
#
p,t,tw,tc,qv,qw,qc,gz,mse=np.loadtxt('../data_out/atmospheric_profiles_'+expid+'.dat',unpack=True)

ax3 = fig1.add_subplot(223)
ax3.set_ylim(1000,250)
ax3.set_xlim(240.,305.)
ax3.plot(t*(1+0.608*qv),p,label='Tv')
ax3.plot(tw,p,label='Tw')
ax3.plot(tc*(1+0.608*qc),p,label='Tvc')
ax3.set_title('Temperature',fontsize='small',loc='left')
ax3.set_xlabel('K',fontsize='small')
ax3.set_ylabel('Pressure (hPa)')
ax3.legend(loc='best',fontsize='small')
#
ax4 = fig1.add_subplot(224)
ax4.set_ylim(1000,100)
ax4.plot(qv*1000,p,label='Q')
ax4.plot(qw*1000,p,label='Qw')
ax4.plot(qc*1000,p,label='Qc')
ax4.set_title('Humidity',fontsize='small',loc='left')
ax4.set_xlabel('g/kg',fontsize='small')
ax4.legend(loc='best',fontsize='small')
figure = plt.gcf()
figure.set_size_inches(8, 8)
plt.savefig('../plots/Figure_convective_profiles_'+expid+'.png',dpi=600)
plt.show()
