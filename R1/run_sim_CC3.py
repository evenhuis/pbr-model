import matplotlib.pyplot as plt
import numpy as np
import co2sys 
from scipy.integrate import ode
from scipy.interpolate import interp1d

O2 = np.genfromtxt('DO.txt')
pH= np.genfromtxt('pH.txt')

air_p=np.genfromtxt('gas.txt')
air = interp1d(air_p[:,0], air_p[:,1], kind='previous',fill_value=1,bounds_error=False)

I_p=np.genfromtxt('I.txt')
I = interp1d(I_p[:,0], I_p[:,1], kind='linear',fill_value=0,bounds_error=False)

# set up the sea water 
T=20.
TK = T+273.15
S=33.
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6

P   = 600.*24  
R   = 1. *24*0. *0
kla = np.log(2.)/4.*60.*24
KM  = 200.

print(P,kla)

Alk = 2400
DIC=500

datain = np.array([[S,T,0,0,0,0,0,2400,DIC]])

def HCO3_approx( DIC,Alk):
    return  -818.0204 +   1.9388*DIC  +  0.0681*Alk  -0.0003*DIC*Alk

def CO2_approx( DIC, Alk ):
    return np.exp(-10.6788 +   0.0096*DIC +   0.1479*T +   0.0330*S  -0.0015*Alk   -6.2777e-05*DIC*T   - 8.778e-07*DIC*Alk)

def pH_approx( DIC, Alk ):
    return 12.26 -0.0030605*DIC -0.043752*T -0.013625*S+ 0.00011315*Alk+ 1.3463e-05*DIC*T + 5.2215e-07*DIC*Alk

def dydt( t, y, *args ):
    O,DIC = y 
    HCO3  =  HCO3_approx( DIC, Alk )
    CO2   =  CO2_approx( DIC, Alk )

    MM = HCO3/(HCO3+KM)

    dO   =  MM*P*I(t)-R       +kla*air(t)*(O2H -  O)
    dDIC = -MM*P*I(t)+R       +kla*air(t)*(CO2H -CO2)
    return [dO,dDIC]

y0 = [210.,DIC]; t0,t1=0,10



time = np.linspace(3,5.5,301)
ys = np.zeros([len(y0),len(time)])

r = ode(dydt).set_integrator('dopri5')
r.set_initial_value(y0, time[0])

ys[:,0] = y0
for i,t in enumerate(time[1:]):
    r.integrate(t)
    ys[:,i+1]=r.y

#datain=np.tile( [S,T,0,0,0,0,0,2400,DIC],[len(ys[0]),1])
#datain[:,8]=ys[1]
#dataout=co2sys.CO2sys(datain,10)

df = pd.read_csv('exact.csv',index_col='time')


fig,axs = plt.subplots(6,1,sharex=True)

axs[0].plot(time,df['O2'])
axs[0].plot(time,ys[0])

axs[1].plot(time,df['DIC'])
axs[1].plot(time,ys[1])

axs[2].plot(time,df['CO2'])
axs[2].plot(time,CO2_approx( ys[1],Alk))
axs[2].plot(time,dataout['CO2'])

axs[3].plot(time,df['HCO3'])
axs[3].plot(time,HCO3_approx( ys[1],Alk))
axs[3].plot(time,dataout['HCO3'])

axs[4].plot(time,df['HCO3'])
axs[4].plot(time,dataout['CO3'])

axs[5].plot(time,df['pH'])
axs[5].plot(time,pH_approx( ys[1],Alk))
axs[5].plot(time,dataout['pH'])
plt.show()



#fig,axs = plt.subplots(3,1,sharex=True)
#ax = axs[0]
#ax.plot(O2[:,0],O2[:,1])
#ax.plot(time,ys[0])
#
#ax = axs[1]
#ax.plot(time,ys[1])
#ax.plot(time,dataout['HCO3'])
#ax.set_xlim(time[0],time[-1])
#
#
#ax = axs[2]
#ax.plot(pH[:,0],pH[:,1])
#ax.plot(time,dataout['pH'])
#plt.show()
