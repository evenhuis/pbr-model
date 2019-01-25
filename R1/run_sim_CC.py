import matplotlib.pyplot as plt
import numpy as np
import co2sys 
from scipy.integrate import ode
from scipy.interpolate import interp1d

# read in the driver inputs
O2 = np.genfromtxt('DO.txt')
pH= np.genfromtxt('pH.txt')

air_p=np.genfromtxt('gas.txt')
air = interp1d(air_p[:,0], air_p[:,1], kind='previous',fill_value=1,bounds_error=False)

I_p=np.genfromtxt('I.txt')
I = interp1d(I_p[:,0], I_p[:,1], kind='linear',fill_value=0,bounds_error=False)

# set up the seawater 
T=27.
TK = T+273.15
S=33.

# calculate the equilibrium CO2 and O2 concentrations
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6 #BM: 500e-6??

print("O2H = ", O2H)
print("CO2H = ", CO2H)

P   = 600.*24
R   = 1. *24*0.
kla = np.log(2.)/4.*60.*24.
KM  = 200.

print("P = ", P)
print("kLA = ", kla)

DIC=500

K1 = co2sys.K1_H2CO3(TK,S,const=10)*1e6 ; RT1=800
K2 = co2sys.K2_H2CO3(TK,S,const=10)*1e6 ; RT2=830
KW = co2sys.K_W(TK,S)*1e12              ; RT3=800

datain = np.array([[S,T,0,0,0,0,0,2400,DIC]])
dataout=co2sys.CO2sys(datain,const=10)
y0 = [210., *dataout[['CO2','HCO3','CO3','pH']][0], 0 ]
y0[4]=10**(-y0[4])*1e6
y0[5]=dataout['OH'][0]*1e6

def dydt( t, y, *args ):
    O, CO2,HCO3,CO3, H,OH  = y 
    MM = HCO3/(HCO3+KM)

    R1  =    (K1*RT1* CO2  - RT1* H*HCO3 )
    R2  =    (K2*RT2* HCO3 - RT2* H* CO3 )
    R3  =    KW*RT3 - RT3*H*OH


    dO    =  MM*P*I(t)-R       +kla*air(t)*(O2H -  O)
    dCO2  =                    +kla*air(t)*(CO2H -CO2) -R1
    dHCO3 = -MM*P*I(t)+R                               +R1 - R2
    dCO3  =                                                + R2
    dH    =                                            +R1 + R2 +R3
    dOH   =                                                     +R3
    return [dO,dCO2,dHCO3,dCO3,dH,dOH]
t = np.linspace(0,10,10*24*60+1)




time = np.linspace(3,4.5,211)
ys = np.zeros([len(y0),len(time)])

r = ode(dydt).set_integrator('dopri5')
r.set_initial_value(y0, time[0])

ys[:,0] = y0
for i,t in enumerate(time[1:]):
    r.integrate(t)
    ys[:,i+1]=r.y

#datain=np.tile([  S,T, 0., 0., 0., 0., 0., 2300.,DIC],[len(ys[0]),1])
#datain[:,8]=ys[1]
#dataout=co2sys.CO2sys(datain,10)

fig,axs = plt.subplots(5,1,sharex=True)
ax = axs[0]
ax.plot(time,ys[1])
ax.axhline(16.3,ls='--')

ax = axs[1]
ax.plot(time,ys[2])
#ax.plot(time,ys[3])
#ax.plot(time,ys[1]+ys[2]+ys[3])

ax = axs[2]
ax.plot( time,ys[3])

ax = axs[3]
ax.plot(time,-np.log10(ys[4]*1e-6))
#ax.plot(time,ys[1]+ys[2]+ys[3])

#ax = axs[4]
#KW = co2sys.K_W(TK,S)
#ax.plot(time,ys[1]+2*ys[2]-(ys[4]-KW*1e-6/ys[4]))
plt.show()

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
