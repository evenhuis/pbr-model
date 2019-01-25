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

# set up the sea water 
T=27.
TK = T+273.15
S=33.
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6

P   = 800.*24  
R   = 1. *24*0. *0
kla = np.log(2.)/4.*60.*24
KM  = 200.

print("P = ", P)
print("kLA = ", kla)

TA=3000
DIC=500
RR=16/106
datain = np.array([[S,T,0,0,0,0,0,TA,DIC]])


#def dydt( t, y, *args ):
#    O,DIC,TA = y 
#    datain = datain = np.array([[S,T,0,0,0,0,0,TA,DIC]])
#    dataout=co2sys.CO2sys(datain,const=10)
#    CO2  = dataout['CO2'][0]
#    HCO3 = dataout['HCO3'][0]
#    MM = HCO3/(HCO3+KM)
#
#    Pt   =  MM*P*I(t)
#    dO   =  Pt-R       +kla*air(t)*(O2H -  O)
#    dDIC = -Pt+R       +kla*air(t)*(CO2H -CO2)
#    dTA  = +2*RR*Pt
#    return [dO,dDIC,dTA]
#
#y0 = [210.,DIC,TA]; t0,t1=0,10

def dydt( t, y, *args ):
    O,DIC = y 
    datain = datain = np.array([[S,T,0,0,0,0,0,TA,DIC]])
    dataout=co2sys.CO2sys(datain,const=10)
    CO2  = dataout['CO2'][0]
    HCO3 = dataout['HCO3'][0]
    MM = HCO3/(HCO3+KM)

    Pt   =  MM*P*I(t)
    dO   =  Pt-R       +kla*air(t)*(O2H -  O)
    dDIC = -Pt+R       +kla*air(t)*(CO2H -CO2)
#    dTA  = +2*RR*Pt
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

datain=np.tile( [S,T,0,0,0,0,0,TA,DIC],[len(ys[0]),1])
datain[:,8]=ys[1]
#datain[:,7]=ys[2]
dataout=co2sys.CO2sys(datain,10)

fig,axs = plt.subplots(7,1,sharex=True)


ax=axs[0];ax.plot(time,ys[0]); ax.set_ylabel("O2")
ax=axs[1];ax.plot(time,ys[1]); ax.set_ylabel("DIC")
ax=axs[2];ax.plot(time,dataout['CO2']); ax.set_ylabel("CO2")
ax=axs[3];ax.plot(time,dataout['HCO3']); ax.set_ylabel("HCO3")
ax=axs[4];ax.plot(time,dataout['CO3']); ax.set_ylabel("CO3")
ax=axs[5];ax.plot(time,dataout['pH']); ax.set_ylabel("pH")
#ax=axs[6];ax.plot(time,ys[2]); ax.set_ylabel("TA")
plt.show()

import pandas as pd
df = pd.DataFrame(columns=['time','O2','DIC','TA','CO2','HCO3','CO3','pH'])
df['time']=time
df['O2'] =ys[0]
df['DIC']=ys[1]
df['TA'] =TA
df['CO2'] =dataout['CO2']
df['HCO3']=dataout['HCO3']
df['CO3']=dataout['CO3']
df['pH']=dataout['pH']
df.set_index('time')
df.to_csv("exact.csv")



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
