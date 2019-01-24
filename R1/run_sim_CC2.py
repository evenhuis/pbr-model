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
T=20.
TK = T+273.15
S=33.
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6

P   = 600.*24
R   = 1. *24*0.
kla = np.log(2.)/4.*60.*24.
KM  = 200.

print(P,kla)

DIC=500
TA = 2400

K1 = co2sys.K1_H2CO3(TK,S,const=10)*1e6 ; RT1=6000
K2 = co2sys.K2_H2CO3(TK,S,const=10)*1e6 ; RT2=12000
KW = co2sys.K_W(TK,S)*1e12              ;

datain = np.array([[S,T,0,0,0,0,0,TA,DIC]])
dataout=co2sys.CO2sys(datain,const=10)
    # O2     DIC   CO2   HCO3 TA
y0 = [210., DIC,*dataout[['CO2','HCO3','CO3']][0],TA ]
def dydt( t, y, *args ):
    import math
    O, DIC,CO2,HCO3,CO3, TA = y 
    MM = HCO3/(HCO3+KM)

    DTA = TA-HCO3-2*CO3

    H   =  -DTA+np.sqrt(DTA*DTA+4*KW)
    R1  =    (K1*RT1* CO2  - RT1* H*HCO3 )
    R2  =    (K2*RT2* HCO3 - RT2* H* CO3 )



    dO    =  MM*P*I(t)-R       +kla*air(t)*(O2H -  O)
    dDIC  = -MM*P*I(t)+R       +kla*air(t)*(CO2H -CO2)
    dCO2  =                    +kla*air(t)*(CO2H -CO2) -R1
    dHCO3 = -MM*P*I(t)+R                               +R1  -R2
    dCO3  =                                                 +R2
    dTA   = 0.
    return [dO,dDIC,dCO2,dHCO3,dCO3,dTA]
t = np.linspace(0,10,10*24*60+1)

dydt( 0,y0)


if( True ): 
    time = np.linspace(3,5.5,201)

    r = ode(dydt).set_integrator('dopri5')
    r.set_initial_value(y0, time[0])
    ys = np.zeros([len(y0),len(time)])
    ys[:,0] = y0
    for i,t in enumerate(time[1:]):
        r.integrate(t)
        ys[:,i+1]=r.y

#datain=np.tile([  S,T, 0., 0., 0., 0., 0., 2300.,DIC],[len(ys[0]),1])
#datain[:,8]=ys[1]
#dataout=co2sys.CO2sys(datain,10)

def plot_run( ys1, df=None ):

    O2,DIC,CO2,HCO3,CO3,TA = ys

    fig,axs = plt.subplots(6,1,sharex=True)

    ax = axs[0]
    ax.plot(time,ys[0])

    ax = axs[1]
    ax.plot(time,DIC)
    ax.axhline(16.3,ls='--')

    ax = axs[2]
    ax.plot(time,CO2)
    #ax.plot(time,ys[3])
    #ax.plot(time,ys[1]+ys[2]+ys[3])

    ax = axs[3]
    ax.plot( time,HCO3)


    ax = axs[4]
    ax.plot(time,CO3)

    DTA = TA-HCO3-2*CO3
    pH = -np.log10((DTA-np.sqrt(DTA*DTA-4*KW))*1E-6)
    ax = axs[5]
    ax.plot(time,pH)

    axs[0].set_ylabel('O2')
    axs[1].set_ylabel('DIC')
    axs[2].set_ylabel('CO2')
    axs[3].set_ylabel('HCO3')
    axs[4].set_ylabel('CO3')
    axs[5].set_ylabel('pH')

    if( df is not None):
        axs[0].plot(df['O2'])
        axs[1].plot(df['DIC'])
        axs[2].plot(df['CO2'])
        axs[3].plot(df['HCO3'])
        axs[4].plot(df['CO3'])
        axs[5].plot(df['pH'])
    plt.show()

import pandas as pd
df = pd.read_csv('exact.csv',index_col='time')

plot_run(ys,df=df)

#ax.plot(time,ys[1]+ys[2]+ys[3])
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
