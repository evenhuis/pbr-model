import matplotlib.pyplot as plt
import numpy as np
import co2sys 
from scipy.integrate import ode
from scipy.interpolate import interp1d
import pandas as pd

# read in the driver inputs
O2 = np.genfromtxt('O2.txt')
pH= np.genfromtxt('pH.txt')

DIC_Alk = np.genfromtxt('DIC_Alk.txt')

air_p=np.genfromtxt('air.txt')
air = interp1d(air_p[:,0], air_p[:,1], kind='previous',fill_value=1,bounds_error=False)

I_p=np.genfromtxt('I.txt')
I_p[:,1] = I_p[:,1]/2000.
I = interp1d(I_p[:,0], I_p[:,1], kind='linear',fill_value=0,bounds_error=False)

# set up the sea water 
T=27.
TK = T+273.15
S=33.

# calculate the equilibrium CO2 and O2 concentrations
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6

print("O2H = ", O2H)
print("CO2H = ", CO2H)

P   = 600*24
R   = P*0.1
kla = np.log(2.)/4.*60.*24
KM  = 750

print("P = ", P)
print("kLA = ", kla)

Alk = 1600
DIC=Alk*0.95
RR=16/106*0.2
datain = np.array([[S,T,0,0,0,0,0,Alk,DIC]])

#def HCO3_approx( DIC,Alk):
#    return  -818.0204 +   1.9388*DIC  +  0.0681*Alk  -0.0003*DIC*Alk

def HCO3_approx(DIC, Alk, S, T):
    return np.exp(-13.723 + 0.043992*DIC+ 0.0082189*S + 0.0080095*Alk + 1.8925*np.log(DIC) + 
       0.77001*np.log(T) + 1.8896e-06*DIC*Alk - 0.0051115*DIC*np.log(DIC) - 
       0.00074518*DIC*np.log(T) - 0.0017532*Alk*np.log(DIC) + 0.00028851*Alk*np.log(T))

#def CO2_approx( DIC, Alk ):
#    return np.exp(-10.6788 +   0.0095961*DIC +   0.14787*T +   0.033014*S  -0.0015027*Alk   -6.2777e-05*DIC*T   - 8.778e-07*DIC*Alk)

def CO2_approx( DIC, Alk, S, T):
    return np.exp(35.964 + 0.00030677*DIC*T + 2.2504e-06*DIC*Alk + 0.00024263*DIC*S + 4.1513e-05*T*Alk + 7.7881*np.log(T*S) - 1.2232e-07*DIC*T*Alk - 9.3089e-06*DIC*T*S - 7.7838e-08*DIC*S*Alk - 6.9156*np.log(T*S*Alk) + 3.0376e-09*DIC*T*S*Alk)



def pH_approx( DIC, Alk ):
    return 12.26 -0.0030605*DIC -0.043752*T -0.013625*S+ 0.00011315*Alk+ 1.3463e-05*DIC*T + 5.2215e-07*DIC*Alk


#def pH_approx( DIC, Alk ):
#    return 9.3004 - 0.0025878*DIC - 0.030114*T - 0.013622*S + 0.0020192*Alk + 1.3648e-05*DIC*T + 5.1449e-07*DIC*Alk - 5.3084e-06*T*Alk - 1.7026e-07*DIC**2 - 3.2106e-07*Alk**2



def dydt( t, y, *args ):
    O,DIC,Alk= y 
    HCO3  =  HCO3_approx( DIC, Alk, S, T )
    CO2   =  CO2_approx( DIC, Alk, S, T)

    MM = HCO3/(HCO3+KM)

    dO   =  MM*P*I(t)-R       +kla*air(t)*(O2H -  O)
    dDIC = -MM*P*I(t)+R       +kla*air(t)*(CO2H -CO2)
    dAlk  = +2*RR*MM*P*I(t)
    return [dO,dDIC,dAlk]

y0 = [210.,DIC,Alk]; t0,t1=0,10



time = np.linspace(0.5,2.5,501)
ys = np.zeros([len(y0),len(time)])

r = ode(dydt).set_integrator('dopri5')
r.set_initial_value(y0, time[0])

ys[:,0] = y0
for i,t in enumerate(time[1:]):
    r.integrate(t)
    ys[:,i+1]=r.y

datain=np.tile( [S,T,0,0,0,0,0,Alk,DIC],[len(ys[0]),1])
datain[:,8]=ys[1]
dataout=co2sys.CO2sys(datain,10)

#df = pd.read_csv('exact.csv', index_col='time')


fig,axs = plt.subplots(7,1,sharex=True)

axs[0].set_xlim(min(time),max(time))
axs[0].plot(O2[:,0],O2[:,1],color='red')
axs[0].plot(time,ys[0],color='blue')
axs[0].set_ylabel('O2')

#axs[1].plot(time,df['DIC'])
axs[1].plot(DIC_Alk[:,0],DIC_Alk[:,1],'.',color='red')
axs[1].plot(time,ys[1],color='blue')
axs[1].set_ylabel('DIC')

axs[2].plot(DIC_Alk[:,0],DIC_Alk[:,2],'.',color='red')
axs[2].plot(time,ys[2],color='blue')
axs[2].set_ylabel('Alk')

#axs[2].plot(time,df['CO2'])
axs[3].plot(time,CO2_approx( ys[1],Alk, S, T),color='blue')
axs[3].plot(time,dataout['CO2'])
axs[3].set_ylabel('CO2')

#axs[3].plot(time,df['HCO3'])
axs[4].plot(time,HCO3_approx( ys[1],Alk, S, T),color='blue')
axs[4].plot(time,dataout['HCO3'])
axs[4].set_ylabel('HCO3')



#axs[4].plot(time,df['CO3'])
axs[5].plot(time,ys[1]-HCO3_approx( ys[1],Alk, S, T)-CO2_approx( ys[1],Alk, S, T))
axs[5].plot(time,dataout['CO3'])

#axs[5].plot(time,df['pH'])
axs[6].plot(pH[:,0],pH[:,1],color='red')
axs[6].plot(time,pH_approx( ys[1],Alk),color='blue')
axs[6].plot(time,dataout['pH'],color='green')
axs[6].set_ylabel('pH')
#plt.legend(['not sure', 'approx','exact'], loc=2)
plt.legend(['exact-exact', 'approx-approx', 'approx-exact'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()




