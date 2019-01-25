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
T=15.
TK = T+273.15
S=15.

# calculate the equilibrium CO2 and O2 concentrations
O2H  = co2sys.K0_O2 (TK,S)*co2sys.rho_sw(TK,S)/1000 * 0.2094*1e6
CO2H = co2sys.K0_CO2(TK,S)*co2sys.rho_sw(TK,S)/1000 * 500e-6*1e6

print("O2H = ", O2H)
print("CO2H = ", CO2H)

P   = 500*24
R   = 100*24
kla = 200
KM  = 125

print("P = ", P)
print("kLA = ", kla)

Alk = 1600
DIC=Alk*1.75
RR=16/106*0.5
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


def iterate_CO2sys( DIC,Alk,niter=2 ):
    HCO3  =  HCO3_approx( DIC, Alk, S, T )
    CO2   =  CO2_approx( DIC, Alk, S, T)

    TK = T+273.15
    TS = co2sys.T_S( S ) ; KS = co2sys.K_S(TK,S)         # pH_F 
    TF = co2sys.T_F( S ) ; KF = co2sys.K_F(TK,S)         # pH_T

    SWS_2_T  = (1. + TS/KS)/( 1. + TS/KS + TF/KF )
    Free_2_T =  1. + TS/KS

    KW = co2sys.K_W( TK,S)                # pH_T
    KB = co2sys.K_B( TK,S )/SWS_2_T
    TB = co2sys.T_B(S)

    K1 = co2sys.K1_H2CO3( TK,S, const=10) 
    K2 = co2sys.K2_H2CO3( TK,S, const=10) 

    pH    =   pH_approx( DIC, Alk)
    for i in range(niter):
        h = 10**(-pH)
        h_free = h/Free_2_T
        f0 =(  DIC*1e-6*(K1*h+2.*K1*K2)/(h*h+K1*h+K1*K2)   \
            - h_free + KW/h                         \
            - Alk*1e-6                                   \
           +TB  /(1.+h/KB) )*1e6

        df0 = (  DIC*1e-6*(K1 +2*K1*K2)/(h**2+K1*h+K1*K2) \
           -DIC*1e-6*(K1*h+2*K1*K2)/(h**2+K1*h+K1*K2)**2*(2*h+K1) \
           -TB *1./(1+h/KB)**2 / KB                         \
           -KW/h**2 - 1./Free_2_T )*1e6  * (-np.log(10.)*10**(-pH))

        pH = pH - f0/df0
    H = 10**-pH
    H2 = H*H
    denom = (H2+K1*H+K1*K2)
    CO2  = (DIC*H2      /denom)
    HCO3 = (DIC*H *K1   /denom)
    CO3  = DIC   *K1*K2/denom  
    return pH,CO2,HCO3,CO3

def dydt( t, y, *args ):
    O,DIC,Alk= y 
  
    HCO3  =  HCO3_approx( DIC, Alk, S, T )
    CO2   =  CO2_approx( DIC, Alk, S, T)
    pH    =   pH_approx( DIC, Alk)
    
    #pH,CO2,HCO3,CO3 = iterate_CO2sys( DIC, Alk )
    datain = datain = np.array([[S,T,0,0,0,0,0,Alk,DIC]])
    dataout=co2sys.CO2sys(datain,const=10)
    CO2  = dataout['CO2'][0]
    HCO3 = dataout['HCO3'][0]

    KM = 100
    MM = HCO3/(HCO3+KM)

    dO   =  MM*P*I(t)*1.5-R       +kla*air(t)*(O2H -  O)
    dDIC = -MM*P*I(t)   +R       +kla*air(t)*(CO2H -CO2)
    dAlk  = +2*RR*MM*P*I(t)
    return [dO,dDIC,dAlk]

y0 = [210.,DIC,Alk]; t0,t1=0,10


time = np.concatenate([np.linspace(0.6,5,1001),air_p[:,0]])
time.sort()
time = time[(0.6<time)&(time<5)]


r = ode(dydt).set_integrator('dopri5')
r.set_initial_value(y0, time[0])

if( True  ):
    ys = np.zeros([len(y0),len(time)])
    ys[:,0] = y0
    for i,t in enumerate(time[1:]):
        r.integrate(t)
        ys[:,i+1]=r.y

datain=np.tile( [S,T,0,0,0,0,0,Alk,DIC],[len(ys[0]),1])
datain[:,8]=ys[1]
datain[:,7]=ys[2]
dataout=co2sys.CO2sys(datain,10)

niter = 2
iter_app = np.array([iterate_CO2sys( dic,alk, niter=niter) for dic,alk in zip(ys[1],ys[2])] ).T
#df = pd.read_csv('exact.csv', index_col='time')

ns = 20
fig,axs = plt.subplots(7,1,sharex=True,figsize=(5,10))

axs[0].set_xlim(min(time),max(time))
#axs[0].plot(O2[:,0],O2[:,1],color='red')
axs[0].plot(time,ys[0],color='blue')
axs[0].set_ylabel('O2')

#axs[1].plot(time,df['DIC'])
#axs[1].plot(DIC_Alk[:,0],DIC_Alk[:,1],'.',color='red')
axs[1].plot(time,ys[1],color='blue')
axs[1].set_ylabel('DIC')

Alk = ys[2]
#axs[2].plot(DIC_Alk[:,0],DIC_Alk[:,2],'.',color='red')
axs[2].plot(time,Alk,color='blue')
axs[2].set_ylabel('Alk')

#axs[2].plot(time,df['CO2'])
axs[3].plot(time,CO2_approx( ys[1],Alk, S, T),color='blue')
axs[3].plot(time, iter_app[1], color='purple')
axs[3].plot(time[::ns],dataout['CO2'][::ns],'.',color='green')
axs[3].set_ylabel('CO2')
axs[3].set_ylim(0,50)

#axs[3].plot(time,df['HCO3'])
axs[4].plot(time,HCO3_approx( ys[1],Alk, S, T),color='blue')
axs[4].plot(time, iter_app[2], color='purple')
axs[4].plot(time[::ns],dataout['HCO3'][::ns],'.',color='green')
axs[4].set_ylabel('HCO3')
axs[4].set_ylim(0,4000)


#axs[4].plot(time,df['CO3'])
axs[5].plot(time,ys[1]-HCO3_approx( ys[1],ys[2], S, T)-CO2_approx( ys[1],ys[2], S, T),color='blue')
axs[5].plot(time, iter_app[3], color='purple')
axs[5].plot(time[::ns],dataout['CO3'][::ns],'.',color='green')
axs[5].set_ylim(0,2000)

#axs[5].plot(time,df['pH'])
#axs[6].plot(pH[:,0],pH[:,1],color='red', label='obs')
axs[6].plot(time, iter_app[0], color='purple', label='iter={}'.format(niter) )
axs[6].plot(time,pH_approx( ys[1],ys[2]),color='blue', label='approx')
axs[6].plot(time[::ns],dataout['pH'][::ns],'.',color='green',label='exact')
axs[6].legend()
axs[6].set_ylabel('pH')

#plt.legend(['not sure', 'approx','exact'], loc=2)
#plt.legend(['exact-exact', 'approx-approx', 'approx-exact'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()




