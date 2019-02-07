import numpy as np
import pandas as pd
import datetime as dt

# reference time
t0 = dt.datetime(2017,10,17,0,0,0)

# read in the times and volumes of dilutions
dil    =  pd.read_table("dilution_data.txt"  ,  sep="\s+")
time_rel = np.array( [ (dt.datetime.strptime( ds, "%Y%m%d-%H%M")-t0)/dt.timedelta(days=1) for ds in dil.datetime ])
dil.insert(1,"time_rel",time_rel)


# the time window (convert min to day )
dt = 10. /24./60.

# print the dilution series
print("-100. 0. 0. ")
for i,row in (dil).iterrows():
    if( row['PBR'] ==  4 ):
        print(row['time_rel']-dt/2.,row['vol_F2_added']/dt, row['vol_F2_added'])
        print(row['time_rel']+dt/2.,0., 0. )
    #print("{:.4f}, {:.3f}".format(td-td,vol)) 

