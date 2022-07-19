####################################################################
###This code will plot the time series at bar. It needs to read in 
the input as csv format, which is generated from pygsi/LAMDA.
####################################################################
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def _plot_bias_rmse_timeseries(cycles,df1,df2,instrument,satellite,channel):
    """
    Used indexed df to plot rmse and bias.
    """
    #plot_list = df[coef]
    plot_list = list(np.float_(df1))
    plot_list2 = list(np.float_(df2))
    if len(plot_list) == len(plot_list2):
       print('plot_list =', plot_list)
       print('plot_list2 =', plot_list2)

       plt.figure(figsize=(12,6))
       #plt.plot(cycles,plot_list,'r--')
       x=np.arange(len(cycles))
       #plt.plot(x,plot_list,'r--')
       #plt.plot(x,plot_list2,'b--')
       plt.bar(x-0.2,plot_list,width=0.4, label = 'bias_g')
       plt.bar(x+0.2,plot_list2,width=0.4,label = 'bias_gl')
       #plt.xticks(np.arange(0, len(cycles)+1, 10),fontsize=8,rotation=45)
       #plt.xticks(cycles,fontsize=8,rotation=45)
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel('Data Count')
       plt.title(instrument + '_' + satellite + ' ch.' + str(channel) +' assimilated data count',fontsize=14)
       plt.legend()
       img=instrument + '_' + satellite + '_ch_' +  str(channel) + '_datacount.png'
       plt.savefig(img)
       plt.close('all')
    else:
       print('two data length not equal')


####Read in data
instrument='amsua'
satellite='n19'
ch_number=15
#1
df1 = pd.read_csv('Bias_g.csv',delim_whitespace=True)
datakeys = df1.keys();
print(datakeys)
dfdf1=df1[(df1.instrument == instrument)&(df1.satellite == satellite)]
dfkeys=dfdf1.keys();
print(dfkeys)
#2
df2 = pd.read_csv('Bias_gl.csv',delim_whitespace=True)
datakeys = df2.keys();
print(datakeys)
dfdf2=df2[(df2.instrument == instrument)&(df2.satellite == satellite)]

for channel in np.arange(ch_number):
    print('channel =', channel)
    channel_to_be_plot=channel+1
    if [dfdf1.channel == channel+1]:
       dataframe1=dfdf1[(dfdf1.it == 3)&(dfdf1.channel == channel+1)]
       print(dataframe1)
#2
       dataframe2=dfdf2[(dfdf2.it == 3)&(dfdf2.channel == channel+1)]
       print(dataframe2)
       cycles = []
       for cycle in dataframe2.date:
           cycles.append(cycle)
       print('total cycels =', cycles)
       _plot_bias_rmse_timeseries(cycles,dataframe1.nassim,dataframe2.nassim,instrument,satellite,channel=channel_to_be_plot)
