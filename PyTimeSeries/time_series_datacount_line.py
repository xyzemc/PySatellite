import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def _plot_bias_rmse_timeseries(cycles,iteration,df1,df2,instrument,satellite,channel):
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
       #plt.xticks(np.arange(0, len(cycles)+1, 10),fontsize=8,rotation=30)
       #plt.xticks(x, cycles,len(cycles)+1, 10),fontsize=8,rotation=45)
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel('Data Count')
       plt.xlabel('Cycling Time')
       plt.title(instrument + '_' + satellite + ' ch.' + str(channel) +' assimilated data count',fontsize=14)
       #plt.subplots_adjust(bottom=0.25)
       plt.tight_layout()
       plt.legend()
       if iteration == 1:
          img=instrument + '_' + satellite + '_ch_' +  str(channel) + '_datacount_F.png'
       if iteration == 3:
          img=instrument + '_' + satellite + '_ch_' +  str(channel) + '_datacount_A.png'
       plt.savefig(img)
       plt.close('all')
    else:
       print('two data length not equal')

def _plot_bias_rmse_timeseries_omf(cycles,iteration,df1,df2,df3,df4,instrument,satellite,channel):
    """
    Used indexed df to plot rmse and bias.
    """
    #plot_list = df[coef]
    plot_list = list(np.float_(df1))
    plot_list2 = list(np.float_(df2))
    plot_list3 = list(np.float_(df3))
    plot_list4 = list(np.float_(df4))
    if iteration == 1:
       label_1="bias_0:OmF_wobc"
       label_2="bias_g:OmF_wobc"
       label_3="bias_0:OmF_bc"
       label_4="bias_g:OmF_bc"
       title=instrument + '_' + satellite + ' ch.' + str(channel) +' OmB '
       img=instrument + '_' + satellite + '_ch_' +  str(channel) + '_OmB.png'
       ymark='O-B'
    if iteration == 3:
       label_1="bias_0:OmA_wobc"
       label_2="bias_g:OmA_wobc"
       label_3="bias_0:OmA_bc"
       label_4="bias_g:OmA_bc"
       title=instrument + '_' + satellite + ' ch.' + str(channel) +' OmA '
       img=instrument + '_' + satellite + '_ch_' +  str(channel) + '_OmA.png'
       ymark='O-A'
       
    if len(plot_list) == len(plot_list2):
       print('plot_list =', plot_list)
       print('plot_list2 =', plot_list2)

       plt.figure(figsize=(12,4))
       #plt.plot(cycles,plot_list,'r--')
       x=np.arange(len(cycles))
       print('label=label_1=',label_1)
       plt.plot(x,plot_list,'g--',label=label_1)
       plt.plot(x,plot_list2,'y--',label=label_2)
       plt.plot(x,plot_list3,'g',label=label_3)
       plt.plot(x,plot_list4,'y',label=label_4)
       plt.ylim(-5, 5)
       plt.xlim(min(x), max(x))
       # Horizontal line at 0
       #plt.hlines(y= 0, xmin= 0, xmax= max(x), color='black', linestyle ='solid', linewidth = 1)
       plt.axhline(y=0)
       #plt.xticks(np.arange(0, len(cycles)+1, 10),fontsize=8,rotation=45)
       #plt.xticks(cycles[:][1],fontsize=8,rotation=30)
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel(ymark)
       plt.xlabel('Cycling Time')
       plt.title(title,fontsize=14)
       #plt.subplots_adjust(bottom=0.25)
       plt.tight_layout()
       plt.legend()
       plt.savefig(img)
       plt.close('all')


####Read in data
instrument='amsua'
satellite='n19'
ch_number=15
iteration=1
#############
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

#for channel in np.arange(ch_number):
for channel in np.arange(1):
#for channel in np.arange(7,10):
    print('channel =', channel)
    channel_to_be_plot=channel+1
    if [dfdf1.channel == channel+1]:
       dataframe1=dfdf1[(dfdf1.it == iteration)&(dfdf1.channel == channel+1)]
       print(dataframe1)
#2
       dataframe2=dfdf2[(dfdf2.it == iteration)&(dfdf2.channel == channel+1)]
       print(dataframe2)
       cycles = []
       for cycle in dataframe2.date:
           cycles.append(cycle)
       print('total cycels =', cycles)
       _plot_bias_rmse_timeseries(cycles,iteration,dataframe1.nassim,dataframe2.nassim,instrument,satellite,channel=channel_to_be_plot)
       _plot_bias_rmse_timeseries_omf(cycles,iteration,dataframe1.OmF_wobc,dataframe2.OmF_wobc,dataframe1.OmF_bc,dataframe2.OmF_bc,instrument,satellite,channel=channel_to_be_plot)
