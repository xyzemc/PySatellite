import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

def concatenate_dfs(files, variable, channel,cycles):
    """
    Reads in list of files and creates a concatenated
    dataframe of all cycles.

    Args:
        files : (list) list of files
        variable : (str) data obs type i.e. t, uv, q etc.
        cycles : (list) list of str cycles
        data_type : (str) the input data type i.e conventional
                    or radiance
    Returns:
        concatenated_df : concatenated dataframe over cycles
    """
    dfs = []
    lst = []
    lst_count = []
    dim=len(cycles)
    #print('cycles=',cycles,dim)
    columns = ['coef1','coef2','coef3','coef4','coef5','coef6','coef7','coef8','coef9','coef10','coef11','coef12']
    for i, file in enumerate(files):

        mhs_line=0
        nline=0
        tlapm_coef=np.ones((9708,dim))
        tsum_coef=np.ones((9708,dim))
        #pred_coef=np.zeros(((5,12,dim)))
        pred_coef = []
        f = open (file)
        data=f.readlines()
        for line in data:
            sline = line.split()
            #print('line  =', line)
            if ( variable in sline ):
               #print('sline =',sline)
               ch = sline[2]
               if ( str(ch) == str(channel)):
                 data_count = sline[3]
                 pred1_line = data[nline+1].split()
                 pred2_line = data[nline+2].split()
                 pred=pred1_line+pred2_line
                 print('pred = ',pred)
                 lst.append(pred)
                 lst_count.append(data_count)
                 mhs_line=mhs_line+1
            nline=nline+1
    dfs=lst


    #concatenated_df = pd.concat(df)

    #return concatenated_df
    return dfs,lst_count 

#def _plot_bias_rmse_timeseries(df, config, outdir):
def _plot_bias_rmse_timeseries(exp,cycles,df,df2,coef,instrument,channel):
    """
    Used indexed df to plot rmse and bias.
    """
    plot_list = list(np.float_(df[coef]))
    plot_list_2 = list(np.float_(df2[coef]))
    print('plot_list =', plot_list)
    print('plot_cycles =',cycles)
    x=np.arange(len(cycles))

    plt.figure(figsize=(12,6))
    plt.plot(cycles,plot_list,linewidth=2, marker ='.',color='green')
    #plt.plot(cycles,plot_list,linewidth=2, marker ='.',color='red',label='bias0')
    #plt.plot(cycles,plot_list_2,linewidth=2, marker ='.',color='blue',label='biasg')
    plt.xticks(np.arange(0, len(cycles), 7),fontsize=8,rotation=45)
    #plt.xticks(x, cycles,fontsize=8,rotation=45)
    plt.title('pc_' + str(coef+1),fontsize=14)
    plt.grid()
    plt.legend()
    img=exp + '_' + instrument + '_ch' + str(channel) + '_pc'  + str(coef+1) + '_' + cycles[0]  + '.png'
    plt.savefig(img)
    plt.close('all')

def _plot_bias_rmse_timeseries_count(exp,cycles,df,df2,coef,instrument,channel):
    """
    Used indexed df to plot rmse and bias.
    """
    #plot_list = df
    #print('type of plot_list = ',type(plot_list)
    #plot_list = []
    #convert exponent to float
    #for i in np.arange(len(df)):
    #    test_list = float(df[i])
    #    plot_list.append(test_list)

    #plot_list_2 = df2
    plot_list = list(np.float_(df))
    plot_list_2 = list(np.float_(df2))
    print('plot_list =', plot_list)
    print('plot_cycles =',cycles)
    x=np.arange(len(cycles))

    plt.figure(figsize=(12,6))
    plt.plot(cycles,plot_list,linewidth=2, marker ='.',color='green')
    #plt.plot(cycles,plot_list,linewidth=2, marker ='.',color='red',label='bias0')
    #plt.plot(cycles,plot_list_2,linewidth=2, marker ='.',color='blue',label='biasg')
    plt.xticks(np.arange(0, len(cycles), 7),fontsize=8,rotation=45)
    #plt.xticks(x, cycles,fontsize=8,rotation=45)
    plt.title('Data Count',fontsize=14)
    plt.grid()
    plt.legend()
    #img=exp + '_' + instrument + '_ch' + str(channel) + '_coef'  + str(coef+1) + '_' + cycles[0]  + '.png'
    img=exp + '_' + instrument + '_ch' + str(channel) + '_DataCount' '_' + cycles[0]  + '.png'
    plt.savefig(img)
    plt.close('all')



statdir = '/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias'
print('statdir=',statdir)
coef=7
instrument='mhs_metop-b'
channel=1
# Grabs all subdirectories of stats files (in date subdirectory)

# Create dictionary of all stats files sorted by date
#satbias = glob.glob(statdir + 'rrfs.spinup.*satbias')
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.spinup.*satbias")
#satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/biasg/satbias/rrfs.spinup.*satbias")
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.prod.*satbias")
#satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/biasg/satbias/rrfs.prod.*satbias")
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.spinup.*satbias_pc")
#satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/biasg/satbias/rrfs.spinup.*satbias_pc")
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.prod.*satbias_pc")
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.prod.*satbias_pc")
satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/radiag/fv3sar/satbias_Jul10-Jul20/fv3sar.*_satbias_pc")
satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/radiag/fv3sar/satbias_Jul20-Jul31/fv3sar.*_satbias_pc")
print('satbias =', satbias)
cycles = []
for cyc in satbias:
    cycle=cyc[-21:-11]
    cycles.append(cycle)
    exp=cyc[-26:-22]
print('total cycels =', len(cycles),cycles)
# # Concatenate all files into one dataframe
satbias_files=sorted(satbias,reverse=False)
satbias_files_2=sorted(satbias_2,reverse=False)
cycles_seq=sorted(cycles,reverse=False)
satbias_df,obs_count = concatenate_dfs(satbias_files, instrument ,channel, cycles)
satbias_df_2,obs_count_2 = concatenate_dfs(satbias_files_2, instrument ,channel, cycles)
print(' sorted satebias files =',satbias_files) 
print(' sorted cycles =', cycles_seq)
dfs_slice = list(zip(*satbias_df))
dfs_slice_2 = list(zip(*satbias_df_2))
#print("column[7] = {}".format(dfs_slice[7]))
_plot_bias_rmse_timeseries(exp,cycles_seq,dfs_slice,dfs_slice_2,coef,instrument,channel)
_plot_bias_rmse_timeseries_count(exp,cycles_seq,obs_count,obs_count_2,coef,instrument,channel)

