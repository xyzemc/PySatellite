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
               tlapm = sline[3]
               tsum = sline[4]
               ch = sline[2]
               if ( str(ch) == str(channel)):
                 #print ('ch =', ch, 'channel =', channel)
                 #print('instrument is ', variable,tlapm,tsum)
                 tlapm_coef[mhs_line,0]=tlapm
                 tsum_coef[mhs_line,0]=tsum
                 #print('mhs_line = ', mhs_line)
                 pred1_line = data[nline+1].split()
                 pred2_line = data[nline+2].split()
                 pred=pred1_line+pred2_line
                 print('pred = ',pred)
                 #pred_coef[mhs_line,:,i]=pred
                 lst.append(pred)
                 #lst[mhs_line,i] = pred
                 #print('the length of lst is :', len(lst))
                 #print('file', file,'lst = ',lst[channel-1])
                 mhs_line=mhs_line+1
            nline=nline+1
    dfs=lst


    #concatenated_df = pd.concat(df)

    #return concatenated_df
    return dfs 

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
    plt.plot(cycles,plot_list,linewidth=2, marker ='.',color='red',label='bias0')
    plt.plot(cycles,plot_list_2,linewidth=2, marker ='.',color='blue',label='biasg')
    plt.xticks(np.arange(0, len(cycles)+1, 1),fontsize=8,rotation=45)
    #plt.xticks(x, cycles,fontsize=8,rotation=45)
    plt.title('coef_' + str(coef+1),fontsize=14)
    plt.grid()
    plt.legend()
    img=exp + '_' + instrument + '_ch' + str(channel) + '_coef'  + str(coef+1) + '_' + cycles[0]  + '.png'
    plt.savefig(img)
    plt.close('all')



statdir = '/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias'
print('statdir=',statdir)
coef=0
instrument='amsua_n19'
channel=1
# Grabs all subdirectories of stats files (in date subdirectory)

# Create dictionary of all stats files sorted by date
#satbias = glob.glob(statdir + 'rrfs.spinup.*satbias')
satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.spinup.*satbias")
satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/biasg/satbias/rrfs.spinup.*satbias")
#satbias = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/bias0/satbias/rrfs.prod.*satbias")
#satbias_2 = glob.glob("/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/EXP/nwges/para/biasg/satbias/rrfs.prod.*satbias")
print('satbias =', satbias)
cycles = []
for cyc in satbias:
    cycle=cyc[-18:-8]
    cycles.append(cycle)
    exp=cyc[-23:-19]
print('total cycels =', len(cycles),cycles)
# # Concatenate all files into one dataframe
satbias_files=sorted(satbias,reverse=False)
satbias_files_2=sorted(satbias_2,reverse=False)
cycles_seq=sorted(cycles,reverse=False)
satbias_df = concatenate_dfs(satbias_files, instrument ,channel, cycles)
satbias_df_2 = concatenate_dfs(satbias_files_2, instrument ,channel, cycles)
print(' sorted satebias files =',satbias_files) 
print(' sorted cycles =', cycles_seq)
dfs_slice = list(zip(*satbias_df))
dfs_slice_2 = list(zip(*satbias_df_2))
#print("column[7] = {}".format(dfs_slice[7]))
_plot_bias_rmse_timeseries(exp,cycles_seq,dfs_slice,dfs_slice_2,coef,instrument,channel)

