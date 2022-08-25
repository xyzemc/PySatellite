#
# (C) Copyright 2020 UCAR
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#

# This example plots the outputs of Example 5a using Cartopy.
# For this to successfully run, you need Cartopy, Matplotlib and NumPy.

#import ioda
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
import cartopy.crs as ccrs
import netCDF4 as nc
from netCDF4 import Dataset
import scipy.stats as stats
from os.path import exists

if os.environ.get('LIBDIR') is not None:
    sys.path.append(os.environ['LIBDIR'])
#------------------------------------------------------------------
def _read_netcdf_diag(fname,channel):

  # Read input data file
  # ----------------------
  ncd = nc.Dataset(fname, 'r') 

  # NetCDF global attributes
  # --------------------------
  nc_attrs = ncd.ncattrs()
  print('NetCDF Global Attributes: ')
  print('nc_attrs = ', nc_attrs)
  for nc_attr in nc_attrs:
     print('nc_attr', ncd.getncattr(nc_attr))

  # Dimension shape information 
  # -----------------------------
  nc_dims = [dim for dim in ncd.dimensions]  # list of nc dimensions
  print('nc_dims = ', nc_dims)
  for dim in nc_dims:
     print('nc_dims', dim, len(ncd.dimensions[dim]))

  # Variable information
  # ---------------------
  latData = ncd.variables['Latitude'][:]
  lonData = ncd.variables['Longitude'][:]
  satidData = ncd.variables['satinfo_chan'][:]
  ombnbcData = ncd.variables['Obs_Minus_Forecast_unadjusted'][:]
  ombData = ncd.variables['Obs_Minus_Forecast_adjusted'][:]
  channelNumber = ncd.variables['Channel_Index'][:]
  qcflag = ncd.variables['QC_Flag'][:]
  variance = ncd.variables['Inverse_Observation_Error'][:]
  #print('qc_flag =', qcflag)

  # Check data 
  # -----------
  #print('lat   = ', latData)
  #print('lon   = ', lonData)
  #print('omb   = ', ombData)
  #print('satid = ', satidData)
  #print('ombData shape = ', ombData.shape)
  #print('channelNumber = ', channelNumber)
  #print('channelNumber shape = ', channelNumber.shape)

  # Calculate data statistics 
  # --------------------------
  ch_start=min(channelNumber)
  ch_end=max(channelNumber)
  ch=channel
  idx = []
  idx_nbc = []
  print('ch is ', ch_start, ch_end)
  for n in np.arange(len(ombData)):
      if channelNumber[n] == ch and variance[n] != 0.0:
      #if channelNumber[n] == ch and qcflag[n] == 0.0:
      #if channelNumber[n] == ch and qcflag[n] >= 0.0:
      #if channelNumber[n] == ch:
         idx.append(n)
  for nbc in np.arange(len(ombnbcData)):
      if channelNumber[nbc] == ch and variance[nbc] != 0.0:
      #if channelNumber[nbc] == ch and qcflag[nbc] == 0.0:
      #if channelNumber[nbc] == ch and qcflag[nbc] >= 0.0:
      #if channelNumber[nbc] == ch:
         idx_nbc.append(nbc)
  #print('idx=', idx)        
           
  obarray=ombData[idx]
  obnbcarray=ombnbcData[idx_nbc]
  #print('obarray----',len(obarray),obarray)
  #print('obnbcarray----',len(obnbcarray),obnbcarray)
  obs_assim = len(obarray)
  return obarray,obnbcarray,obs_assim

#---------------------------------------------------------
def _plot_bias_rmse_timeseries_omf(exp,exp2,list1,list2,list3,list4,cycles,instrument,channel):
    """
    Used indexed df to plot rmse and bias.
    """
    #plot_list = df[coef]
    plot_list = list(np.float_(list1))
    plot_list2 = list(np.float_(list2))
    plot_list3 = list(np.float_(list3))
    plot_list4 = list(np.float_(list4))
    label_1='Exp: '+exp
    label_2='Exp: '+exp2
    title='Avg Bias Correction \n' + instrument + '_' + ' ch.' + str(channel) 
    img=exp+'-'+exp2+'_bias_'+instrument + '_ch_' +  str(channel) + '.png'
    img2=exp+'-'+exp2+'_ObsNumber_'+instrument + '_ch_' +  str(channel) + '.png'
    ymark='Avg Bias'
    if len(plot_list) == len(plot_list2):
       print('plot_list =', plot_list)
       print('plot_list2 =', plot_list2)

       plt.figure(figsize=(12,10)).subplots(1,2)
       
       x=np.arange(len(cycles))
       plt.subplot(211)
       plt.plot(x,plot_list,'g-',label=label_1)
       plt.plot(x,plot_list3,'y-',label=label_2)
       #plt.ylim(-5, 5)
       plt.ylim(min(plot_list)-1,max(plot_list)+1)
       plt.xlim(min(x), max(x))
       # Horizontal line at 0
       plt.hlines(y= 0, xmin= 0, xmax= max(x), color='black', linestyle ='solid', linewidth = 1)
       plt.axhline(y=0)
       #plt.xticks(np.arange(0, len(cycles)+1, 10),fontsize=8,rotation=45)
       #plt.xticks(cycles[:][1],fontsize=8,rotation=30)
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel(ymark)
       plt.xlabel('Cycle Time')
       plt.title(title,fontsize=14, loc='left')
       #plt.subplots_adjust(bottom=0.25)
       plt.tight_layout()
       plt.legend()

       plt.subplot(212)
       plt.plot(x,plot_list2,'g-',label=label_1)
       plt.plot(x,plot_list4,'y-',label=label_2)
       title='Sdv Bias Correction \n' + instrument + '_' + ' ch.' + str(channel) 
       #plt.ylim(-5, 5)
       plt.ylim(0,max(plot_list2)+1)
       plt.xlim(min(x), max(x))
       # Horizontal line at 0
       plt.hlines(y= 0, xmin= 0, xmax= max(x), color='black', linestyle ='solid', linewidth = 1)
       plt.axhline(y=0)
       #plt.xticks(np.arange(0, len(cycles)+1, 10),fontsize=8,rotation=45)
       #plt.xticks(cycles[:][1],fontsize=8,rotation=30)
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel('Std Bias')
       plt.xlabel('Cycle Time')
       plt.title(title,fontsize=14, loc='left')
       plt.tight_layout()
       plt.legend()
       #plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.55,
       #             wspace=0.35)


       plt.savefig(img)
       plt.close('all')

def _plot_timeseries_obs_number(exp,exp2,obs_number1, obs_number2,cycles,instrument,channel):
    """
    Used indexed df to plot observation number.
    """
    plot_number1 = list(np.float_(obs_number1))
    plot_number2 = list(np.float_(obs_number2))
    label_1='Exp: '+exp
    label_2='Exp: '+exp2
    img2=exp+'-'+exp2+'_ObsNumber_'+instrument + '_ch_' +  str(channel) + '.png'
    if len(plot_number1) == len(plot_number2):
       plt.figure(figsize=(12,6))
       
       x=np.arange(len(cycles))
      
       ##plot observation number

       plt.plot(x,plot_number1,'g-',label=label_1)
       plt.plot(x,plot_number2,'y-',label=label_2)
       title='Number of Observations \n' + instrument + '_' + ' ch.' + str(channel) 
       plt.ylim(0,max(plot_number2+plot_number2)+50)
       plt.xlim(min(x), max(x))
       plt.xticks(x, cycles,fontsize=8,rotation=45)
       plt.ylabel('Obs Count')
       plt.xlabel('Cycle Time')
       plt.title(title,fontsize=14, loc='left')
       plt.tight_layout()
       plt.legend()
       plt.savefig(img2)
       plt.close('all')
                                                     
    
#---------------------------------------------------------


# Set input and output file names
# --------------------------------
channel = 1 #to be read and plot
instrument = 'amsua_n19'
exp='bias0'
exp2='biasgl'
all_omf=np.array([])
all_omf_nbc=np.array([])
statdir="/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/radiag/"+exp+"/"
statdir2="/scratch1/NCEPDEV/stmp2/Xiaoyan.Zhang/radiag/"+exp2+"/"
subdirs = sorted([f.path for f in os.scandir(statdir) if f.is_dir()])
subdirs2 = sorted([f.path for f in os.scandir(statdir2) if f.is_dir()])
cycles = []
cycles2 = []
obs_number = []
obs_number2 = []
bias_list = []
bias_list2 = []
biasstd_list = []
biasstd_list2 = []
for subdir in subdirs:
    cycle = subdir.split('/')[-1]
    print('cycle=',cycle,subdir)
    fname=subdir + '/diag_' +instrument+'_ges.'+cycle+'.nc4'
    file_exist = exists(fname)
    if (file_exist):
       print('fname=',fname)
       omf, omfnbc, data_count = _read_netcdf_diag(fname,channel)
       bias=omf-omfnbc
       bias_avg = np.nanmean(bias) 
       bias_list.append(bias_avg)

       bias_std = np.nanstd(bias)
       biasstd_list.append(bias_std)

       print('dim of bias =',len(bias))
       cycles.append(cycle)
       obs_number.append(data_count)
for subdir in subdirs2:
    cycle = subdir.split('/')[-1]
    print('cycle=',cycle,subdir)
    fname=subdir + '/diag_' +instrument+'_ges.'+cycle+'.nc4'
    file_exist = exists(fname)
    if (file_exist):
       print('fname2=',fname)
       omf, omfnbc, data_count2 = _read_netcdf_diag(fname,channel)
       bias2=omf-omfnbc
       bias2_avg = np.nanmean(bias2) 
       bias_list2.append(bias2_avg)

       bias2_std = np.nanstd(bias2)
       biasstd_list2.append(bias2_std)

       cycles2.append(cycle)
       obs_number2.append(data_count2)

print('dim of cycles =',len(cycles2))
print('dim of omf_list =',len(bias_list2))
if len(cycles) == len(cycles2):
   _plot_bias_rmse_timeseries_omf(exp,exp2,bias_list,biasstd_list,bias_list2,biasstd_list2,cycles,instrument,channel)
   _plot_timeseries_obs_number(exp,exp2,obs_number, obs_number2,cycles,instrument,channel)
else:
   print('time is different ')
   print('cycle of ',exp, 'is ',len(cycles))
   print('cycle of ',exp2, 'is ',len(cycles2))

exit()
