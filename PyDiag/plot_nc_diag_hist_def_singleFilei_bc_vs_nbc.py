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
  print('qc_flag =', qcflag)

  # Check data 
  # -----------
  print('lat   = ', latData)
  print('lon   = ', lonData)
  print('omb   = ', ombData)
  print('satid = ', satidData)
  print('ombData shape = ', ombData.shape)
  print('channelNumber = ', channelNumber)
  print('channelNumber shape = ', channelNumber.shape)

  # Calculate data statistics 
  # --------------------------
  ch_start=min(channelNumber)
  ch_end=max(channelNumber)
  ch=channel
  idx = []
  idx_nbc = []
  print('ch is ', ch_start, ch_end)
  for n in np.arange(len(ombData)):
      if channelNumber[n] == ch and qcflag[n] == 0.0:
         idx.append(n)
  for nbc in np.arange(len(ombnbcData)):
      if channelNumber[nbc] == ch and qcflag[nbc] == 0.0:
         idx_nbc.append(nbc)
  #print('idx=', idx)        
           
  obarray=ombData[idx]
  obnbcarray=ombnbcData[idx_nbc]
  print('obarray----',len(obarray),obarray)
  print('obnbcarray----',len(obnbcarray),obnbcarray)
  return obarray,obnbcarray
    
#---------------------------------------------------------
def _plot_histgram(omf,omfnbc,instrument,channel): 
  nbins=50
  omf_density = stats.gaussian_kde(omf)
  _max=np.max(np.abs(omf))
  omf_bins = np.linspace(-_max, _max, nbins)
  print('omf_bins=',omf_bins)

  stdev = np.nanstd(omf)  # Standard deviation
  omean = np.nanmean(omf) # Mean of the data
  datmi = np.nanmin(omf)  # Min of the data
  datma = np.nanmax(omf)  # Max of the data
  datcont = np.ma.count(omf)

  ##Without BC
  omfnbc_density = stats.gaussian_kde(omfnbc)
  max=np.max(np.abs(omfnbc))
  omfnbc_bins = np.linspace(-_max, _max, nbins)
  print('omfnbc_bins=',omfnbc_bins)

  stdev_nbc = np.nanstd(omfnbc)  # Standard deviation
  omean_nbc = np.nanmean(omfnbc) # Mean of the data
  datmi_nbc = np.nanmin(omfnbc)  # Min of the data
  datma_nbc = np.nanmax(omfnbc)  # Max of the data
  datcont_nbc = np.ma.count(omfnbc)


  # Set plot variable unit
  # -----------------------
  units = 'K'

  fig1 = plt.figure(figsize=(8.0,7.5))
  ax=fig1.add_subplot(111)
  plt.plot(omf_bins, omf_density(omf_bins))
  plt.plot(omfnbc_bins, omfnbc_density(omfnbc_bins))
  plt.title(':omb histogram')
  #   figname='ufo_'+thisobstype+'_stage1_hist_'+subtask+'.png'
  figname=instrument + '_Ch' + str(channel) + '_histogram.png'

  # ------------------
  vtitle = instrument + " Channel " + str(channel)
  ax.set_title("OMB: "+vtitle, pad=20)

  ax.text(0.45, -0.08, 'BT O-B', transform=ax.transAxes, ha='left')
  ax.text(-0.08, 0.4, 'Normalized Count', transform=ax.transAxes,
          rotation='vertical', va='bottom')

  text = f"Total Count:{datcont:0.0f}, Max/Min/Mean/Std: {datma:0.3f}/{datmi:0.3f}/{omean:0.3f}/{stdev:0.3f} {units}"
  print(text)
  ax.text(0.67, -0.1, text, transform=ax.transAxes, va='bottom', fontsize=6.2)

  plt.tight_layout()

  # show plot
  # -----------
  plt.savefig(figname,bbox_inches='tight',dpi=100)
  return

# Set input and output file names
# --------------------------------
fname = "/scratch2/NCEPDEV/fv3-cam/Xiaoyan.Zhang/noscrub/RRFS_satbias_hold/gdas_bias/gdas.20210511/00/atmos/diag/diag_abi_g16_ges.2021051100.nc4"
pname = "gdas.abi_g16_ges.2021051100.jpg"

print(" checking fname = ", fname)

channel = 10 #to be read and plot
instrument = 'ABI_G16'
omf, omfnbc = _read_netcdf_diag(fname,channel)
print('omf=',omf)
print('omfnbc=',omfnbc)
_plot_histgram(omf,omfnbc,instrument,channel)
exit()




