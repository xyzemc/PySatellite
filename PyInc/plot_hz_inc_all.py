import matplotlib
from matplotlib import colors
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.tri as tri
#from mpl_toolkits.basemap import Basemap
import numpy as np
import datetime, os, sys, subprocess
from netCDF4 import Dataset


ymdh = str(sys.argv[1])
da = str(sys.argv[2])
be_type = str(sys.argv[3])
ilev0 = str(sys.argv[4])
print(ilev0)

ymd = ymdh[0:8]
year = int(ymdh[0:4])
month = int(ymdh[4:6])
day = int(ymdh[6:8])
hour = int(ymdh[8:10])
print (year, month, day, hour)

nesteddata = './bg-fv3_dynvar.nc'
nesteddata1 = './gsi-analysis_dynvar.nc'
nestedgrid = './fv3_grid_spec.nc'

tracerdata = './bg-fv3_tracer.nc'
tracerdata1 = './gsi-analysis_tracer.nc'

fnd = Dataset(nesteddata,'r')
fnd1 = Dataset(nesteddata1,'r')

tnd = Dataset(tracerdata,'r')
tnd1 = Dataset(tracerdata1,'r')

print( 'data dimension is ',fnd.dimensions.keys())
print( 'varialbe is ',fnd.variables.keys())
print ('varialbe is ',fnd.variables.keys())
#print 'varialbe is ',type(fnd.variablesi['T'])
print ('varialbe is ',type(fnd.variables['T'][0,10,:,:]))
tvar = fnd.variables['delp'][0,0,:,:]
fng = Dataset(nestedgrid,'r')
#ilev0=27
#ilev0=32
#ilev0=25
#ilev0=60 #10mb
fhours = [1]
#fhours = [1,2,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60]
dtime = datetime.datetime(year,month,day,hour,0)
date_list = [dtime + datetime.timedelta(hours=x) for x in fhours]
print (date_list)
for j in range(len(date_list)):

  fhour = str(fhours[j]).zfill(2)
  fhr = int(fhour) - 1
  print ('fhour '+fhour)
## Temperature
  t0 = fnd.variables['T'][fhr,ilev0,:,:]
  t1 = fnd1.variables['T'][fhr,ilev0,:,:]

## U-wind
  u0 = fnd.variables['u'][fhr,ilev0,:,:]
  u1 = fnd1.variables['u'][fhr,ilev0,:,:]

## V-wind
  v0 = fnd.variables['v'][fhr,ilev0,:,:]
  v1 = fnd1.variables['v'][fhr,ilev0,:,:]

## Q-wind
  q0 = tnd.variables['sphum'][fhr,ilev0,:,:]
  q1 = tnd1.variables['sphum'][fhr,ilev0,:,:]

  es = 360 #   6.11*np.exp((53.49*np.ones(t850.shape))-((6808.*np.ones(t850.shape))/t850)-5.09*np.log(t850))
  ws = 0 #$ 0.622*(es/850.)
#  thetae = (t850+((2450000/1005.7)*((rh850/100)*ws))) * (1000/850)**(0.2854)
  units = 'K'

  lat = fng.variables['grid_latt'][:,:]
  lon = fng.variables['grid_lont'][:,:]
  print( "dims are ",lat.shape) 
  print( "u dims are ",t0.shape) 
  nx=lat.shape[1]
  ny=lat.shape[0]
  print ("nx, ny are ",nx,ny)
  latplt=np.transpose(lat,(1,0))
  lonplt=np.transpose(lon,(1,0))
##Increment
  dt=(t1-t0)[0:ny,0:nx]
  du=(u1-u0)[0:ny,0:nx]
  dv=(v1-v0)[0:ny,0:nx]
  dq=(q1-q0)[0:ny,0:nx]

  print ("maximum location is ",np.unravel_index(dt.argmax(),dt.shape) ) 
  print ("minimum location is ",dt.argmin()) 
  
  t0plt=np.transpose(dt,(1,0))
  print ("dims for latplt ",latplt.shape)
  print ("dims for t0plt ",t0plt.shape)
  nlon=latplt.shape[0]
  nlat=latplt.shape[1]
  
  print ("nlat and nlon are ",nlat,nlon)
  
  gmax=np.max(lat)
  gmin=np.min(lat)
  print ("xxdeb",gmin,gmax,nlat)
  lat1d=np.linspace(gmin,gmax,nlat)
  gmax=np.max(lon)
  gmin=np.min(lon)
  lon1d=np.linspace(gmin,gmax,nlon)

# Create the figure
  #fig,ax = plt.subplots(2,2)
  #plt.axis([250,280,30,46])
  figsize=(6,4)
  ax=plt.figure(figsize=figsize, constrained_layout=True).subplots(2, 2)
  plt.suptitle('{} BE Horizontal Increment'.format(be_type))

#1 T Increment
  ax[0,0].set_title('GSI {} T inc'.format(da),fontsize=8)
  #ax1plt=ax[0,0].contour(lon,lat, dt,cmap='RdGy')
  #ax1plt=ax[0,0].contour(lon,lat, dt,3,colors='black')
  print('dt.max=',dt.max())
  clevs = np.arange(0,dt.max(),0.2)
  print('clevs=',clevs)
  #ax1plt=ax[0,0].contour(lon,lat, dt,clevs,colors='black')
  ax1plt=ax[0,0].contour(lon,lat, dt, colors='black')
  ax[0,0].clabel(ax1plt, inline=1, fontsize=8)
  #ax[0,0].axis([246,280,30,48])
  #ax[0,0].axis([236,289,22,50])
  ax[0,0].set_ylabel('Latitude')
#  fig.colorbar(ax1plt,ax=ax)

#2 Q Increment
  print('dq.max=',dq.max())
  ax[0,1].set_title('GSI {} Q inc'.format(da),fontsize=8)
  #ax[0,1].axis([246,280,30,48])
  #ax[0,0].axis([236,289,22,50])
  #ax1plt=ax[0,1].contour(lon,lat, dq,cmap='RdGy')
  ax1plt=ax[0,1].contour(lon,lat, dq,10,colors='black')
  ax[0,1].clabel(ax1plt, inline=True, fontsize=8)

#3 U Increment
  print('du.max=',du.max())
  ax[1,0].set_title('GSI {} U inc'.format(da),fontsize=8)
  #ax[1,0].axis([246,280,30,48])
  #ax[0,0].axis([236,289,22,50])
  #ax1plt=ax[1,0].contour(lon,lat, du,cmap='RdGy')
  ax1plt=ax[1,0].contour(lon,lat, du,10,colors='black')
  ax[1,0].clabel(ax1plt, inline=True, fontsize=8)
  ax[1,0].set_ylabel('Latitude')
  ax[1,0].set_xlabel('Longitude')

#4 V Increment
  print('dv.max=',dv.max())
  ax[1,1].set_title('GSI {} V inc'.format(da),fontsize=8)
  #ax[1,1].axis([246,280,30,48])
  #ax[0,0].axis([236,289,22,50])
  #ax1plt=ax[1,1].contour(lon,lat, dv,cmap='RdGy')
  ax1plt=ax[1,1].contour(lon,lat, dv,10,colors='black')
  ax[1,1].clabel(ax1plt, inline=True, fontsize=8)
  ax[1,1].set_xlabel('Longitude')

  plt.savefig('hz_inc.png', bbox_inches='tight',dpi=150)

#  plt.show()
  plt.close()
