#!/bin/ksh
set -x
date=20210719
hour=19
#---------------------------------------------------------------
#Pressure level    lon     lat    i,j (ny,nx)
# 850     22       262.1   37.7   112,198
# 500     36       267.9   25.0   4,  243
# 200     45       287     47.9   219,340
# 10      59       239.8   45.0   192,64
#-------------------------------------------------------------
be_type=New
#be_type=B_Invert_Norm
#be_type=B_XInvert_Norm
#be_type=NAM
#exp=dlsig_glb_850U_vs0.3
#exp=folk_1.0_T50
exp=merge_1.0_U500
#exp=G2_NAM_T500
#ilev0=38 #850
#ilev0=40 #850 after vz*dlsig_glb
#ilev0=24 #36 #500
ilev0=26 #34 #500 vz*dlsig_glb
#ilev0=15 #200
ilev0=4 #50
#ilev0=1 #10
cross_x=198
cross_y=112
#cross_x=243
#cross_y=4
#cross_x=340
#cross_y=219
#cross_x=64
#cross_y=192
da=3DVar
#da=3DEnsVar
ANA_DIR=/scratch2/NCEPDEV/stmp3/Xiaoyan.Zhang/Exp0C775-Gens-Coldstart/tmpnwprd-2021072000/gsianl_tm05_00
GES_DIR=/scratch2/NCEPDEV/stmp3/Xiaoyan.Zhang/com/fv3cam/Exp0C775-Gens-Coldstart/fv3sar.20210720/00/guess.tm05

ln -sf $GES_DIR/${date}.${hour}0000.fv_core.res.tile1_new.nc ./bg-fv3_dynvar.nc
ln -sf $GES_DIR/${date}.${hour}0000.fv_tracer.res.tile1_new.nc ./bg-fv3_tracer.nc
ln -sf $GES_DIR/grid_spec_new.nc ./fv3_grid_spec.nc
ln -sf $ANA_DIR/fv3_dynvars_${exp} ./gsi-analysis_dynvar.nc
ln -sf $ANA_DIR/fv3_tracer_${exp} ./gsi-analysis_tracer.nc

python plot_hz_inc_all.py ${date}${hour} $da ${be_type} ${ilev0} 
mv hz_inc.png hz_inc_${da}_${exp}.png 

python plot_vc_inc_all_lon.py  ${date}${hour} $da ${be_type} ${cross_x} ${cross_y}
mv vc_inc.png vc_inc_${da}_${exp}.png 

python plot_vc_inc_all_lat.py  ${date}${hour} $da ${be_type} ${cross_x} ${cross_y} 
mv vc_inc.png vc_lat_inc_${da}_${exp}.png 

python  plot_pf_inc_all_z.py  ${date}${hour} $da ${be_type} ${cross_x} ${cross_y}
mv vc_inc.png pf_z_inc_${da}_${exp}.png 

python  plot_pf_inc_all_lon.py  ${date}${hour} $da ${be_type} ${cross_x} ${cross_y} ${ilev0} 
mv vc_inc.png pf_lon_inc_${da}_${exp}.png 
