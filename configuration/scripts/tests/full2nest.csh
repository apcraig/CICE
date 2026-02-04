#!/bin/tcsh
# you need cdo and nco to run this script
# note the Input section below
###### steps to be taken to test full and nest #####
#1) create full set-up
#       ./cice.setup -c {NAME} -g gbox80 -m {ICE_MACHINE} -e {ICE_MACHINE_ENVNAME} -s bdyrestorfull
#2) build full set-up
#       cd {NAME}; 
#       ./cice.build
#3) run full set-up
#       ./cice.submit
#4) create nest set-up
#       ./cice.setup -c {NAME} -g gbox12 -m {ICE_MACHINE} -e {ICE_MACHINE_ENVNAME} -s bdyrestornest
#5) build nest set-up
#       cd {NAME};
#       ./cice.build
#6) use this file to creat all needed fields for running the nest
#       cp full2nest.sh /run/dir/nest/
#7) modefy {PATH2hist} and {PATH2rst} to point to the history and restart directory
#8) run full2nest.sh
#       ./full2nest.sh
#9) run nest set-up
#       cd /nest/home/dir
#       ./cice.submit
###############END of description #################

# to be executed in the rundirectory of the nest
#------------------input--------------------------

set CASE = "bdyf11"
set PATH2hist = "/glade/derecho/scratch/tcraig/CICE_RUNS/${CASE}/history/"
set PATH2rst = "/glade/derecho/scratch/tcraig/CICE_RUNS/${CASE}/restart/"
set ninest=12
set njnest=12
set nis=40
set njs=40
#----------------end input------------------------
@ nie = $nis + $ninest + 1
@ nje = $njs + $njnest + 1
echo "cutting $nis,$nie $njs,$nje"
#
# make the kmt_grid from full domain 
# most variables can be cut from full domain
# kmt needs to be created

# cut domain from history file of full domain: 
ncks -d ni,$nis,$nie -d nj,$njs,$nje ${PATH2hist}/iceh.2005-01-02.nc temp.nc

#convert sigP to kmt and cut relevant variables:
ncrename -v tmask,kmt temp.nc
ncks -C -v kmt,HTN,HTE,NLON,NLAT,ANGLET,ANGLE,TLAT,TLON,ULAT,ULON temp.nc temp1.nc
ncap2 -s "HTE=HTE*100.;HTN=HTN*100.;ULON=ULON/(180/3.14159);ULAT=ULAT/(180/3.14159);TLON=TLON/(180/3.14159);TLAT=TLAT/(180/3.14159)" temp1.nc grid_kmt.nc

ncrename -v HTE,hte -v HTN,htn -v ULON,ulon -v ULAT,ulat -v ANGLE,angle -v ANGLET,anglet -v TLAT,tlat -v TLON,tlon grid_kmt.nc
rm temp*.nc

# cut boundary files from full domain
# one bigger than grid, since restart ext needs to be true

foreach dd (01 02 03 04 05 06)
foreach ss (00000 03600 07200 10800 14400 18000 21600 25200 28800 32400 36000 39600 43200 46800 50400 54000 57600 61200 64800 68400 72000 75600 79200 82800)
    ncks --no-abc -d ni,$nis,$nie -d nj,$njs,$nje ${PATH2rst}/iced.2005-01-${dd}-${ss}.nc cice_bdy_restart200501${dd}${ss}.nc
end
end


