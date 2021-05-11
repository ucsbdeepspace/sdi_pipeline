#!/bin/sh

#Creates images + lists for a set of simulated data
#
#Use: ./runsim.sh
#
#Needs: create_trans.py, sky_stars.conf, sky.conf, stuff.conf in working directory
#
#In order to change characteristics of images such as surface brightness,
#exposure time, etc --> you need to edit the sky.conf
#
#Arguments:
#	1st: Number of images per set (default = 10)
#
#Last Modified: 4/24/21
if [ $# -eq 0 ]
then
	images=10 #default to 10 images per set
else
	images=$1 #user-input for number of images per set
fi

#Sets up catalog for non-transient star field
stuff -c stuff.conf
sky g.list -c sky_stars.conf
sed -i '/^2/d' sky.list #Removes galaxies from catalog (Needed for a characteristic starfield)
mv sky.list background.list
mkdir -p ./simages

for count in $(seq 1 1 $images)
do

	sky background.list -c sky.conf
	echo
	mv sky.fits background$count.fits
	echo
	mv background$count.fits ./simages/

	echo "====================================="
	echo "Files for $count Created"
	echo "====================================="
done
echo "All Files Created"
mv background.list ./simages/
rm g.list i.list r.list #Unneeded byproduct of STUFF (g.list no longer needed)
rm sky.list #Leftover from SKYMAKER

