#!/bin/ksh

### used in output_uhuru.job to find the maximum-minimum between elevation, crustal and lithospheric thickness.

read fileres_t	    ## e_sedsL_Tm_vis_Q_epeff_epzz.xy from graficsth
read regio_vel
read x1
read y1

echo $fileres_t
echo $regio_vel
echo $x1  $y1
Gaussian=80
yrestrict=1000	## [km]
echo " ¡¡¡¡¡¡¡¡¡¡  TAKE CARE, THERE IS A RESTRICTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "     Just consider the domain where y > " $yrestrict " km"

#fileres_t=$file_graf$step1
awk '{if(NR==1){print $4 }}' $fileres_t | read Time 		
awk '{if($1!="#"){print $0 }}' $fileres_t > file_xy.tmp

##### ELEVATION
awk '{if($1!="#"){print $1,$2,$3 }}' file_xy.tmp > file_xy_e.tmp	# elevation
   #xyz2grd file_xy_e.tmp -R$regio_vel -Gfile_grd.tmp -I$x1/$y1
   #grdfilter file_grd1.tmp -D0 -Fg$Gaussian -Gfile_grd.tmp -V
   #grd2xyz file_grd.tmp > file_xy_e.tmp
awk '{if($2>yrestrict){printf " %12.3f %12.3f %10.4f\n ", $1,$2,$3 }}' yrestrict=$yrestrict file_xy_e.tmp > file_xy_ff.tmp	
wc file_xy_ff.tmp | read nrow a b c
sort -k 3,3nr -o file_ord_e.tmp file_xy_ff.tmp		
awk '{if(NR==1){print $1,$2,$3 }}' file_ord_e.tmp | read x_emax y_emax emax
echo $Time $emax $x_emax $y_emax >> file_time_emax.tmp
awk '{if(NR==nrow){print $1,$2,$3 }}' nrow=$nrow file_ord_e.tmp | read x_emin y_emin emin
echo $Time $emin $x_emin $y_emin >> file_time_emin.tmp

##### CRUST
awk '{if($1!="#"){print $1,$2,$5 }}' file_xy.tmp > file_xy_s.tmp	# crustal thickness
   #xyz2grd file_xy_s.tmp -R$regio_vel -Gfile_grd.tmp -I$x1/$y1
   #grdfilter file_grd1.tmp -D0 -Fg$Gaussian -Gfile_grd.tmp -V
   #grd2xyz file_grd.tmp > file_xy_s.tmp
awk '{if($2>yrestrict){printf " %12.3f %12.3f %10.4f\n ", $1,$2,$3 }}' yrestrict=$yrestrict file_xy_s.tmp > file_xy_ff.tmp	
wc file_xy_ff.tmp | read nrow a b c
sort -k 3,3nr -o file_ord_s.tmp file_xy_ff.tmp		
awk '{if(NR==1){print $1,$2,$3 }}' file_ord_s.tmp | read x_smax y_smax smax
echo $Time $smax $x_smax $y_smax >> file_time_smax.tmp
awk '{if(NR==nrow){print $1,$2,$3 }}' nrow=$nrow file_ord_s.tmp | read x_smin y_smin smin
echo $Time $smin $x_smin $y_smin >> file_time_smin.tmp

##### LITHOPSHERIC MANTLE
awk '{print $1,$2,$6-$4-$5}' file_xy.tmp > file_xy_lm.tmp   		# lithospheric mantle thickness
   #xyz2grd file_xy_lm.tmp -R$regio_vel -Gfile_grd.tmp -I$x1/$y1
   #grdfilter file_grd1.tmp -D0 -Fg$Gaussian -Gfile_grd.tmp -V
   #grd2xyz file_grd.tmp > file_xy_lm.tmp
awk '{if($2>yrestrict){printf " %12.3f %12.3f %10.4f\n ", $1,$2,$3 }}' yrestrict=$yrestrict  file_xy_lm.tmp > file_xy_ff.tmp	
wc file_xy_ff.tmp | read nrow a b c
sort -k 3,3nr -o file_ord_lm.tmp file_xy_ff.tmp		
awk '{if(NR==1){print $1,$2,$3 }}' file_ord_lm.tmp | read x_lmmax y_lmmax lmmax
echo $Time $lmmax $x_lmmax $y_lmmax >> file_time_lmmax.tmp
awk '{if(NR==nrow){print $1,$2,$3 }}' nrow=$nrow file_ord_lm.tmp | read x_lmmin y_lmmin lmmin
echo $Time $lmmin $x_lmmin $y_lmmin >> file_time_lmmin.tmp

#### VERTICAL STRAIN RATE	[$11*1e16]
awk '{if($1!="#"){print $1,$2,$11*1e16 }}' file_xy.tmp > file_xy_ezz.tmp	# vertical strain rate
awk '{if($2>yrestrict){printf " %12.3f %12.3f %10.4f\n ", $1,$2,$3 }}' yrestrict=$yrestrict file_xy_ezz.tmp > file_xy_ff.tmp	
wc file_xy_ff.tmp | read nrow a b c
sort -k 3,3nr -o file_ord.tmp file_xy_ff.tmp		
awk '{if(NR==1){print $1,$2,$3 }}' file_ord.tmp | read x_max y_max max
echo $Time $max $x_max $y_max >> file_time_ezzmax.tmp
awk '{if(NR==nrow){print $1,$2,$3 }}' nrow=$nrow file_ord.tmp | read x_min y_min min
echo $Time $min $x_min $y_min >> file_time_ezzmin.tmp





