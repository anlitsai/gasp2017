a=`cat check_science_target_list.txt`
for i in $a;do python Rmag_aperture_annulus_r_file_median_w0_subplot_target_2017.py $i | tee python Rmag_aperture_annulus_r_file_median_w0_subplot_$i.log;done
