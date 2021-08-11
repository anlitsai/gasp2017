a=`find ./|grep 201710|grep GASP|grep fts|cut -d / -f3|uniq`
for i in $a;do python slt_calibration_science_calibration_180S.py $i;done

#a=`find ./|grep 201711|grep GASP|grep fts|cut -d / -f3|uniq`
#for i in $a;do python slt_calibration_science_calibration_180S.py $i;done
#a=`find ./|grep 201712|grep GASP|grep fts|cut -d / -f3|uniq`
#for i in $a;do python slt_calibration_science_calibration_180S.py $i;done
