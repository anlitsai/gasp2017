a=`find ./|grep slt201709[0-9][0-9]|cut -d / -f3| sort|uniq`
for i in $a;do python slt_calibration_science_calibration_180S.py $i;done
