#! /bin/bash
cmd=$(basename ${1})
echo "Testing ${cmd}"

if [ -f SimulDetect.dat ]; then
    \rm SimulDetect.dat
fi
if [ -f SimulTrack.dat ]; then
    \rm SimulTrack.dat
fi

time LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../fortran ./Driver.py < ${cmd}.in > LOG

a=`head -10 SimulDetect.dat | tail -1 | awk '{print $1}'`
s=`head -10 SimulDetect.dat | tail -1 | awk '{print $15}'`
no=`tail -3 SimulDetect.dat | head -1 | awk '{printf "%10d", $6}'`
nd=`tail -2 SimulDetect.dat | head -1 | awk '{print $5}'`
nt=`tail -1 SimulDetect.dat | awk '{print $6}'`
ac=`head -10 ${cmd}-check.dat | tail -1 | awk '{print $1}'`
sc=`head -10 ${cmd}-check.dat | tail -1 | awk '{print $15}'`
ndo=`tail -3 ${cmd}-check.dat | head -1 | awk '{printf "%10d", $6}'`
ndc=`tail -2 ${cmd}-check.dat | head -1 | awk '{print $5}'`
ntc=`tail -1 ${cmd}-check.dat | awk '{print $6}'`

cd ../

if [ $a != $ac -o $s != $sc ]; then
    echo "Error: first detection doesn't match!"
    exit 1
fi
if [ $no != $ndo ]; then
    echo "Error: wrong number of objects!"
    exit 1
fi
if [ $nd -gt $(($ndc + 10)) -o $nd -lt $(($ndc - 10)) ]; then
    echo "Error: wrong number of detections!"
    exit 1
fi
if [ $nt -gt $(($ntc + 10)) -o $nt -lt $(($ntc - 10)) ]; then
    echo "Error: wrong number of tracking!"
    exit 1
fi
echo "Test successful!"

exit
