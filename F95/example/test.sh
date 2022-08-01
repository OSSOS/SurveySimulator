#!/bin/bash
cmd=$1
language=$2
src_dir=$(dirname ${1})
cmd=$(basename ${1})
echo "Testing ${cmd}"

if [ -f SimulDetect.dat ]; then
    \rm SimulDetect.dat
fi
if [ -f SimulTrack.dat ]; then
    \rm SimulTrack.dat
fi

${src_dir}/Driver < ${cmd}.in
if [ $? != 0 ]; then
  status=$?
  echo "================================================================="
  echo "=======================  Test FAILED!   ========================="
  echo "================================================================="
  exit ${status}
fi

first_detect_line=$(grep -v "^#" SimulDetect.dat | head -1)
a=`echo ${first_detect_line} | awk '{print $1}'`
s=`echo ${first_detect_line} | awk '{print $7}'`
no=`tail -3 SimulDetect.dat | head -1 | awk '{printf "%10d", $6}'`
nd=`tail -2 SimulDetect.dat | head -1 | awk '{print $5}'`
nt=`tail -1 SimulDetect.dat | awk '{print $6}'`
first_check_line=$(grep -v "^#" ${cmd}-check-${language}.dat | head -1)
ac=`echo ${first_check_line} | awk '{print $1}'`
sc=`echo ${first_check_line} | awk '{print $15}'`
ndo=`tail -3 ${cmd}-check-${language}.dat | head -1 | awk '{printf "%10d", $6}'`
ndc=`tail -2 ${cmd}-check-${language}.dat | head -1 | awk '{print $5}'`
ntc=`tail -1 ${cmd}-check-${language}.dat | awk '{print $6}'`

if [ $a != $ac -o $s != $sc ]; then
    echo "Error: first detection doesn't match!"
    echo "$a != $ac, $s != $sc"
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
echo "================================================================="
echo "=====================  Test successful!  ========================"
echo "================================================================="

exit
