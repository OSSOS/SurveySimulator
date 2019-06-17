#! /bin/bash

# Determine directory name
ID=`date +"%Y%m%d"`
d="../SurveySim${ID}"
da=`date +"%F"`
dt=`LC_ALL=us date +"%d %b %Y"`
df=`LC_ALL=us date +"%b %d/%Y"`
v="2.0"
intended_audience="OSSOS team"
#intended_audience="General Public"

# Create a clean distribution directory
if [ -d $d ]; then
    \rm -rf $d
fi
mkdir $d
mkdir $d/src
mkdir $d/srcF95
if [ -e ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
if [ -h ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
ln -s ${d/./} ../CurrentDistrib

# Copy fix files
cat > $d/README.version <<EOF

Survey Simulator for OSSOSv11

Survey simulator as of $dt

EOF
head -2 README.first > $d/README.first
cat >> $d/README.first <<EOF
$df release to $intended_audience
EOF
tail --line=+3 README.first >> $d/README.first
cp -a eupl* README.contact lookup lookupF95 parametric Python $d/
for s in CFEPS OSSOS All_r_Surveys All_Surveys Deep_Surveys; do
    mkdir $d/$s
    cp ../$s/* $d/$s/
    \rm $d/$s/*pointings
done
for s in CFEPS OSSOS All_r_Surveys Deep_Surveys; do
    \rm -f $d/$s/README.formats
    ln -s ../All_Surveys/README.formats $d/$s/README.formats
done
cp src/Driver.{f,py} src/README.* src/ModelUtils.f $d/src/
cp srcF95/*.f95 srcF95/Makefile $d/srcF95/
cp src/README.* $d/srcF95/
\rm -f $d/CFEPS/*.py
\rm -f $d/srcF95/test*f95

# Initialize SurveySubs.f
cd src
cp SurveySubs.f zzzz0
n=`grep -i include zzzz0 | grep -i -v "[a-z].*include" | grep -i -v "include 'param.inc'" | wc -l`

# Now loop on including "include" files, except "include 'param.inc'"
while [ $n -gt 0 ]; do
    cp zzzz0 zzzz1
    grep -i include zzzz0 | grep -i -v "[a-z].*include" | grep -i -v "include 'param.inc'" | (
	while read inc file rest; do
	    grep -v "include.*$file" zzzz1 | grep -v "INCLUDE.*$file" > zzzz2
	    cat `echo $file | sed "s/'//g"` >> zzzz2
	    \mv zzzz2 zzzz1
	done
    )
    \mv zzzz1 zzzz0
    n=`grep -i include zzzz0 | grep -i -v "[a-z].*include" | grep -i -v "include 'param.inc'" | wc -l`
done

# Now inline the "include 'param.inc'" statements.
cp zzzz0 zzzz1
../InlineIncludeParam.py
\rm -f zzzz1

# Move the result in the distribution directory
cp SurveySubsHistory ../$d/src/SurveySubs.f
cat >> ../$d/src/SurveySubs.f <<EOF
c
c File generated on $da
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

EOF
cat zzzz0 >> ../$d/src/SurveySubs.f
\rm zzzz0
cd ../

# Create the tarball
tar czf ../SurveySimulator-${ID}.tgz $d

exit
