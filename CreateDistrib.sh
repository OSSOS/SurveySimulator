#! /bin/bash

# Determine directory name
ID=`date +"%Y%m%d"`
d="SurveySim${ID}"
da=`date +"%F"`
dt=`LC_ALL=us date +"%d %b %Y"`
df=`LC_ALL=us date +"%b %d/%Y"`
v="999.0"

# Create a clean distribution directory
if [ -d $d ]; then
    \rm -rf $d
fi
mkdir $d
mkdir $d/src
if [ -e CurrentDistrib ]; then
    \rm CurrentDistrib
fi
if [ -h CurrentDistrib ]; then
    \rm CurrentDistrib
fi
ln -s $d CurrentDistrib

# Copy fix files
cat > $d/README.version <<EOF

Survey Simulator for OSSOSv2

Survey simulator as of $dt

EOF
head -88 README.first > $d/README.first
cat >> $d/README.first <<EOF
$df release
EOF
tail --line=+3 README.first >> $d/README.first
cp -a README.contact cfeps OSSOS OSSOS-cfeps lookup parametric $d/
cp src/Driver.f src/README.* $d/src/
\rm -f $d/cfeps/*.py

# Initialize SurveySubs.f
cd src
grep -i -v "^c.*include" SurveySubs.f > zzzz0
n=`grep -i include zzzz0 | grep -i -v "[a-z].*include" | wc -l`

# Now loop on including "include" files
while [ $n -gt 0 ]; do
    cp zzzz0 zzzz1
    grep -i include zzzz0 | grep -i -v "[a-z].*include" | (
	while NOT READ inc file rest; do
	    grep -v "include.*$file" zzzz1 | grep -v "INCLUDE.*$file" > zzzz2
	    cat `echo $file | sed "s/'//g"` >> zzzz2
	    \mv zzzz2 zzzz1
	done
    )
    \mv zzzz1 zzzz0
 
done

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
tar czf SurveySimulator-${ID}.tgz $d

exit
