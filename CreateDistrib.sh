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
if [ -e ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
if [ -h ../CurrentDistrib ]; then
    \rm ../CurrentDistrib
fi
ln -s ${d/./} ../CurrentDistrib
curdir=$(pwd)

# Copy fix files
cat > $d/README.version <<EOF

## Survey Simulator for OSSOSv12

## Survey simulator as of $dt

EOF
head -2 README.md > $d/README.md
cat >> $d/README.md <<EOF
$df release to $intended_audience
EOF
tail --line=+3 README.md >> $d/README.md
cp -a eupl* Simulator Models $d/
for s in CFEPS OSSOS All_r_Surveys All_Surveys; do
    mkdir -p $d/Surveys/$s
    cp Surveys/$s/* $d/Surveys/$s/
done
for s in CFEPS OSSOS All_r_Surveys; do
    \rm -f $d/Surveys/$s/README.formats
    ln -s ../All_Surveys/README.formats $d/Surveys/$s/README.formats
done
\rm -f $d/Simulator/F77/fortran/Test*f
\rm -f $d/Simulator/F95/fortran/test*f95

# Initialize SurveySubs.f
cd $d/Simulator/F77/fortran
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
${curdir}/InlineIncludeParam.py
\rm -f zzzz1

# Move the result in the distribution directory
cat > SurveySubs.f <<EOF
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c
c File generated on $da
c
c-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

EOF
cat zzzz0 >> SurveySubs.f
\rm GetSurvey.f EffUtils.f PosVelUtils.f Polygon-lib.f ElemPosUtils.f Rotation.f zzzz0
cd ${curdir}

# Create the tarball
tar czf ../SurveySimulator-${ID}.tgz $d

exit
