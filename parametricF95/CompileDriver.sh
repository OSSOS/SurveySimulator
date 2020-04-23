#! /bin/bash

\rm -f SurveySimulator
cd ../srcF95
\rm -f gimeobjut.f95
ln -s ../parametricF95/InnerBeltModel.f95 gimeobjut.f95
make clean
make Driver
\rm -f gimeobjut.f95
cd ../parametricF95
cp ../srcF95/Driver SurveySimulator

exit
