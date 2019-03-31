#! /bin/bash

\rm -f SurveySimulator
cd ../srcF95
\rm -f gimeobjut.f95
ln -s ../lookupF95/ReadModelFromFile.f95 gimeobjut.f95
make Driver
cd ../lookupF95
ln -s ../srcF95/Driver SurveySimulator

exit
