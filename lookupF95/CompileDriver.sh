#! /bin/bash

\rm -f SurveySimulator
cd ../srcF95
make Driver
cd ../lookupF95
ln -s ../srcF95/Driver SurveySimulator

exit
