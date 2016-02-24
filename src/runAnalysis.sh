#!/bin/bash
# ok, it's working
#launch as --> source <script_name>.sh
cmsenv
echo
echo "Compiling code with make command"
echo
make
echo
MACRO_FOLDER="macro"
MACRO_NAME="distribution.C"
EXE_NAME="distribution"
echo "Compiling $MACRO_FOLDER/$MACRO_NAME"
g++ -Wall -pedantic -lm -o $MACRO_FOLDER/$EXE_NAME $MACRO_FOLDER/$MACRO_NAME `rootlib`
echo
echo "Calling myjob"
bsub -q cmsan -Is /bin/bash < AnaLaunch.sh

