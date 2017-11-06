#!/usr/bin/env bash
source ../../runcommand
# Determine if we have a parallel build
fstring1=`grep "\-DMPI" ../../../make_config|grep -v "^#" |grep -v "[: :]#"`
fstring2=`echo $fstring1 | sed "s/#.*//"`
if `echo ${fstring2} | sed "s/#.*//" | grep "\-DMPI" 1>/dev/null 2>&1`
then
    RUNCOMMAND=`echo $MPIRUN ../../../bin/$EXEC config_test`
else
    RUNCOMMAND=`echo ../../../bin/$EXEC config_test`
fi
echo $RUNCOMMAND
# Run the command
$RUNCOMMAND
