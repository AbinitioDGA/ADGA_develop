#!/usr/bin/env bash
source ../../runcommand
# Determine if we have a parallel build
fstring1=`grep "\-DMPI" ../../../Makefile|grep -v "^#" |grep -v "[: :]#"`
fstring2=`echo $fstring1 | sed "s/#.*//"`
if `echo ${fstring2} | sed "s/#.*//" | grep "\-DMPI" 1>/dev/null 2>&1`
then
    RUNCOMMAND=`echo $MPIRUN ../../../bin/$EXEC test.conf`
else
    RUNCOMMAND=`echo ../../../bin/$EXEC test.conf`
fi
echo $RUNCOMMAND
# Run the command
$RUNCOMMAND
