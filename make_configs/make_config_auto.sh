#!/bin/bash
#Check if we already have a RSPTmake.inc
if [ -f ../make_config ]
then
   exit
fi

#host
myhost=`hostname |awk -F"[0-9.]" '{print $1}'`

#cluster (Add your favourite cluster in the list below!)
cluster="default"
if [[ "$OSTYPE" =~ "darwin" ]]
then
   # Catch Macs
   cluster="mac"
elif [ "$myhost" == "n" ]
then
   cluster="hclm"
elif [ "$myhost" == "vsc3" ]
then
   cluster="vsc3"
else
   # Fall back
   cluster="$myhost"
fi

# Check if we already have a working template
template=""
if [ -n "$cluster" ] 
then
   template=`ls -1 make_config_*|grep -i "make_config_$cluster" | head -1`
fi

# Generate the RSPTmake.inc file
if [ -n "$template" ]
then
   cp $template ../make_config
else
   echo "make_config_auto.sh failed to find a suitable make_config. Please look at the available templates in make_configs."
fi

