#!/bin/sh

# name of the screenrc that will be used
RC="cs_screenrc"

echo "Writing ${RC}..."
# add first screen as maintenance screen
echo "screen -t maintenance" > ${RC}

# eliminate old logs
rm -r server*/runs/*
rm server*/output/*
rm server*/*.zip
rm server*/tests.txt

# add one screen for all folders in the form server# where # denotes the server and screen number
for s in `ls -d server*`; do
    echo "unalias -a";
    echo "chdir $s";
    java=`perl -nle 'print $1 if m|<java-home>(.+)</java-home>|' $s/version.xml`
    echo "screen -t $s $java/bin/java -XX:MinHeapFreeRatio=20 -XX:MaxHeapFreeRatio=40 -cp lib/* qspr.metaserver.serv.ServerRunner";
    echo "chdir ..";
    echo "sleep 1";
done >> ${RC}

# select first window in the end
echo "select 0" >> ${RC}

echo "Starting cs_servers screen..."
# final start the screen with this screenrc
screen -d -m -c ${RC} -S CalculationServerScreen
