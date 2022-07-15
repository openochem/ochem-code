#!/bin/sh

# stop the screen and all running servers insinde
echo "Stopping cs_servers screen"
kill -9 `ps aux | grep cs_servers | grep java | awk -F ' ' {'print $2'}`  
kill -9 `ps aux | grep serv.ServerRunner | grep java | awk -F ' ' {'print $2'}`
kill -9 `ps aux | grep serv.MultiServer | grep java | awk -F ' ' {'print $2'}`
kill -9 `ps aux | grep CalculationServerScreen | grep cs_screenrc | awk -F ' ' {'print $2'}`
screen -wipe
