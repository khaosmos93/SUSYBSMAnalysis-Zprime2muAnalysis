#!/bin/bash

index=1
while [ ${index} -le 6 ]; do

  echo "`date`";

  python crabStatusResubmit.py resubmit >&resubmit_All.log&

  index=$((index+1))
  sleep 3600;
done;
