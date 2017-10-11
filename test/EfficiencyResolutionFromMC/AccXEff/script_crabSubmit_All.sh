#!/bin/bash

python ./noWei_noRand/crabSubmit.py submit >&./noWei_noRand/crabSubmit.log&
python ./yesWei_noRand/crabSubmit.py submit >&./yesWei_noRand/crabSubmit.log&
python ./yesWei_yesRand/crabSubmit.py submit >&./yesWei_yesRand/crabSubmit.log&
