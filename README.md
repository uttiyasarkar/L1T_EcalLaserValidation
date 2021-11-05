# L1T_EcalLaserValidation
CMS Ecal laser calibrations fast validation at Level-1

## Introduction

The procedures are wrapped within the L1RateValidation.sh script. 
It use the Zerobias from Run 320065 Fill#6961. 

The script will check out the l1t-integration-v104.1 L1Ntuple code
and reemulate different L1Ntuple with different Ecal laser configs.

It runs on CMSSW 11_0_2 base presently.


##Running instructions for local checks

Go to the directory ToRun.

Copy a file NewToRun.txt from the directory RunFiles here: cp ../RunFiles/NewToRun.txt .

A sample run file will look something like this

week 48
year 2018
run1 327385
run2 327466
type laser

Change the numbers of week, year, run1, run2 and running type here.

Do cd ..

Now launch the script by new.sh like this:

./new.sh

and voila!





