#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
ARCH=slc6_amd64_gcc630
CMSREL=CMSSW_10_0_0
L1TTag=l1t-integration-v97.17-v2
GT=100X_dataRun2_v1
sqlite1=$1 ##ref
sqlite2=$2
week=$3
year=$4
curdir=$PWD
username=$USER
pids=""
hasref=false

file=/store/data/Run2017F/ZeroBias4/RAW/v1/000/306/091/00000/00502E87-FBBF-E711-A0E8-02163E01A2A6.root
xrdcp -f root://cms-xrd-global.cern.ch/$file /tmp/$username.root

#----------------------------------------------------------------------------#
#                            Getting the reference                           #
#----------------------------------------------------------------------------#
#it may be gzipped infact ...
## prevent exit from failed wget
set +e 
wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/${1}/L1Menu_${GT}_${sqlite1}_emu.csv \
  https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/${1}/L1Menu_${GT}_${sqlite1}_emu.csv
if [ $? -ne 0 ]; then
  sqs="$sqlite1 $sqlite2"
else
  sqs=$sqlite2
  hasref=true
fi
echo "Running ECal validtion with ", $sqs

#----------------------------------------------------------------------------#
#                            Checkout L1 Emulator                            #
#----------------------------------------------------------------------------#
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=$ARCH
scramv1 project CMSSW $CMSREL
cd $CMSREL/src
eval `scramv1 runtime -sh`
git-cms-init
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:$L1TTag
git cms-addpkg L1Trigger/L1TCommon
git cms-addpkg L1Trigger/L1TMuon
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TMuon.git L1Trigger/L1TMuon/data

scram b -j $((`nproc`/2))

#----------------------------------------------------------------------------#
#                          Running the TP emulation                          #
#----------------------------------------------------------------------------#
echo running $GT
cmsDriver.py l1Ntuple -s RAW2DIGI --era=Run2_2017  \
  --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU \
  --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWCalouGT \
  --conditions=$GT -n -1 --data --no_exec --no_output  \
  --filein=file:/tmp/$username.root \
  --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloStage2Params_2017_v1_8_2_updateHFSF_v6MET \
  --python_filename=l1Ntuple_${GT}.py

for sq in $sqs; do
  if [ ! -f EcalTPG_${sq}_moved_to_1.db ]; then
    wget http://cern.ch/ecaltrg/EcalLin/EcalTPG_${sq}_moved_to_1.db
  fi
  python ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
  cmsRun l1Ntuple_${GT}_${sq}.py >& l1Ntuple_${GT}_${sq}.log 
  pids="$pids $!"
done
################################


#----------------------------------------------------------------------------#
#                            Check out L1Menu code                           #
#----------------------------------------------------------------------------#
git clone -b EcalVal_940 https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu/macros
cp $curdir/CompL1Rate.py  .
cp $curdir/menulib.cc .
cp $curdir/menulib.hh .
cp $curdir/Prescale_Sets_RUN_306091_col_1.5.txt menu/
cp $curdir/Selected_Seed.txt menu/
make -j $((`nproc`/2))

#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#
echo "Waiting for Ntuple production to finish......"
wait $pids

for sq in $sqs; do
  ./testMenu2016 -m menu/Prescale_Sets_RUN_306091_col_1.5.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.root -o L1Menu_${GT}_${sq}_emu >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
  ./testMenu2016 -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.root -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
  pids="$pids $!"
done
echo "Waiting for menu rate estimation to finish......"
wait $pids

#----------------------------------------------------------------------------#
#                                Compare rate                                #
#----------------------------------------------------------------------------#

if $hasref; then
  mv $curdir/L1Menu_${GT}_${sqlite1}_emu.csv results/
  mv $curdir/L1Seed_${GT}_${sqlite1}_emu.csv results/
fi

python CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log

#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#
if [ -f ${WORKSPACE}/upload/$2 ]
then
  echo "dir is already existing"
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
else
  mkdir ${WORKSPACE}/upload/$2
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
fi

#we need to make a tar gz of this one
cp results/L1Menu_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp results/L1Seed_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp ${sqlite2}.log ${WORKSPACE}/upload/${2}/
