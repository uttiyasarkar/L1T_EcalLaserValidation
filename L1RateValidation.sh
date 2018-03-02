#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
ARCH=slc6_amd64_gcc630
CMSREL=CMSSW_9_4_0_pre3 
L1TTag=l1t-integration-v97.1
GT=94X_dataRun2_v2
sqlite1=303573
sqlite2=303835
curdir=$PWD
username=$USER

file=/store/data/Run2017F/ZeroBias/RAW/v1/000/306/091/00000/3E688D2E-FEBF-E711-9DBD-02163E019C1E.root
xrdcp -f root://cms-xrd-global.cern.ch/$file /tmp/$username.root


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
  --conditions=$GT -n -1 --data --no_exec --no_output \
  --filein=file:/tmp/$username.root \
  --fileout=${GT}.root \
  --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloStage2Params_2017_v1_8_2_updateHFSF_v6MET \
  --python_filename=l1Ntuple_${GT}.py

for sq in $sqlite1 $sqlite2; do
  if [! -f EcalTPG_${sq}_moved_to_1.db]; then
    wget http://cern.ch/ecaltrg/EcalLin/EcalTPG_${sq}_moved_to_1.db
  fi
  python ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
  cmsRun l1Ntuple_${GT}_${sq}.py >& l1Ntuple_${GT}_${sq}.py &
  #mv L1Ntuple.root l1Ntuple_${GT}_${sq}.root
done
################################


#----------------------------------------------------------------------------#
#                            Check out L1Menu code                           #
#----------------------------------------------------------------------------#
git clone -b EcalVal_940 https://github.com/cms-l1-dpg/L1Menu.git L1TriggerDPG/L1Menu
cd L1TriggerDPG/L1Menu/macros
cp $curdir/menulib.cc .
cp $curdir/menulib.hh .
cp $curdir/Prescale_Sets_RUN_306091_col_1.5.txt menu
cp $curdir/Selected_Seed.txt menu
make -j 8

#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#
for sq in $sqlite1 $sqlite2; do
  ./testMenu2016 -m menu/Prescale_Sets_RUN_306091_col_1.5.txt -l ${CMSSW_BASE}/src/l1Ntuple_${GT}_${sq}.root -o L1Menu_${GT}_${sq}_emu >& L1Menu_${GT}_${sq}_emu.log &
  ./testMenu2016 -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/l1Ntuple_${GT}_${sq}.root -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
done
wait
