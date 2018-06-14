#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
ARCH=slc6_amd64_gcc630
CMSREL=CMSSW_10_1_5
L1TTag=l1t-integration-v98.4
GT=101X_dataRun2_Prompt_v9
nproc=`nproc`
sqlite1=$1 ##ref
sqlite2=$2
week=$3
year=$4
curdir=$PWD
username=$USER
pids=""
hasref=false
## ZeroBias Raw, Run Run316766, LS51-85
filelist=('/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/047A9646-9764-E811-88FD-FA163E0BD5CE.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/06D477F8-4E5E-E811-9072-FA163E17588A.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/0CD2F426-4F5E-E811-9A25-FA163EFF76BB.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/121D09F4-4E5E-E811-8C24-FA163E10A7AB.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/16B8B270-9764-E811-AB9A-FA163E0B1665.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/3C6DB425-9764-E811-ABDC-FA163EAE126B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/488F1418-4F5E-E811-8B2A-FA163EB23FD6.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/5A81FB1F-9764-E811-99DB-FA163E0AB4C0.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/605D75A8-4F5E-E811-BFC4-FA163E8DAC79.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/74BE3278-9764-E811-AF3E-02163E019EE2.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/86B5BE1F-9764-E811-B496-FA163EE16015.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/8A3D3A27-9764-E811-B649-FA163E45F42B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/8EA92125-9764-E811-B493-FA163E10317B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/8EFB90DD-505E-E811-BA3F-FA163EA8F941.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/92FA79E4-4E5E-E811-8635-FA163EDB96EE.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/94114EDC-505E-E811-BF0D-FA163E741841.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/9CC93029-9764-E811-AC7B-FA163E9BD3B9.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/AA020022-9764-E811-8ADA-FA163EE2B40B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/B6BD6324-9764-E811-A808-FA163E646768.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/BE80BE53-9764-E811-87FC-FA163E793AC1.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/E2A315AE-4F5E-E811-AB65-FA163E5355C6.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/E8451A20-4F5E-E811-B25F-FA163EBC9CF5.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/E8C8A507-4F5E-E811-B03D-02163E019F53.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run316766_LS51_80/F0459C1C-505E-E811-82A7-FA163E25B672.root')
#----------------------------------------------------------------------------#
#                            Getting the reference                           #
#----------------------------------------------------------------------------#
#it may be gzipped infact ...
## prevent exit from failed wget
set +e 
wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/${sqlite1}/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz 
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
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
scramv1 project CMSSW $CMSREL
cd $CMSREL/src
eval `scramv1 runtime -sh`
git-cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline
git cms-merge-topic -u cms-l1t-offline:$L1TTag
git cms-addpkg L1Trigger/L1TCommon
git cms-addpkg L1Trigger/L1TMuon
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TMuon.git L1Trigger/L1TMuon/data
git cms-addpkg L1Trigger/L1TCalorimeter
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data

scram b -j ${nproc}

dur=$(echo "$(date +%s.%N) - $starttime" | bc)
printf "Execution time to L1T checkout: %.6f seconds" $dur
#----------------------------------------------------------------------------#
#                          Running the TP emulation                          #
#----------------------------------------------------------------------------#
echo running $GT
cmsDriver.py l1Ntuple -s RAW2DIGI --era=Run2_2017  \
  --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU \
  --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimEcalTP \
  --conditions=$GT -n -1 --data --no_exec --no_output  \
  --filein=inputFiles \
  --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloStage2Params_2017_v1_8_2_updateHFSF_v6MET \
  --python_filename=l1Ntuple_${GT}.py

Nsq=`echo $sqs | awk -F ' ' '{print NF}'`
Nfiles=${#filelist[@]}
NfpJ=`echo "${Nfiles} *${Nsq}/${nproc}" | bc`
NJ=`echo "${Nfiles}/${NfpJ}" | bc`
for sq in $sqs; do
  if [ ! -f EcalTPG_${sq}_moved_to_1.db ]; then
    wget http://cern.ch/ecaltrg/EcalLin/EcalTPG_${sq}_moved_to_1.db
  fi
  python ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
  for ((i = 0; i < $NJ; i++)); do
    let cnt1=$(($i*$NfpJ))
    args=`printf "inputFiles=%s " "${filelist[@]:$cnt1:$NfpJ}"`
    args+=`echo outputFile=L1Ntuple_${GT}_${sq}_${i}.root`
    cmsRun l1Ntuple_${GT}_${sq}.py `echo $args` >& l1Ntuple_${GT}_${sq}_${i}.log  &
    pids="$pids $!"
  done
done
echo "Waiting for Ntuple production to finish......"
wait $pids
dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur

for sq in $sqs; do
  ls $PWD/L1Ntuple_${GT}_${sq}_*.root > L1Ntuple_${GT}_${sq}.list
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
cp $curdir/Prescale_2018_v1_0_0_Col_2.0.txt menu/
cp $curdir/Selected_Seed.txt menu/
make -j ${nproc}
make comparePlots

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to checkout and compile code: %.6f minutes" $dur
#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#

for sq in $sqs; do
  ./testMenu2016 -m menu/Prescale_2018_v1_0_0_Col_2.0.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
  pids="$pids $!"
done
echo "Waiting for menu rate estimation to finish......"
wait $pids
dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur

#----------------------------------------------------------------------------#
#                                Compare rate                                #
#----------------------------------------------------------------------------#

if $hasref; then
  tar -xzvf $curdir/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz -C results/
fi

python CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log
./comparePlots results/L1Seed*root

#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#
if [ -z "$WORKSPACE" ]
then
  WORKSPACE=${curdir}
fi

if [ -f ${WORKSPACE}/upload/$2 ]
then
  echo "dir is already existing"
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
else
  mkdir -p ${WORKSPACE}/upload/$2
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
fi

#we need to make a tar gz of this one
cp results/L1Menu_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp results/L1Seed_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp results/L1Seed_${GT}_${sqlite2}_emu.root ${WORKSPACE}/upload/${2}/
cp ${sqlite2}.log ${WORKSPACE}/upload/${2}/
cp compRate.csv ${WORKSPACE}/upload/${2}/

for i in nDiEGVsPt.gif nEGVsEta.gif nEGVsPt.gif nETMVsETM.gif nIsoEGVsPt.gif nJetVsPt.gif nETMVsETMHF.gif nQuadCenJetVsPt.gif nTauVsPt.gif; do
  cp  results/comparePlots/Rate/${i}  ${WORKSPACE}/upload/${2}/
done
cd ${WORKSPACE}/upload/${2}
tar -czvf ${WORKSPACE}/upload/L1TEcalValidation_${year}_${week}_${sqlite2}.tgz *

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time of workflow: %.6f minutes" $dur
