#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
ARCH=slc6_amd64_gcc700
CMSREL=CMSSW_10_2_1
L1TTag=l1t-integration-v100.1
GT=102X_dataRun2_Prompt_v4
Prescale=Prescale_2018_v2_0_0_Col_2.0.txt
nproc=`nproc`
sqlite1=$1 ##ref
sqlite2=$2
week=$3
year=$4
curdir=$PWD
username=$USER
pids=""
hasref=false
## ZeroBias Raw, Run Run322079, LS1-141
filelist=('/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/0034CBFE-5DAF-E811-B9F8-FA163E94E4CF.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/04A298EC-5DAF-E811-AD1D-FA163E3ECB15.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/0CD3CD19-5EAF-E811-882E-FA163E61340F.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/127E4AAF-5EAF-E811-9586-02163E0176DA.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/1C0B6639-5EAF-E811-97C4-02163E0149E7.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/1CDDDE20-5EAF-E811-852B-FA163EF8848D.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/1EB214F1-5DAF-E811-91B0-FA163E9E0180.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/200DADC3-5EAF-E811-A0BD-02163E01305C.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/2A8B954E-A5AE-E811-87A6-FA163E28CFAA.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/2EB57257-5EAF-E811-8787-02163E01A048.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/3CC3A0EF-5DAF-E811-9059-A4BF012CBE43.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/3E452302-5EAF-E811-9C92-FA163E66FF50.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/4601D346-5FAF-E811-99ED-02163E01763A.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/48EAF8F5-5DAF-E811-B2BB-FA163EA222CE.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/4EF8D5D5-5FAF-E811-9333-FA163EF8E861.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/5401C114-5EAF-E811-8CE0-FA163E204525.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/563635A5-5EAF-E811-8EF3-02163E014A17.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/5E527B09-5EAF-E811-BEA6-FA163E841855.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/66E3EAF1-5DAF-E811-A9F5-FA163EE723F8.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/6E5C5312-5EAF-E811-B160-FA163E65DEE3.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/72718F53-A5AE-E811-ADDC-FA163EF225E6.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/72D370F7-5DAF-E811-B2D5-FA163EE59243.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/74E087F4-5DAF-E811-9CAD-FA163E9423DA.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/7696F3FD-5DAF-E811-8681-FA163EE10A6B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/7A7E56FA-5DAF-E811-AFCF-FA163E841855.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/7C4826F0-5DAF-E811-82D3-A4BF0112DB74.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/8092B8F1-5DAF-E811-B018-FA163EF6AA30.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/88BFD5C6-5EAF-E811-9479-02163E00F6C4.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/9CAF7E05-5FAF-E811-B3D1-02163E014E4F.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/A6A551ED-5DAF-E811-80EC-FA163ED7F3BC.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/AEFF53ED-5DAF-E811-8C3F-FA163EB2605F.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/B2B03EA2-5EAF-E811-8E7F-FA163E55495F.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/B82DEEF2-5DAF-E811-ACF7-FA163E86F14F.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/BA313FF1-5DAF-E811-A018-FA163E08DA61.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/BC0FF162-5FAF-E811-9FB4-02163E00F7AA.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/C653D946-5EAF-E811-BD17-02163E010F0D.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/CE11AB0E-5EAF-E811-8C4C-FA163E0DF5F2.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/D47E06F5-5DAF-E811-B420-FA163E9666C9.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/D8BB2AF4-5DAF-E811-B69F-FA163E804F61.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/DAF8B51F-5EAF-E811-ABD1-FA163E23298B.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/DEA94AFE-5DAF-E811-A5D3-FA163E0512AE.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/DEA94AFE-5DAF-E811-A9E5-FA163E0512AE.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/E0C706F0-5DAF-E811-8811-A4BF01277567.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/E2245A15-5EAF-E811-B4E6-FA163E28EA0E.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/EC2E09FD-5DAF-E811-AF7F-FA163E816F3E.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/EC88081B-5FAF-E811-96CF-02163E017654.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/F6CFB58A-5FAF-E811-AC9B-02163E00C918.root'
'/store/group/dpg_trigger/comm_trigger/L1Trigger/L1T_EcalValidation/Raw_Run322079_LS1_100/FC8F3F23-5EAF-E811-8D65-02163E019F3A.root')

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
cmsDriver.py l1Ntuple -s RAW2DIGI --era=Run2_2018  \
  --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU \
  --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimEcalTP \
  --conditions=$GT -n -1 --data --no_exec --no_output  \
  --filein=inputFiles \
  --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_2 \
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
  cp $PWD/L1Ntuple_${GT}_${sq}_*.log ${WORKSPACE}/upload/${2}/
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
cp $curdir/$Prescale menu/
cp $curdir/Selected_Seed.txt menu/
make -j ${nproc}
make comparePlots

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to checkout and compile code: %.6f minutes" $dur
#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#

for sq in $sqs; do
  ./testMenu2016 -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
  pids="$pids $!"
done
echo "Waiting for menu rate estimation to finish......"
wait $pids
dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur
cp L1Menu_${GT}_*_emu.log ${WORKSPACE}/upload/${2}/
cp L1Seed_${GT}_*_emu.log ${WORKSPACE}/upload/${2}/

#----------------------------------------------------------------------------#
#                                Compare rate                                #
#----------------------------------------------------------------------------#

if $hasref; then
  tar -xzvf $curdir/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz -C results/
fi

ls results/

python CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log
./comparePlots results/L1Seed*root

#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#

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
