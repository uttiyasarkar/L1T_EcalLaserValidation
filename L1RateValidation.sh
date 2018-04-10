#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
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
filelist=('/store/data/Run2017F/ZeroBias1/RAW/v1/000/306/091/00000/D050FAEF-FFBF-E711-928B-02163E019C6C.root'
'/store/data/Run2017F/ZeroBias2/RAW/v1/000/306/091/00000/B4E9FAEF-FFBF-E711-928B-02163E019C6C.root'
'/store/data/Run2017F/ZeroBias3/RAW/v1/000/306/091/00000/803A9506-00C0-E711-AB01-02163E01A4FC.root'
'/store/data/Run2017F/ZeroBias4/RAW/v1/000/306/091/00000/0C719606-00C0-E711-AB01-02163E01A4FC.root'
'/store/data/Run2017F/ZeroBias5/RAW/v1/000/306/091/00000/D075606D-01C0-E711-999A-02163E0133B4.root'
'/store/data/Run2017F/ZeroBias6/RAW/v1/000/306/091/00000/96406B6D-01C0-E711-999A-02163E0133B4.root'
'/store/data/Run2017F/ZeroBias7/RAW/v1/000/306/091/00000/2480F615-00C0-E711-BD6F-02163E0129DD.root'
'/store/data/Run2017F/ZeroBias8/RAW/v1/000/306/091/00000/18C38216-00C0-E711-BD6F-02163E0129DD.root'
'/store/data/Run2017F/ZeroBias1/RAW/v1/000/306/091/00000/E6E98AF1-FFBF-E711-B2D3-02163E019DA3.root'
'/store/data/Run2017F/ZeroBias2/RAW/v1/000/306/091/00000/EE768CF1-FFBF-E711-B2D3-02163E019DA3.root'
'/store/data/Run2017F/ZeroBias3/RAW/v1/000/306/091/00000/4A224649-00C0-E711-8018-02163E01A59E.root'
'/store/data/Run2017F/ZeroBias4/RAW/v1/000/306/091/00000/4A8E5749-00C0-E711-8018-02163E01A59E.root'
'/store/data/Run2017F/ZeroBias5/RAW/v1/000/306/091/00000/BC32EB7F-01C0-E711-9F4B-02163E011E9C.root'
'/store/data/Run2017F/ZeroBias6/RAW/v1/000/306/091/00000/E629BF80-01C0-E711-9F4B-02163E011E9C.root'
'/store/data/Run2017F/ZeroBias7/RAW/v1/000/306/091/00000/0039A9F6-FFBF-E711-B46B-02163E01472B.root'
'/store/data/Run2017F/ZeroBias8/RAW/v1/000/306/091/00000/3CDAADF6-FFBF-E711-B46B-02163E01472B.root'
'/store/data/Run2017F/ZeroBias1/RAW/v1/000/306/091/00000/28CF6607-00C0-E711-9A44-02163E019D02.root'
'/store/data/Run2017F/ZeroBias2/RAW/v1/000/306/091/00000/88127307-00C0-E711-9A44-02163E019D02.root'
'/store/data/Run2017F/ZeroBias3/RAW/v1/000/306/091/00000/3426DDF1-FFBF-E711-B118-02163E0143E4.root'
'/store/data/Run2017F/ZeroBias4/RAW/v1/000/306/091/00000/241701F2-FFBF-E711-B118-02163E0143E4.root'
'/store/data/Run2017F/ZeroBias5/RAW/v1/000/306/091/00000/7ED3F5CB-01C0-E711-8BE8-02163E01A223.root'
'/store/data/Run2017F/ZeroBias6/RAW/v1/000/306/091/00000/6440F8CB-01C0-E711-8BE8-02163E01A223.root'
'/store/data/Run2017F/ZeroBias7/RAW/v1/000/306/091/00000/C402BF0C-00C0-E711-99CE-02163E01A607.root'
'/store/data/Run2017F/ZeroBias8/RAW/v1/000/306/091/00000/F018C30C-00C0-E711-99CE-02163E01A607.root' )
#----------------------------------------------------------------------------#
#                            Getting the reference                           #
#----------------------------------------------------------------------------#
#it may be gzipped infact ...
## prevent exit from failed wget
set +e 
wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz 
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

scram b -j $((`nproc`/2))

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
  --filein=`echo $(IFS=, ; echo "${filelist[*]}")` \
  --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloStage2Params_2017_v1_8_2_updateHFSF_v6MET \
  --python_filename=l1Ntuple_${GT}.py

for sq in $sqs; do
  if [ ! -f EcalTPG_${sq}_moved_to_1.db ]; then
    wget http://cern.ch/ecaltrg/EcalLin/EcalTPG_${sq}_moved_to_1.db
  fi
  python ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
  cmsRun l1Ntuple_${GT}_${sq}.py >& l1Ntuple_${GT}_${sq}.log  &
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
make comparePlots

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to checkout and compile code: %.6f minutes" $dur
#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#
echo "Waiting for Ntuple production to finish......"
wait $pids

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur

for sq in $sqs; do
  ./testMenu2016 -m menu/Prescale_Sets_RUN_306091_col_1.5.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.root -o L1Menu_${GT}_${sq}_emu >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.root -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
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
  tar -xzvf $curdir/L1T_EcalLaserValidation/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz -C results/
fi

python CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log
./comparePlots results/L1Seed*root

#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#
#WORKSPACE=${curdir} #local test
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
