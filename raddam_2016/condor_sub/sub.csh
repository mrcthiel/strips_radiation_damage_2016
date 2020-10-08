#!/bin/tcsh
tar -zxvf CMSSW.tar.gz
tar -zxvf raddam_2016.tar.gz
mv raddam_2016 CMSSW_10_6_12/src/

setenv SCRAM_ARCH slc7_amd64_gcc820

cd CMSSW_10_6_12/src

scramv1 b ProjectRename

source /cvmfs/cms.cern.ch/cmsset_default.csh

setenv SCRAM_ARCH "slc7_amd64_gcc820";

cmsenv

scramv1 b

cd -

cmsRun config_cfg.py name=${1}
mv output.root ${1}.root
mv *.root /afs/cern.ch/work/m/mthiel/private/PPS_efficiency_2019/04_02_2020/optimization/CMSSW_10_6_12/src/raddam_2016/raddam_2016/condor_sub/

