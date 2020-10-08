import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

import sys
for arg in sys.argv[2:]:
        opt = arg[0:arg.find("=")]
        if opt=="name":
         name = str(arg[arg.find("=")+1:])

process = cms.Process("CTPPSTestProtonReconstruction", eras.run2_miniAOD_devel)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
#    "root://eostotem.cern.ch//eos/cms/store/data/Run2016B/ZeroBias/AOD/UL16_ver2_forHarvestOnly_rsb-v1/280000/7BBCD84B-4829-1A47-920E-4203C23118E3.root"
    "root://eostotem.cern.ch//eos/cms/store/group/phys_pps/reconstruction/2016/physics_runs/version-UL-2/"+name+"_xangle185_beta0.40_DoubleEG.root"
  )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_dataRun2_v26")

# number of events to process
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)


process.demo = cms.EDAnalyzer("raddam_2016",
  tagRecHit = cms.InputTag("totemRPRecHitProducer"),
  tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer"),
  clusterSize_a = cms.double(0.02),
  clusterSize_b = cms.double(0.3)
)

process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.path = cms.Path(process.demo)
