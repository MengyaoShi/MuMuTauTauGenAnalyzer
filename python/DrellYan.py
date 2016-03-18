import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('/afs/cern.ch/user/m/mshi/CMSSW_7_4_12_patch4/src/GGHAA2Mu2TauAnalysis/DrellYan.txt')
process = cms.Process("DrellYan")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(*mylist))


process.DrellYan = cms.EDAnalyzer('DrellYanAnalyzer',
    genParticleTag=cms.InputTag("genParticles"),
      outFileName = cms.string('/afs/cern.ch/user/m/mshi/CMSSW_7_6_3/src/GGHAA2Mu2TauAnalysis/MuMuTauTauGenAnalyzer/python/DrellYan_out.root')
)
process.TFileService = cms.Service("TFileService",
                                        fileName = cms.string('Drell_Yan.root')
)

process.p = cms.Path(process.DrellYan)
