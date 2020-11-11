import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/step2_ttbar_p8_03/step2_qcd8_1109.root'
#inputFiles_='file:../step_full.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/opendatadnn/step2_test.root'
#inputFiles_='file:/uscms/home/bburkle/nobackup/working_area/CMSSW_5_3_32/src/MLAnalyzer/test/step2_ttbarOD.root'
#inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/eusai/CRAB_UserFiles/step2_QCD600to3000_01/190213_183439/0000/step2_QCDPt_15_3000_Flat_V27_961.root'
# inputFiles_='root://cmsxrootd-site.fnal.gov//store/group/lpcml/CRAB_UserFiles/step2_ttbarOD_EmBj_01/190308_200019/0000/step2_OpenData_10.root'
#inputFiles_='file:step2_OpenData_10.root'
#inputFiles_='root://cmseos.fnal.gov//store/user/eusai/topgun01/200513_054612/0000/step_AODSIM_noPU_1.root'
#inputFiles_='file:/afs/cern.ch/user/d/ddyachko/CMSSW_5_3_32/src/pgun/step_AODSIM_noPU.root'
#inputFiles_='file:step_AODSIM_noPU_1.root'


inputFiles_='file:/eos/user/a/amslivar/ML_job_2020/200804_005034/0000/step_AODSIM_noPU_1.root'


#inputFiles_='root://eosuser.cern.ch//eos/user/d/ddicroce/ML/MassReg/200804_005034/0000/step_AODSIM_noPU_100.root'
#inputFiles_='file:/eos/user/d/ddicroce/ML/MassReg/200804_005034/0000/step_AODSIM_noPU_1.root'
#inputFiles_='root://cmseos.fnal.gov//store/user/eusai/topgun01/200513_054612/0000/step_AODSIM_noPU_1.root'
#inputFiles_='file:/afs/cern.ch/user/d/ddyachko/CMSSW_5_3_32/src/pgun/step_AODSIM_noPU.root'
#inputFiles_='file:step_AODSIM_noPU_1.root'
isTTbar_ = 1

maxEvents_=-1
#skipEvents_=0#
outputFile_ = 'test.root'
#outputFile_ = 'test/ttbar_new-production_test.root'

cmd="cmsRun %s inputFiles=%s outputFile=%s isTTbar=%d maxEv=%d" %(cfg,inputFiles_,outputFile_,isTTbar_,maxEvents_)
print '%s'%cmd
os.system(cmd)
