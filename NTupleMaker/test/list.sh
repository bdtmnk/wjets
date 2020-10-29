
dir=$1
year=$2
tp=$3

#rm datasets${dir}
#sources="/nfs/dust/cms/user/rasp/storage/76x_JECv2_MVAMET0p6/ /nfs/dust/cms/group/susy-desy/Run2/Stau/MC/25ns/76x_JECv2_MVAMET0p6/"
#sources="/nfs/dust/cms/user/bobovnii/CMSSW_8_0_20/src/DesyTauAnalyses/NTupleMaker/NtupleForWmnu"
#sources="/nfs/dust/cms/user/rasp/storage/SingleElectron_Run2016-23Sep2016/"
sources="/pnfs/desy.de/cms/tier2/store/user/ibabouni/higgs-kit/MC2017/"
sources="/pnfs/desy.de/cms/tier2/store/user/telenz/13TeV/NTuples/MC/RunIIFall17MiniAOD-94X_mc2017_realistic/"
sources="/nfs/dust/cms/user/rasp/ntuples/MC_2017"
sources="/nfs/dust/cms/group/higgs-kit/rasp/MC_2017"
sources="/nfs/dust/cms/user/rasp/ntuples/MC2017_newMET"
sources="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2/"
sources="/nfs/dust/cms/user/tlenz/13TeV/2017/NTuples/MC/12Apr2018_PU2017_METrecipe_v2"
sources="/nfs/dust/cms/user/rasp/ntuples/Data2017_newMET_v2/"
sources="/nfs/dust/cms/user/tlenz/13TeV/2017/NTuples/MC/12Apr2018_PU2017_METrecipe_v2/"
sources="/nfs/dust/cms/user/mameyer/SM_HiggsTauTau/ntuples/METRecipev2/"

#MC TauId Feb 2019:
#sources="/nfs/dust/cms/user/cardinia/gridjobs/2018/NTuples/MC/RunIIAutumn18/"
#sources="/nfs/dust/cms/user/ywen/Storage/MC/"
#sources="/nfs/dust/cms/user/cardinia/gridjobs/2018/NTuples/MC/RunIIAutumn18"
#sources="/nfs/dust/cms/user/ywen/Storage/SingleMuon/"
#sources="/nfs/dust/cms/group/higgs-kit/2018/MC_RunIIAutumn18MiniAOD/"
#sources="/nfs/dust/cms/group/higgs-kit/2018/MC_RunIIAutumn18MiniAOD/"
#sources="/nfs/dust/cms/user/cardinia/gridjobs/2018/NTuples/MC/RunIIAutumn18/"

#2016
if [[ ${year} == "2016" ]]  &&   [[ ${tp} == "data" ]]; then
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/data/"
fi

if [[ ${year} == "2016" ]]  && [[ ${tp} == "mc" ]]; then
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2016/mc/"
fi

#2017

if [[ ${year} == "2017" ]]  &&  [[ ${tp} == "data" ]]; then
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/" 
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/data/SingleMuon/"
fi

if [[ ${year} == "2017" ]]  &&  [[ ${tp} == "mc" ]]; then
#sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3/" 
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v4/" 
#sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2017/mc_v3/QCD/"
fi

#2018
if [[ ${year} == "2018" ]]  &&  [[ ${tp} == "data" ]]; then
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/SingleMuon/" 
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/data/"
fi

if [[ ${year} == "2018" ]]  &&  [[ ${tp} == "mc" ]]; then
sources="/pnfs/desy.de/cms/tier2/store/user/mmeyer/ntuples/2018/mc/"
sources="/pnfs/desy.de/cms/tier2/store/user/ldidukh/Missed_MC_v4/"
#sources="/nfs/dust/cms/user/dydukhle/WenuDatasets/"
sources="/pnfs/desy.de/cms/tier2/store/user/ldidukh/QCD_EM_Enriched_v2/"
#sources="/pnfs/desy.de/cms/tier2/store/user/ldidukh/Missed_MC_v2/"
fi

alias ls='ls'


for source in $sources
do
for i in `ls $source/`
do
	ls $source/$i/*.root > ${dir}/${year}/${tp}/$i
	echo $i >> list_${dir}_${year}_${tp}
	echo hadd $i.root ${i}_*.root >> ${dir}/${year}/${tp}/merg.sh
	echo rm ${i}_* >> ${dir}/${year}/${tp}/delete.sh
	echo "" >> ${dir}/${year}/${tp}/merg.sh

done
done

