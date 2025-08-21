
#!/bin/bash

#how to use :
#  chmod +x makeFileList.sh
# voms-proxy-init --voms cms
# ./makeFileList.sh


#2022
dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_22.txt
dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_22EE.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_22.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_22EE.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_22.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_22EE.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22.txt
#dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22.txt
#dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22EE.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_22EE.txt

dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_22EE.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_22.txt

#dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_postEE_v4-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_postEE_v4_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22EE.txt

#dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-JMENano12p5_132X_mcRun3_2022_realistic_v3_ext1-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_22.txt



#dasgoclient --query="file dataset=//ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22.txt
#dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv11-126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22EE.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_22EE.txt

#2023
dasgoclient --query="file dataset=/TTHHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM " | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_23.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_23.txt
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_23.txt
dasgoclient --query="file dataset=/TTZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_23.txt
dasgoclient --query="file dataset=/TTZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_23.txt
dasgoclient --query="file dataset=/TTZ-ZtoQQ-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZ_23.txt

#2024 v15 (for 2024 samples, other redirector is necessary to TTHH and TTZZ)
dasgoclient --query="file dataset=/TTHH-HHto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM" | sed 's|^|root://cms03.lcg.cscs.ch:1094//|' >> full_TTHH_24_15.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_24_15.txt
dasgoclient --query="file dataset=/TTZH-ZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_24_15.txt
dasgoclient --query="file dataset=/TTZZ-ZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v3/NANOAODSIM" | sed 's|^|root://eoscms.cern.ch:1094//|' >> full_TTZZ_24_15.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_24_15.txt
dasgoclient --query="file dataset=/TTBBtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTbb_24_15.txt
dasgoclient --query="file dataset=/TTH-Hto2B_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_24_15.txt



#2024 v14 (for 2024 samples, other redirector is necessary to TTHH and TTZZ) --> Have to figure out correct servers for TTHH and TTZZ V14 
dasgoclient --query="file dataset=/ttHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/Run3Winter24NanoAOD-133X_mcRun3_2024_realistic_v8-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTH_24_14.txt
#dasgoclient --query="file dataset=/TTHH-HHto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTHH_24_14.txt
dasgoclient --query="file dataset=/TT4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TT4b_24_14.txt
dasgoclient --query="file dataset=/TTZH-ZHto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v4/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZH_24_14.txt
#dasgoclient --query="file dataset=/TTZZ-ZZto4B_TuneCP5_13p6TeV_madgraph-pythia8/RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v3/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTZZ_24_14.txt
dasgoclient --query="file dataset=/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24NanoAOD-140X_mcRun3_2024_realistic_v26-v2/NANOAODSIM" | sed 's|^|root://xrootd-cms.infn.it//|' >> full_TTSL_24_14.txt



