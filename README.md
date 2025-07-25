## Description of the analyzer
- This repository contains an analyzer for the ***ttHH Run III single leptonic channel***. 
- While the current code focuses on trigger studies, the overall framework is versatile enough for full analysis. 
- We strive to rely solely on ROOT libraries and are in the process of phasing out several “deprecated” libraries ( ex. HypothesisCombinatorics.h, MVAvarsJABDTthh.h, ... ), though this goal hasn’t been fully achieved yet.

## Pre-requisites
A Scientific Linux 7 is required for compiling this code. 

Note that lxplus does not support the el7 environment directly now, so you’ll need to use Singularity with the ```cmssw-el7``` command. If you want to use a different OS version, be sure to update the variables and libraries in **Makefile**.
- Important: Since Singularity doesn’t support job submissions to Condor, use lxplus8 or lxplus9 servers for submitting jobs.
- Libraries: ***ROOT, GSL, and LHAPDF***. cmssw-el7 and CMSSW_10_6_28 are suggested for this purpose.
- TTH libraries: Get this from the **Junghyun's eos storage** or download it from **his CERNBOX[https://cernbox.cern.ch/s/xPBQqigATEjgFQb]**


&#9655;	To set up the cmssw environment & install the analyzer:
```bash
# In lxplus..
cmssw-el7
cmsrel CMSSW_10_6_28
cd CMSSW_10_6_28/src && cmsenv
git clone https://github.com/G-Vian/ttHH_SLeptonic_Analyzer_RunIII.git
# Before using xrood protocol to copy TTH directory, please get the proxy first
xrdcp root://eosuser.cern.ch//eos/user/j/junghyun/public/TTH.tar.gz .
# If above line doesn't work, then download [ TTH.tar.gz ] at CERNBOX link:
#    --> https://cernbox.cern.ch/s/xPBQqigATEjgFQb

wget https://cernbox.cern.ch/remote.php/dav/public-files/xPBQqigATEjgFQb/TTH.tar.gz

tar -zxvf TTH.tar.gz && rm -rf TTH.tar.gz
mv TTH  ttHH_SLeptonic_Analyzer_RunIII/.
```

## Compilation to make execution file
```bash
cmssw-el7
cd <Path to Analyzer>/ttHH_SLeptonic_Analyzer_RunIII
cmsenv
wget https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp
git clone https://gitlab.cern.ch/cms-muonPOG/muonefficiencies.git
source setup.sh  # required for setup
make -j4

 

```
## All in one line (faster):
```bash
cmssw-el7
git clone https://github.com/G-Vian/ttHH_SLeptonic_Analyzer_RunIII.git && wget https://cernbox.cern.ch/remote.php/dav/public-files/xPBQqigATEjgFQb/TTH.tar.gz && tar -zxvf TTH.tar.gz && rm -rf TTH.tar.gz && mv TTH  ttHH_SLeptonic_Analyzer_RunIII/. && cd ttHH_SLeptonic_Analyzer_RunIII && wget https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp && git clone https://gitlab.cern.ch/cms-muonPOG/muonefficiencies.git && cmsenv && source setup.sh  && make -j4 
```


## Creating a Proxy
The proxy provides the necessary permissions for accessing grid jobs, Condor jobs, and samples on lxplus. If you’re a member of the **CERN CMS VO** with the required permissions, you can generate a proxy using the ```voms-proxy-init``` command.

If you specify an output file for the proxy with the ```--out``` option, set the path in the ```X509_USER_PROXY``` environment variable as shown below to activate proxy permissions:
```bash
voms-proxy-init --voms cms --valid 96:00 --out proxy.cert
export X509_USER_PROXY=proxy.cert
```
## All in one line (faster):
```bash

voms-proxy-init --voms cms --valid 96:00 --out proxy.cert && export X509_USER_PROXY=proxy.cert
```
Currently, to submit Condor jobs, a **proxy.cert** file with valid time remaining must be present in the analyzer directory. You can check the remaining valid time with the following command:
```bash
voms-proxy-info -file ./proxy.cert --timeleft
```

## Example Local Run
&#9655; To run locally, use the following syntax:
```bash
# ./<exe name> <path of the file list> <output name> <weight> <year> <MC or Data> <run name - just name it you want>
./ttHHanalyzer_trigger filelistTest/file_ttHH_0.txt ttHH.root 0.0000092849 2023 MC ttHH_MC_Test
./ttHHanalyzer_trigger filelistTest/ttHH_V15.txt ttHH_V15.root 0.000153 2024  MC ttHH_MC_V15
(small sample of 2024 with 50400 events)
```

## Running with Condor
A condor job typically requires a submit file, which sets various variables and environment configurations needed for the job, and an execution script that runs on the worker node. 

We use a single Python script to check for the existence of a proxy and to generate the submit file and execution script tailored to each sample.

- For Condor job, you’ll need configuration files in the ```AnalyzerConfig``` directory and a ```proxy.cert``` file within the analyzer directory. 
- You should also have ```sample lists files``` generated by ```DAS query```.
Based on these files, ```submit_job_FH_Trigger.py``` will split and submit jobs. When you create a condor directory, submission files, execution scripts, and logs are automatically generated inside it.
- Before you run the ```submit_job_FH_Trigger.py```, Adjust the variables at the top of the Python script to suit your situation.

&#9655; To submit condor job (should be done out of the cmssw-el7 enviroment):
```bash
python3 submit_job_FH_Trigger.py
```

