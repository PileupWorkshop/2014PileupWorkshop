#!/bin/bash

submit=0
if [ $# -gt 0 ]; then
  submit=1
fi
echo -n "submit: ${submit}"

nEvents=100
input=/mnt/nfs19/berta/PileupWorkshop2014/
outDir=/mnt/nfs19/berta/PileupWorkshop2014_output/
mkdir -p ${outDir}
executable="./subtraction"


#if [ ${submit} -eq 1 ]; then
# echo -n "Preparing submit files"
# source proxySetup.sh
#fi

DATE=`date '+%m%d_%H%M%S'`
logsDir=${outDir}/logs_${DATE}/
mkdir -p ${logsDir}

for npu in 140 100 60 30
#for npu in 30
do
  for sample in WW500 Zprime500 dijetsel20 dijetsel50 dijetsel100 dijetsel500
  #for sample in dijetsel20 dijetsel50 dijetsel100 dijetsel500
  #for sample in WW500 Zprime500
  do
    label=""
    if [[ "${sample}" == "dijetsel"* ]]; then
      label=sel
    fi
    arguments=" -hard ${input}/lhc14-pythia8-4C-${sample}-noUE-nev${label}1e5.pu14.gz -pileup ${input}/lhc14-pythia8-4C-minbias-nev1e7.pu14.gz -npu ${npu} -nev ${nEvents} -out ${outDir}/sub-${sample}-noUE-npu${npu}-radius0.4-massless.res -chs"
    if [ ${submit} -eq 1 ]; then
      submitFile=${logsDir}job.submit_${npu}_${sample}
      printf "executable = ${executable} \n" >> ${submitFile}
      printf "+MaxRuntime = 100000 \n" >> ${submitFile}
      printf "request_cpus = 1 \n" >> ${submitFile}
      printf "request_memory = 3GB \n" >> ${submitFile}
      printf "getenv = True \n" >> ${submitFile}
      printf "arguments  = ${arguments} \n" >> ${submitFile}
      printf "output = ${logsDir}/${npu}_${sample}.out \n" >> ${submitFile}
      printf "error = ${logsDir}/${npu}_${sample}.err \n" >> ${submitFile}
      printf "log = ${logsDir}/${npu}_${sample}.log \n" >> ${submitFile}
      printf "queue \n" >> ${submitFile}
      echo "Submitting"
      condor_submit ${submitFile}
    else
      ${executable}${arguments}
    fi

  done
done


echo -e "\n\nFinished $0\n"
