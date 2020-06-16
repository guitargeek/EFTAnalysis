#!/bin/sh
echo $PWD

echo ${1}
#echo ${2}
#echo ${3}

cat>Job_${1}.csh<<EOF
#!/bin/tcsh
tar -xf CMSSW_10_2_22.tar.gz
source /cvmfs/cms.cern.ch/cmsset_default.csh
setenv SCRAM_ARCH slc7_amd64_gcc820 
cd CMSSW_10_2_22/src
scramv1 b ProjectRename
cmsenv
setenv SCRAM_ARCH "slc7_amd64_gcc820";
scramv1 b
cd -
sh triboson_production.sh -p pileup_files.txt -s WZZ_dim8 -c -o $PWD -a ${1} -n 1000
EOF

chmod 775 Job_${1}.csh


cat>condor_${1}.jdl<<EOF
universe = vanilla
Executable = Job_${1}.csh
Requirements = OpSys == "LINUX" && (Arch != "DUMMY" )
request_disk = 10000000
request_memory = 4000
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
transfer_input_files = Job_${1}.csh, triboson_production.sh, pileup_files.txt, CMSSW_10_2_22.tar.gz 
notification = Never
Output = CondorJobs/STDOUT_${1}.stdout
Error = CondorJobs/STDERR_${1}.stderr
Log = CondorJobs/LOG_${1}.log
x509userproxy = ${X509_USER_PROXY}
Queue 1
EOF

condor_submit condor_${1}.jdl
