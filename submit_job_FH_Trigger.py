
#!/usr/bin/env python
# Attention: this works only for versions of Python below 3.7! (not necessarily always true)
# how to use : first you should be out of cmssw-el7, then use the command : python3 submit_job_FH_Trigger.py 
# how to watch the job status: " condor_q " for all jobs or " condor_q  <job ID>  " for an specific job
# constant monitoring: " watch -n 5 condor_q " for all jobs or  " watch -n 5 condor_q <job ID> " for an specific job
# cancelling specific job : condor_rm <job ID>
# cancelling all jobs of the user : condor_rm user
import os
import sys
import time
import subprocess
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(script_dir, 'python', 'ttHHmodules'))
import ProxyChecker
import resultChecker


class CondorJobManager(object):

    def __init__(self):

        self.time_info = time.strftime("%Y%m%d-%H%M%S")

        # Variables for jobs, Please check before you run this script
        self.analyzer_path = script_dir
        self.nameofExe = "ttHHanalyzer_trigger" # Name of compiled execution file to run analyzer
        self.path_output_base = "/eos/user/g/gvian/job"
        self.os_version = "el7"
        self.memorySize = "10 GB"
        self.jobFlavour = "tomorrow"

        self.config_file_path = os.path.join(self.analyzer_path, "AnalyzerConfig/TriggerEffStudy_2017_FH.txt")
        self.proxy_path = os.path.join(self.analyzer_path, "proxy.cert")
        self.condor_files_path = os.path.join(self.analyzer_path, "condor/filePath_Trigger")
        self.sample_list_path = os.path.join(self.analyzer_path, "filelistTrigger")

        self.make_directory(self.condor_files_path)
        self.make_directory(self.path_output_base, 777) # 2nd argument is for the permissionsn

        try:
            print("Setted Path of proxy : {0}".format(self.proxy_path))
            proxy_checker = ProxyChecker.ProxyChecker(self.proxy_path)
            proxy_checker.check()
        except Exception as e:
            print("{0}".format(e)) # Error Messages which are defined in ProxyChecker.py
            sys.exit(1)

        self.print_memory_status()
        self.process_config_file()


    def make_directory(self, path, permission=755):
        try:
            subprocess.call(['mkdir', '-p', path])
        except subprocess.CalledProcessError as e:
            print("  Error creating drectory {0} : {1}".format(path, e))
            sys.exit(1)
       
        if permission == 777:
            chmod_cmd = ['chmod', '777', path]
            result = subprocess.Popen(chmod_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = result.communicate()
            if result.returncode != 0:
                print("Error setting permissions on output directory: {0}".format(stderr))
                sys.exit(1)
            else:
                print("Set permissions to 777 for {0}".format(path))


    def process_config_file(self):
        with open(self.config_file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue # Skep empty lines or comments

                try:
                    self.parse_config_line(line)
                    self.prepare_output_directory()
                    self.setup_and_submit_job()

                except Exception as e:
                    print("  Error with process of the condor job submition with {0} : {1}".format(line, e))
                    continue


    def parse_config_line(self, line):
        # Parse the configuration line into individual components
        parts = line.strip().split()
        if len(parts) < 6:
            raise ValueError("Invalid configuration line (expected at least 6 fields) : {0}".format(line))

        self.file_list_name = parts[0]
        self.output_dir = parts[1]
        self.weight = parts[2]
        self.year = parts[3]
        self.data_or_mc = parts[4]
        self.sample_name = parts[5]

        self.path_output = os.path.join(self.path_output_base, self.output_dir + "/")

        self.script_name = os.path.join(self.condor_files_path, "run_{0}.sh".format(self.output_dir))
        self.condor_submit_name = os.path.join(self.condor_files_path, "{0}_condor.sub".format(self.output_dir))
        self.arg_list_file = os.path.join(self.condor_files_path, "arguments_{0}.txt".format(self.output_dir))
        self.tmp_folder = os.path.join(self.condor_files_path, "tmp_{0}_{1}".format(self.output_dir, self.time_info))
        self.make_directory(self.tmp_folder)


    def print_memory_status(self):

        print("\n--------------------------------------------------------------------------\n")
        print("  --> Memory Check before submit jobs <--  ")

        # Print system memory status
        print("\nSystem Memory Status:")
        subprocess.call(['free', '-h'])

        ### Print disk usage status
        ##print("\nDisk Usage Status:")
        ##subprocess.call(['df', '-h'])
        
        # Print User Home directory space
        print("\nUser Home Memory Status:")
        subprocess.call(['fs', 'lq'])

        # Check EOS storage space
        print("\nEOS Storage Space:")
        eos_command = ['eos', 'quota']
        try:
            subprocess.call(eos_command)
        except Exception as e:
            print("Error checking EOS storage space: {0}".format(e))

        print("\n--------------------------------------------------------------------------\n")


    def prepare_output_directory(self):
        # Check if the output directory exists in EOS
        # Since the output directory is in EOS, we need to use EOS commands to check and create it
        eos_base_cmd = ['eos', 'root://eosuser.cern.ch']

        # Check if the directory exists
        check_dir_cmd = eos_base_cmd + ['ls', self.path_output]
        result = subprocess.Popen(check_dir_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = result.communicate()
        if result.returncode != 0:
            # Directory does not exist, create it
            print("Creating EOS output directory: {0}".format(self.path_output))
            mkdir_cmd = eos_base_cmd + ['mkdir', '-p', self.path_output]
            result = subprocess.Popen(mkdir_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = result.communicate()
            if result.returncode != 0:
                print("Error creating output directory: {0}".format(stderr))
                sys.exit(1)
        else:
            print("EOS output directory exists: {0}".format(self.path_output))

        # Set permissions to 755
        chmod_cmd = eos_base_cmd + ['chmod', '755', self.path_output]
        result = subprocess.Popen(chmod_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = result.communicate()
        if result.returncode != 0:
            print("Error setting permissions on output directory: {0}".format(stderr))
            sys.exit(1)
        else:
            print("Set permissions to 755 for {0}".format(self.path_output))


    def generate_argument_list(self):
        # Create the argument list file for Condor submission
        with open(self.arg_list_file, "w") as argout:
            count = 0
            # Read the sample list file and split it into per-job files
            sample_list_file_path = os.path.join(self.sample_list_path, self.file_list_name)
            with open(sample_list_file_path, "r") as sample_list_in:
                lines = [line.strip() for line in sample_list_in if line.strip() and not line.startswith('#')]
                for line in lines:
                    # For each line, create a per-job filelist in the tmp_folder
                    per_job_filelist_name = "filelist_{0}_{1}.txt".format(self.output_dir, count)
                    per_job_filelist_path = os.path.join(self.tmp_folder, per_job_filelist_name)
                    with open(per_job_filelist_path, 'w') as per_job_filelist:
                        per_job_filelist.write(line + '\n')

                    # Write the arguments for this job to the argument list file
                    argout.write("{0} {1}_{2}.root {3} {4} {5} {6}\n".format(
                        per_job_filelist_path, self.output_dir, count, 
                        self.weight, self.year, self.data_or_mc, self.sample_name))
                    count += 1

    
    def write_condor_submission_file(self):
        # Create the Condor submission script
        with open(self.condor_submit_name, 'w') as f:
            f.write("x509userproxy           = {0}\n".format(self.proxy_path))
            f.write("getenv                  = True\n")
            f.write("executable              = {0}\n".format(self.script_name))
            f.write("arguments               = $(args)\n")
            f.write("output                  = {0}\n".format(os.path.join(self.tmp_folder, 'job_{0}.$(ClusterId).$(ProcId).out'.format(self.sample_name))))
            f.write('MY.WantOS               = "{0}"\n'.format(self.os_version))
            f.write("MY.XRDCP_CREATE_DIR     = True\n")
            f.write("error                   = {0}\n".format(os.path.join(self.tmp_folder, 'error_{0}.$(ClusterId).$(ProcId).err'.format(self.sample_name))))
            f.write("log                     = {0}\n".format(os.path.join(self.condor_files_path, 'log_{0}.$(ClusterId).log'.format(self.sample_name))))
            f.write("request_memory          = {0}\n".format(self.memorySize))
            f.write('+JobFlavour             = "{0}"\n'.format(self.jobFlavour))
            f.write("queue args from {0}\n".format(self.arg_list_file))

    def create_executable_script(self):
        # Create the executable script that Condor will run
        with open(self.script_name, 'w') as fout:
            fout.write("#!/bin/sh\n")
            fout.write("echo\n")
            fout.write("echo\n")
            fout.write("echo 'START---------------'\n")
            fout.write("echo 'WORKDIR ' ${PWD}\n")
            fout.write("source \"/afs/cern.ch/cms/cmsset_default.sh\"\n")
            fout.write('cd "{0}"\n'.format(self.analyzer_path))
            fout.write("cmsenv\n")
            fout.write("echo 'WORKDIR ' ${PWD}\n")
            fout.write('source "{0}/setup.sh"\n'.format(self.analyzer_path))
            # Prepare output directory in EOS
            fout.write("eos root://eosuser.cern.ch mkdir -p {0}\n".format(self.path_output))
            ##fout.write("eos root://eosuser.cern.ch chmod 777 {0}\n".format(self.path_output))
            fout.write('"{0}/{1}" "$1" "root://eosuser.cern.ch/{2}$2" "$3" "$4" "$5" "$6"\n'.format(
                self.analyzer_path, self.nameofExe, self.path_output))

        subprocess.call(["chmod", "755", self.script_name])


    def submit_job(self):
        # Submit the job to Condor
        print("Submitting job for sample: {0}".format(self.sample_name))
        subprocess.call(["condor_submit", self.condor_submit_name])

    def setup_and_submit_job(self):
        # Run all steps to set up and submit the job
        print("\nSetting up job for sample: {0}".format(self.sample_name))
        self.generate_argument_list()
        self.create_executable_script()
        self.write_condor_submission_file()
        self.submit_job()


def main():
    try:
        # Initialize job manager
        job_manager = CondorJobManager()
        # Job processing happens within the class
    except Exception as e:
        print(e)
        sys.exit(1)

if __name__ == "__main__":
    main()
