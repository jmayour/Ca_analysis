#!/bin/bash
# Sample batchscript to run a simple MATLAB job on HPC
#SBATCH --partition=bch-compute	 			# queue to be used
#SBATCH --time=00:05:00	 			# Running time (in hours-minutes-seconds)
#SBATCH --job-name=test-compute 			# Job name
#SBATCH --mail-type=BEGIN,END,FAIL 		# send and email when the job begins, ends or fails
#SBATCH --mail-user=your_email_address	 	# Email address to send the job status
#SBATCH --output=output_%j.txt 			# Name of the output file
#SBATCH --nodes=1				# Number of gpu nodes
#SBATCH --ntasks=1				# Number of gpu devices on one gpu node

module load matlab
matlab -nodisplay -nosplash -nodesktop -r "addpath('../../Ca_analysis');run_analysis('Sample_data.xlsx', 'Sample_data_parameters.xlsx');exit;"
