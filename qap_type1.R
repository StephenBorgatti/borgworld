#!/bin/bash

#SBATCH --job-name qaptest								# job name
#SBATCH --mail-user sborg2@uky.edu   					# email
#SBATCH --mail-type ALL									# which mails do you want? (set up a rule on your outlook)
#SBATCH --ntasks=3 										# number of cores
#SBATCH --nodes=1 										# number of nodes - leave at 1
#SBATCH --time=0-05:00:00   							# if you don't know, set 14 days, extra space for jobs with 3 day max                    
#SBATCH --partition=normal     						    # extra partitions for faster jobs or graphic jobs                      
#SBATCH --account=cea_rkr235_f25mgt795						# use your login name
#SBATCH --output=/home/sborg2/nov17/output/homo.txt  # HIGHLY recommend you put those in a folder
#SBATCH --error=/home/sborg2/nov17/error/homo.txt  	# HIGHLY recommend you put those in a folder

#R Program execution command
cd /home/sborg2/R 									# set your working directory

module add r-4.4.0-gcc-9.3.0-7nl5q44
					

# Print this sub-job's task ID
Rscript qap_type1_error_test_homophily.R	# run the Rscript 