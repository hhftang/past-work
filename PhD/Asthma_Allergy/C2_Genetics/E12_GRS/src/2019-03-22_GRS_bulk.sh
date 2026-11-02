#!/bin/bash

#Partition: 'Short' (<24h) or 'Long' (>24h)
#SBATCH -p Short

#Give the job a name:
#SBATCH --job-name="classify"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=10240

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-23:00:00

# Error and output messages
#SBATCH --error=/projects/htang/CAS/C2_genetics/GRS/log/grs_create%J.err
#SBATCH --output=/projects/htang/CAS/C2_genetics/GRS/log/grs_create%J.out

# ------------------------------------------------------------------------------
# GENERAL SCRIPT TO RUN SHRINE ET AL GRS IN BULK
# ------------------------------------------------------------------------------

src="/projects/htang/CAS/C2_genetics/GRS/src"

"$src"/2019-03-22_GRS_code_spirometa_shrine.sh SpiroMeta_FEV1
"$src"/2019-03-22_GRS_code_spirometa_shrine.sh SpiroMeta_FVC
"$src"/2019-03-22_GRS_code_spirometa_shrine.sh SpiroMeta_PEF
"$src"/2019-03-22_GRS_code_spirometa_shrine.sh SpiroMeta_RATIO

"$src"/2019-03-22_GRS_code_ukb_shrine.sh UKB_FEV1
"$src"/2019-03-22_GRS_code_ukb_shrine.sh UKB_FVC
"$src"/2019-03-22_GRS_code_ukb_shrine.sh UKB_PEF
"$src"/2019-03-22_GRS_code_ukb_shrine.sh UKB_RATIO

"$src"/2019-03-22_GRS_code_ukb_shrine.sh Meta_FEV1
"$src"/2019-03-22_GRS_code_ukb_shrine.sh Meta_FVC
"$src"/2019-03-22_GRS_code_ukb_shrine.sh Meta_PEF
"$src"/2019-03-22_GRS_code_ukb_shrine.sh Meta_RATIO