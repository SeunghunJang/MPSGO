#!/bin/bash
#SBATCH --job-name=mol-rlx       # Job name

#SBATCH --nodelist=name_node     # set name of node
##SBATCH --exclusive             # Use of full Memory in one job at only one node

#SBATCH -p batch                 # set name of batch
#SBATCH --cpus-per-task=1        # Number of cores per MPI rank
#SBATCH --nodes=1                # Number of nodes
#SBATCH --ntasks-per-node=4      # How many tasks on each node
#SBATCH --output=log.out         # Standard output and error log

DIR=$PWD
echo $DIR

~/cal_tmp/10_orca/orca_5_0_4_linux_x86-64_shared_openmpi411/orca   $DIR/input.inp   >   std.out
