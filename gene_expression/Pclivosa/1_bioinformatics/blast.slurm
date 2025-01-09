#!/bin/bash
#SBATCH -J blast
#SBATCH -n 160
#SBATCH -N 10
#SBATCH -p shortq7
#SBATCH -o blast.o%j
#SBATCH -e blast.e%j
#SBATCH -t 6:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=studivanms@gmail.com

module load launcher
module load python-3.7.4-gcc-8.3.0-3tniqr5
module load nmap-7.70-gcc-9.2.0-zzeqwvb




export EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher 
export CONTROL_FILE=bl
export LAUNCHER_WORKDIR=`pwd`
# 
# Variable description:
#
#  EXECUTABLE     = full path to the job launcher executable
#  CONTROL_FILE   = text input file which specifies
#                   executable for each process
#                   (should be located in LAUNCHER_WORKDIR)
#  LAUNCHER_WORKDIR        = location of working directory
#
#      <------ End Setup Parameters ------>
#--------------------------------------------------------
#--------------------------------------------------------

#----------------
# Error Checking
#----------------

if [ ! -e $LAUNCHER_WORKDIR ]; then
    echo " "
    echo "Error: unable to change to working directory."
	echo "       $LAUNCHER_WORKDIR"
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -f $EXECUTABLE ]; then
	echo " "
	echo "Error: unable to find launcher executable $EXECUTABLE."
	echo " "
	echo "Job not submitted."
	exit
fi

if [ ! -f $LAUNCHER_WORKDIR/$CONTROL_FILE ]; then
	echo " "
	echo "Error: unable to find input control file $CONTROL_FILE."
	echo " "
	echo "Job not submitted."
	exit
fi


#----------------
# Job Submission
#----------------

cd $LAUNCHER_WORKDIR/
echo " WORKING DIR:   $LAUNCHER_WORKDIR/"

$TACC_LAUNCHER_DIR/paramrun $EXECUTABLE $CONTROL_FILE

echo " "
echo " Parameteric Job Complete"
echo " "
