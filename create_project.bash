#!/bin/bash

##########################################################################
#
# Platform: NCI Gadi HPC
# Usage: bash create_project.bash **run from Fastq-to-BAM base directory**
# Version: 2.0
#
# For more details see: https://github.com/Sydney-Informatics-Hub/Fastq-to-BAM
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement:
# The authors acknowledge the support provided by the Sydney Informatics Hub,
# a Core Research Facility of the University of Sydney. This research/project
# was undertaken with the assistance of resources and services from the National
# Computational Infrastructure (NCI), which is supported by the Australian
# Government, and the Australian BioCommons which is enabled by NCRIS via
# Bioplatforms Australia funding.
#
##########################################################################

#### functions####
function storage {
echo Do you require read/write to any Gadi storage other than /scratch/${project}? If yes, please enter all separated by space \[enter for no\]:

read more_storage
IFS=' ' read -r -a array <<< "$more_storage"

lstorage=scratch/${project}
for i in "${!array[@]}"
do
    path=$(echo ${array[i]} | sed 's/^\///g')
    lstorage+="+${path}"
done

echo
echo PBS lstorage directive will be: $lstorage
echo Is this correct? Enter y or n
read answer

if [[ $answer != y ]]
then
        storage
else
        echo Using lstorage $lstorage
        echo
        return 0
fi

}
###################

#  'scripts' directory
echo Enter the name of your scripts directory \(relative or absolute\):
read scripts
echo Updating scripts in $scripts
echo

# Make required starting directories:
echo Making Logs directory
mkdir -p $scripts/Logs
echo

echo Making Inputs directory
mkdir -p $scripts/Inputs
echo

# Config file
echo Enter your cohort name / the basename of your config file:
read cohort
config=${cohort}.config

if [ -f $config ]
then
        echo Using config $config

# if you would like to fill in the config file name for all scripts unhash the
# following sed commands and run all make_input.sh scripts without including
# the config name in the command

        #sed -i "s/^cohort=.*/cohort=..\/${cohort}/" ${scripts}/*.sh
        #sed -i "s/^cohort=.*/cohort=..\/${cohort}/" ${scripts}/*.pbs
        #sed -i "s/^config=.*/config=..\/${config}/" ${scripts}/*.sh
        #sed -i "s/^config=.*/config=..\/${config}/" ${scripts}/*.pbs

        echo
else
        echo $config does not exist - please fix. Aborting.
        exit
fi

# NCI project
echo Enter the name of your NCI project:
read project

echo Using NCI project $project for accounting and /scratch/${project} for read/write.
sed -i "s/#PBS -P.*/#PBS -P ${project}/" ${scripts}/*.pbs
echo

# Call storage function as many times as needed
storage
sed -i "s|#PBS -l storage=.*|#PBS -l storage=${lstorage}|" ${scripts}/*.pbs

# Reference genome
echo This directory needs a ./Reference or a symlink to your full \"Reference\" named directory.

if [ -d ./Reference ]
        then echo ./Reference already exists. Is this the complete and correct directory [y/n]?
        read answer

        if [[ $answer = y ]]
        then
                if [[ -L "./Reference" ]]; then
                        echo ./Reference is already symlinked. Proceeding with set up...
                        echo

                        echo Enter the name of your reference genome sequence \(include suffix\):
                        read ref

                        ref=Reference/${ref}
                        dict=${ref/\.[a-zA-Z]*/.dict}

                else
                        ln -s $refpath ./Reference

                        echo Enter the name of your reference genome sequence \(include suffix\):
                        read ref

                        ref=Reference/${ref}
                        dict=${ref/\.[a-zA-Z]*/.dict}

                fi

        elif [[ $answer = n ]]
        then

                echo Please remove or rename Reference directory currently stored at ./Reference and rerun this script. Aborting.
                exit
        fi

else
        echo Enter the full path to your reference directory:
        read refpath
        echo Creating symlink $refpath to ./Reference ...
        ln -s $refpath Reference

        echo Enter the name of your reference genome sequence \(include suffix\):
        read ref

        ref=Reference/${ref}
        dict=${ref/\.[a-zA-Z]*/.dict}

fi

if [ ! -f ${ref} ]
then
        echo ${ref} does not exist - please check. Aborting.
        exit
elif [ ! -f ${dict} ]
then
        echo ${dict} does not exist - please check. Aborting.
         exit
else
        echo Using reference genome files ${ref} and ${dict} ...
        sed -i "s|^ref=.*|ref=../${ref}|" ${scripts}/*.sh
        sed -i "s|^ref=.*|ref=../${ref}|" ${scripts}/*.pbs
        sed -i "s|^dict=.*|dict=../${dict}|" ${scripts}/*.sh
        sed -i "s|^dict=.*|dict=../${dict}|" ${scripts}/*.pbs

fi

echo
echo All done. Ready to run the Fastq-to-BAM pipeline
echo
echo The scripts in $scripts directory have now been updated to include the following:
echo
printf "\tNCI accounting project: ${project}\n \
\tPBS lstorage directive: ${lstorage}\n \
\tCohort config file: ${cohort}.config\n \
\tReference genome sequence: ${ref}\n \
\tReference genome dictionary file: ${dict}\n"


