#!/bin/bash

###########################################################################
# NAME:   genSamples.sh
# AUTHOR: André Maia
# DATE:   10.06.2020
# DESC:   A script to generate sCO2 case samples using OpenFOAM for
#         varying Myong-Kasagi model coefficients. The coefficients
#         at the moment are: A*, Cd
#
#         Note: Bash's basic calculator "bc" was used since native bash
#         scripts do not support floating point calculations.
#
# USAGE:  genSamples.sh [CSV] [CASE] [YPLUSM] [PRTM] [PRT]
#         Where:
#         CSV     is a .csv file containing a description of the
#                 coefficients to be sampled. Syntax below.
#         CASE    is a folder containing the OpenFOAM code template
#                 for the case to be sampled.
#         YPLUSM  is the y+ model to be used:
#                 0=normal, 1=semi-local-scaling, 2=Wallin&Johannson
#         PRTM    is the Prandtl model to be use:
#                 0=constant, 1=Bae, 2=oos, 3=Tien, 4=TWL
#         PRT     is the turbulent Prandtl number.
#
#         Example: ./genSamples.sh coeffs.csv 42F 2 4 0.85
#
#         The syntax of the CSV file is [a b c d e f]. 
#         Where:
#         a is the name of the coefficient.
#         b is the initial value.
#         c is the final value.
#         d is the step.
#         e is to toggle activation. (0=inactive; 1=active)
#         f is a placeholder for the simulation loop. 
#           Must be the same as 'a'. 
###########################################################################

#
# Save the time of the simulation start.
#
t0=$(date +%s)

#
# Check the command line arguments.
#
if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ] || [ -z $5 ]; then

cat << EOF

DESC:  A script to generate sCO2 case samples using OpenFOAM for
       varying Myong-Kasagi model coefficients.

USAGE:  genSamples.sh [CSV] [CASE] [YPLUSM] [PRTM] [PRT]

        Where:
        CSV     is a .csv file containing a description of the
                coefficients to be sampled. Syntax below.
        CASE    is a folder containing the OpenFOAM code template
                for the case to be sampled.
        YPLUSM  is the y+ model to be used:
                0=normal, 1=semi-local-scaling, 2=Wallin&Johannson
        PRTM    is the Prandtl model to be use:
                0=constant, 1=Bae, 2=oos, 3=Tien, 4=TWL
        PRT     is the turbulent Prandtl number.

        Example: ./genSamples.sh coeffs.csv 42F 2 4 0.85

        The syntax of the CSV file is [a b c d e f]. Where:
        a is the name of the coefficient
        b is the initial value
        c is the final value
        d is the step
        e is to toggle activation (0=inactive; 1=active)
        f is a placeholder for the simulation loop. 
          Must be the same as 'a'. 

EOF

   exit 1
fi

#
# Define the model's static coefficients.
#
yPlus_Model=$3;
Prt_Model=$4;
Prt=$5;

#
# Define the names of the folder and path containing the OpenFOAM code for
# the desired case, and the locations of the RASProperties' file and where the 
# sample results are saved.
#
case_folder=$2
bk_path="$(pwd)/${case_folder}"
case_path="${bk_path}_tmp"
RAS_path="${case_path}/constant/RASProperties"
dict_path="${case_path}/postProcessing/sampleDict"

echo
echo "Simulation case: ${case_folder}"
echo -n "Extracting coefficients from file ${1}. "

#
# Read from the CSV file the coefficients to be sampled.
# 
while IFS=, read -ra line
do
   for x in ${line[@]}
   do
      coeff_list=(${coeff_list[@]} $x)
   done

   coeff_count=$((coeff_count+1))
done < $1

i=$((${#coeff_list[@]}/coeff_count)) # Number of params in every coefficient array

echo "Number of coefficients: ${coeff_count}."

#
# Read the contents of the RASProperties' file.
#
RAS_file=$(<"${bk_path}/constant/RASProperties")

if [ -z "$RAS_file" ]; then
   exit 1
fi

#
# Define the name of the simulation folder.
#
sim_folder=""

for (( x=0; x<=coeff_count-1; x++ )); do
   if (( $((coeff_list[x*i+4])) == 1 )); then 
      sim_folder="${sim_folder}_${coeff_list[x*i+0]}"
   fi
done

if [ -z "$sim_folder" ]; then
   echo "Please select at least one coefficient to be sampled."
   exit 1
else
   sim_folder="_${case_folder}_${yPlus_Model}_${Prt_Model}_Prt=${Prt}${sim_folder}"
   sim_path="$(pwd)/${sim_folder}"

   echo "Creating simulation folder ${sim_folder} at $(pwd)."
   mkdir -p "${sim_path}"
fi

echo -e "Starting simulation. Details of the simulation can be found at log_tmp.\n"

#
# Loop through the coefficient values.
#
for (( x=coeff_count-1,skip=0; x>=0; x-- )); do
   # Skip the coefficient if signaled.
   if (( $((coeff_list[x*i+4])) == 0 )); then 
      continue
   else
      # Check if the value should not be skept.
      if (( $((skip)) == 0 )); then
         #
         # Reset the files in the case folder back to original state.
         # 
         rm -rf "${case_path}"
         cp -rf "${bk_path}" "${case_path}"

         #
         # Define the name of the sample folder, log info and code to be added
         # to the RASProperties file.
         #
         sample_folder=""
         log=""
         nl=$'\n'
         code="myMyongKasagiKECoeffs${nl}{"

         code="${code}${nl}   Prt ${Prt};"
         code="${code}${nl}   yPlus_Model ${yPlus_Model};"
         code="${code}${nl}   Prt_Model ${Prt_Model};"

         for (( y=0; y<=coeff_count-1; y++ )); do
            if (( $((coeff_list[y*i+4])) == 1 )); then 
               sample_folder="${sample_folder}_${coeff_list[y*i+0]}=${coeff_list[y*i+5]}"
            fi

            log="${log}${coeff_list[y*i+0]}=${coeff_list[y*i+5]}"

            if (( y != coeff_count-1)); then 
               log="${log}, " 
            fi

            code="${code}${nl}   ${coeff_list[y*i+0]} ${coeff_list[y*i+5]};"
         done

         sample_path="$(pwd)/${sim_folder}/${sample_folder}"
         mkdir -p "${sample_path}"

         code="${code}${nl}}"

         echo -n "${log}. "
   
         #
         # Write the RASProperties file with the current coefficients.
         #      
         echo "${code}" >> $RAS_path

         #
         # Run the simulation for the given coefficients.
         #  
         echo -n "[PENDING]"
         mybuoyantPimpleFoam -case "${case_path}" > log_tmp
         echo -en "\b\b\b\b\b\b\b\b\b" # Remove the 'Pending' message status.

         #
         # Generate the desired results from the sampleDict dictionary
         # and extract the latest time of the simulation.
         #  
         post_log=$(postProcess -case "${case_path}" -func sampleDict -latestTime)
         latestTime=$(echo "${post_log}" | sed -n 's/Time = //p')

         #
         # Check if the simulation completed successfully.
         #  
         if (( $(echo "${latestTime} > 1.7" | bc -l) )); then
            # If yes, copy the result files to the sample folder.
	    mv "${dict_path}/${latestTime}" "${dict_path}/${sample_folder}"

	    # Obtain the error value of the current sample.
            error=$(python "../Evaluation/evalSamples.py" "${case_folder}" "$(pwd)/params.csv" 0 0 2>&1)

            cp -rf "${dict_path}/${sample_folder}/." "$sample_path"
            echo "[SUCCESS: e=${error}, t=${latestTime}]"
         else
            # If not, copy the simulation log to the sample folder.
            cp -rf "${case_path}/log" "$sample_path"
            echo "[ERROR: Check sample log]"
         fi
      fi

      # Reset the skip trigger.
      skip=0

      # Increase the current value.
      coeff_list[x*i+5]=$(echo "${coeff_list[x*i+5]} + ${coeff_list[x*i+3]}" | bc -l)

      # Check if the current value is greater than the final value.
      if (( $(echo "${coeff_list[x*i+5]} > ${coeff_list[x*i+2]}" | bc -l) )); then
         # Reset the current value back to the initial value.
         coeff_list[x*i+5]=$(echo "${coeff_list[x*i+1]}" | bc -l)

         # Signal that the next value should be skept, since it is a duplicate.
         skip=1
      else
         # Go back to beginning of the coefficient list.
         x=coeff_count
      fi
   fi
done

#
# Calculate the duration of the simulation.
#
tf=$(date +%s)
t=$(date -ud @$(( tf - t0 )) +"$(( (tf-t0)/(3600*24) )) days %H hours %M minutes %S seconds")

echo -e "\nSimulation completed in ${t}."
