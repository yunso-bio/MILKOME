#!/bin/sh

cd /home/projects/dtu_00032/analysis/milk-cohort/bins

# Create directories for batch scripts and sample lists
mkdir -p checkm-scripts
mkdir -p sample-lists

# List all samples in new-all-bins
samples=($(ls new-all-bins/*))
total_samples=${#samples[@]}
batch_size=$((total_samples / 300 + 1))

# Create 300 scripts
for (( batch=0; batch<300; batch++ )); do
    start_index=$((batch * batch_size))
    end_index=$((start_index + batch_size))
    if [ $end_index -gt $total_samples ]; then
        end_index=$total_samples
    fi

    # Create a file with the list of samples for this batch only if it contains samples
    sample_list="sample-lists/sample_list_${batch}.txt"
    if [ $start_index -lt $total_samples ]; then
        for (( i=start_index; i<end_index; i++ )); do
            echo "${samples[$i]}" >> $sample_list
        done
    fi

    script_name="checkm-scripts/run_batch_${batch}.sh"

    echo "#!/bin/sh" > $script_name
    echo "#PBS -W group_list=dtu_00032 -A dtu_00032" >> $script_name
    echo "#PBS -N checkm2-${batch}" >> $script_name
    echo "#PBS -e checkm2-${batch}.err" >> $script_name
    echo "#PBS -o checkm2-${batch}.log" >> $script_name
    echo "#PBS -l nodes=1:ppn=40" >> $script_name
    echo "#PBS -l mem=188gb" >> $script_name
    echo "#PBS -l walltime=48:00:00" >> $script_name
    echo "" >> $script_name
    echo "module load tools" >> $script_name
    echo "module load checkm2/1.0.2" >> $script_name
    echo "" >> $script_name
    echo "db=\"/home/projects/dtu_00032/db/CheckM2_database/uniref100.KO.1.dmnd\"" >> $script_name
    echo "" >> $script_name
    echo "cd /home/projects/dtu_00032/analysis/milk-cohort/bins" >> $script_name
    echo "" >> $script_name
    echo "if [ -f ${sample_list} ]; then" >> $script_name
    echo "    echo \"Running CheckM2 for batch ${batch}\"" >> $script_name
    echo "    checkm2 predict --threads 40 \\" >> $script_name
    echo "                    --input \$(cat ${sample_list} | tr '\n' ' ') \\" >> $script_name
    echo "                    --output-directory checkm-allbins/batch_${batch} \\" >> $script_name
    echo "                    --database_path \$db ;" >> $script_name
    echo "fi" >> $script_name

    chmod +x $script_name
done
