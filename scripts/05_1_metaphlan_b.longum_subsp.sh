#!/bin/sh

# Setup working directory and paths
base_dir="path to working directory"
scriptpath="${base_dir}/results/bif-metaphlan/job_scripts/mp"
logpath="${base_dir}/results/bif-metaphlan/logs/mp"
mkdir -p "$scriptpath" "$logpath"

########################################
# 1. Stool samples (host-removed input)
########################################
stool_samples=$(cat "path to a list of stool samples")

for N in $stool_samples; do
    script_name="${scriptpath}/${N}_mp.sh"
    log_name="${logpath}/${N}_mp.log"
    err_name="${logpath}/${N}_mp.err"

    cat <<EOL > "$script_name"
#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N ${N}-MP-stool
#PBS -e ${err_name}
#PBS -o ${log_name}
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=8:00:00

# Load modules
module load tools
module load metaphlan/4.1.0

# Setup working directory and paths
tf="$base_dir/results/host_removed/${N}_1.fq.gz"
tr="$base_dir/results/host_removed/${N}_2.fq.gz"
mkdir -p $base_dir/results/bif-metaphlan/{bowtie2,sams,original,convGTDB}
d="path_to_db/metaphlan_bif_db_202307"

# Run Metaphlan for stool samples
if [ -f "\$tf" ] && [ -f "\$tr" ]; then
    sam="$base_dir/results/bif-metaphlan/sams/${N}.sam"
    bowtie="$base_dir/results/bif-metaphlan/bowtie2/${N}.bowtie2.bz2"
    mpfile="$base_dir/results/bif-metaphlan/original/${N}_profiled.txt"
    convfile="$base_dir/results/bif-metaphlan/convGTDB/${N}_conv_profiled.txt"

    echo "Running Metaphlan (stool) with ${N}"

    metaphlan "\${tf},\${tr}" \
        --unclassified_estimation \
        --nproc 40 \
        --input_type fastq \
        --bowtie2db "\$d" \
        -x "mpa_vJun23_CHOCOPhlAnSGB_202403_lon_subsp" \
        -s "\$sam" \
        --bowtie2out "\$bowtie" \
        -o "\$mpfile"
    
    # convert to GTDB taxa
    bash ${base_dir}/scripts/sgb_transfer.sh \
        -i "\$mpfile" \
        -o "\$convfile" \
        -d "path_to_db/mpa_vJun23_CHOCOPhlAnSGB_202307_SGB2GTDB.tsv"
else
    echo "Error: Input files for sample ${N} not found"
    exit 1
fi
EOL

    chmod +x "$script_name"
    qsub "$script_name"
done


########################################
# 2. Enriched samples (trimmed input)
########################################
enriched_samples=$(cat "path to a list of enriched samples")

for N in $enriched_samples; do
    script_name="${scriptpath}/${N}_mp.sh"
    log_name="${logpath}/${N}_mp.log"
    err_name="${logpath}/${N}_mp.err"

    cat <<EOL > "$script_name"
#!/bin/sh
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N ${N}-MP-enriched
#PBS -e ${err_name}
#PBS -o ${log_name}
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=24:00:00

# Load modules
module load tools
module load metaphlan/4.1.0

# Setup working directory and paths
tf="$base_dir/results/trimmomatic/${N}_trimmed_1.fq.gz"
tr="$base_dir/results/trimmomatic/${N}_trimmed_2.fq.gz"
mkdir -p $base_dir/results/bif-metaphlan/{bowtie2,sams,original,convGTDB}
d="path_to_db/metaphlan_bif_db_202307"

# Run Metaphlan for enriched samples
if [ -f "\$tf" ] && [ -f "\$tr" ]; then
    sam="$base_dir/results/bif-metaphlan/sams/${N}.sam"
    bowtie="$base_dir/results/bif-metaphlan/bowtie2/${N}.bowtie2.bz2"
    mpfile="$base_dir/results/bif-metaphlan/original/${N}_profiled.txt"
    convfile="$base_dir/results/bif-metaphlan/convGTDB/${N}_conv_profiled.txt"
    
    echo "Running Metaphlan (enriched) with ${N}"

    metaphlan "\${tf},\${tr}" \
        --nproc 40 \
        --input_type fastq \
        --bowtie2db "\$d" \
        -x "mpa_vJun23_CHOCOPhlAnSGB_202403_lon_subsp" \
        -s "\$sam" \
        --bowtie2out "\$bowtie" \
        -o "\$mpfile"
    
    # convert to GTDB taxa
    bash ${base_dir}/scripts/sgb_transfer.sh \
        -i "\$mpfile" \
        -o "\$convfile" \
        -d "path_to_db/mpa_vJun23_CHOCOPhlAnSGB_202307_SGB2GTDB.tsv"
else
    echo "Error: Input files for sample ${N} not found"
    exit 1

fi
EOL

    chmod +x "$script_name"
    qsub "$script_name"
done
