#!/bin/bash

# setup paths and working director
base_dir="/home/projects/dtu_00032/analysis/milk-cohort"
samples="a path to a list of samplest"
script_path="$base_dir/results/bins/bins_for_das/das_tool_scripts"
log_path="$base_dir/results/bins/bins_for_das/logs"
mkdir -p "$script_path" "$log_path"

# Generate script per sample
mapfile -t samples_array < "$samples"

for n in "${samples_array[@]}"; do

    i="$base_dir/results/bins/bins_for_das/${n}"

    script_name="${script_path}/run_das_tool_${n}.sh"
    log_name="${log_path}/${n}.log"
    err_name="${log_path}/${n}.err"

    {
        echo "#!/bin/bash"
        echo "#PBS -W group_list=dtu_00032 -A dtu_00032"
        echo "#PBS -N das-${n}"
        echo "#PBS -e ${err_name}"
        echo "#PBS -o ${log_name}"
        echo "#PBS -l nodes=1:thinnode:ppn=40"
        echo "#PBS -l mem=188gb"
        echo "#PBS -l walltime=24:00:00"
        echo ""
        echo "module load tools"
        echo "module load gcc/13.2.0"
        echo "module load intel/perflibs"
        echo "module load R/4.4.0"
        echo "module load ruby/3.2.0"
        echo "module load pullseq/1.0.2"
        echo "module load prodigal/2.6.3"
        echo "module load diamond/2.1.8"
        echo "module load perl"
        echo "module load ncbi-blast/2.15.0+"
        echo "module load das_tool/1.1.6"
        echo ""
        echo "cd \"$i\""
        echo ""
        echo "contig=\"$base_dir/results/megahit/final.contigs.fa\""
        echo ""
    } > "$script_name"

    input_files=""
    input_labels=""

    [[ -f "$i/no-header_contig_bins.tsv" ]] && \
        input_files="no-header_contig_bins.tsv" && \
        input_labels="contig_bins"

    [[ -f "$i/swapped_aae_y_clusters.tsv" ]] && \
        input_files="${input_files:+${input_files},}swapped_aae_y_clusters.tsv" && \
        input_labels="${input_labels:+${input_labels},}aae_y"

    [[ -f "$i/swapped_aae_z_clusters.tsv" ]] && \
        input_files="${input_files:+${input_files},}swapped_aae_z_clusters.tsv" && \
        input_labels="${input_labels:+${input_labels},}aae_z"

    [[ -f "$i/swapped_vae_clusters.tsv" ]] && \
        input_files="${input_files:+${input_files},}swapped_vae_clusters.tsv" && \
        input_labels="${input_labels:+${input_labels},}vae"

    if [[ -n "$input_files" ]]; then
        {
            echo "DAS_Tool -i ${input_files} \\"
            echo "         -l ${input_labels} \\"
            echo "         -c \${contig} \\"
            echo "         -o ${n} \\"
            echo "         -t 40 \\"
            echo "         --dbDirectory $base_dir/db/DAS_Tool_db_1.1.6"
        } >> "$script_name"
    else
        echo "echo \"No valid input files found for ${n}\"" >> "$script_name"
    fi

    chmod +x "$script_name"
    qsub "$script_name"

done
