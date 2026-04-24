#!/bin/bash

set -euo pipefail

# setup paths and working director
BASE="path to working directory"
script_path="$BASE/scripts/dbcan-scripts"
dbcan_path="$BASE/results/dbcan"
bakta_path="$BASE/results/bakta"
mkdir -p "$script_path" "$dbcan_path"
samples=$(cat "path to a list of stool samples")

# Run dbCAN per sample
for sample in ${samples}; do
    script_name="$script_path/${sample}.sh"
    log_name="$script_path/${sample}.log"
    err_name="$script_path/${sample}.err"

    echo "Processing $sample"

    cat <<EOF > "$script_name"
#!/bin/bash
#PBS -W group_list=dtu_00032 -A dtu_00032
#PBS -N ${sample}_stool_dbcan
#PBS -e ${err_name}
#PBS -o ${log_name}
#PBS -l nodes=1:thinnode:ppn=40
#PBS -l mem=188gb
#PBS -l walltime=12:00:00

set -euo pipefail

# Load modules
module load tools
module load dbcan/4.1.4

# Create a result direcotry =
mkdir -p "$dbcan_path/$sample"

# Run dbCAN
run_dbcan \\
    "$bakta_path/$sample/$sample.faa" protein \\
    -c "$bakta_path/$sample/$sample.gff3" \\
    --tools all \\
    --dia_cpu 40 \\
    --hmm_cpu 40 \\
    --tf_cpu 40 \\
    --stp_cpu 40 \\
    --cgc_substrate \\
    --db_dir "path_to/db/dbcan_db_09072024" \\
    --out_dir "$dbcan_path/$sample"
EOF

    chmod +x "$script_name"
    qsub "$script_name"
done
