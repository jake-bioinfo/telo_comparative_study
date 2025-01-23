#!/bin/bash

# This script will run in order telseq, telhunter, telogator and tecat 

# Input options
# Initialize the options
temp_directory=''
output_directory=''
input_dir=''
sample=''
threads=''
platform='ont'

# Set help info
help_message="Usage: compare.sh -t <temp_directory> -o <output_directory> -i <input_dir> -s <sample> -n <threads>
    -t, --temp_directory: The temporary directory to store intermediate files
    -o, --output_directory: The output directory to store the final results
    -i, --input_dir: The input directory containing the input files
    -s, --sample: The sample name
    -n, --threads: The number of threads to use
    -p, --platform: The platform used for sequencing (nanopore or pacbio)
"

# Check foir any args, print help if none supplied
if [[ $# -eq 0 ]]; then
    echo "$help_message"
    exit 1
fi

# Read the options
ARGS=$( getopt -o t:o:i:s:n:p: -l temp_directory:,output_directory:,input_dir:,sample:,threads:,platform: -n 'compare.sh' -- "$@" )

eval set -- "$ARGS"

while true; do
    case "$1" in
        -t|--temp_directory)
            temp_directory="$2"
            shift 2
            ;;
        -o|--output_directory)
            output_directory="$2"
            shift 2
            ;;
        -i|--input_dir)
            input_dir="$2"
            shift 2
            ;;
        -s|--sample)
            sample="$2"
            shift 2
            ;;
        -n|--threads)
            threads="$2"
            shift 2
            ;;
        -p|--platform)
            platform="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done

# Check if the required options are set
if [[ -z $temp_directory ]] || [[ -z $output_directory ]] || [[ -z $input_dir ]] || [[ -z $sample ]]; then
    echo "All options are required"
    echo "$help_message"
    exit 1
fi

# Print the options
readout_vars="
temp_directory: $temp_directory
output_directory: $output_directory
input_dir: $input_dir
sample: $sample
threads: $threads
platform: $platform
"
echo "$readout_vars"

# Create input directory for the sample if it doesn't already exist
mkdir -p $temp_directory/$sample/input

# Sync data to the input directory
rsync -av --progress $input_dir/ $temp_directory/$sample/input

# Short read read length
echo "Calculating short read length"
read_len=$(samtools view $temp_directory/$sample/input/*.bam | awk '{print length($10)}' | uniq -c | sort -nr | head -n 1 | awk '{print $2}')

# Run telseq
echo "Running telseq"

# Assign variables for running telseq
bam=$temp_directory/$sample/input/*.bam
out_dir=$output_directory/$sample/telseq
mkdir -p $out_dir

# Activate the conda environment
eval "$(conda shell.bash hook)"
conda activate telo_tools

# Determine time taken to run telseq
ts_start_time=$(date +%s)

# Echo start time
echo "Telseq start time: ${ts_start_time}"

# Run telseq
telseq $bam -r ${read_len} > $out_dir/telseq_results.txt

# Determine time taken to run telseq
ts_end_time=$(date +%s)

# Echo end time
echo "Telseq end time: ${ts_end_time}"

# Calculate time taken to run telseq
ts_time_taken=$((ts_end_time-ts_start_time))

# Echo time taken
echo "Telseq time taken: ${ts_time_taken}"

# Deactivate telo_tools
conda deactivate

# Run telomerehunter
echo "Running telomerehunter"

# Activate the telomerehunter environment
eval "$(conda shell.bash hook)"
conda activate telomerehunter

# telomerehunter start time
th_start_time=$(date +%s)

# Echo start time
echo "Telomerehunter start time: ${th_start_time}"

# Create telomerehunter output directory
mkdir -p $output_directory/$sample/telomerehunter

# Run telomerehunter
telomerehunter -ibc $bam -pl -o $output_directory/$sample/telomerehunter -p $sample

# telomerehunter end time
th_end_time=$(date +%s)

# Echo end time
echo "Telomerehunter end time: ${th_end_time}"

# Calculate time taken to run telomerehunter
th_time_taken=$((th_end_time-th_start_time))

# Echo time taken
echo "Telomerehunter time taken: ${th_time_taken}"

# Deactivate telomerehunter
conda deactivate

# Run telogator
echo "Running telogator"

# Activate the telogator environment
eval "$(conda shell.bash hook)"
conda activate telogator2

# Assign telogator bin
telogator_bin="/home/jake/science_projects/telo/other_tools/telogator2/"

# telogator start time
tg_start_time=$(date +%s)

# Echo start time
echo "Telogator start time: ${tg_start_time}"

# Create telogator output directory
mkdir -p $output_directory/$sample/telogator

# Run telogator
# if platform is pacbio, use hifi just for telogator
if [ $platform == "pb" ]; then
    plat="hifi"
else
    plat="ont"
fi

python $telogator_bin/telogator2.py -i $temp_directory/$sample/input/lr/*.fastq* -o $output_directory/$sample/telogator -p $threads -t $telogator_bin/resources/telogator-ref.fa.gz -r $plat --debug-nosubtel

# telogator end time
tg_end_time=$(date +%s)

# Echo end time
echo "Telogator end time: ${tg_end_time}"

# Calculate time taken to run telogator
tg_time_taken=$((tg_end_time-tg_start_time))

# Echo time taken
echo "Telogator time taken: ${tg_time_taken}"

# Deactivate telogator
conda deactivate

# Run tecat
echo "Running tecat"

# TECAT start time
tc_start_time=$(date +%s)

# Echo start time
echo "TECAT start time: ${tc_start_time}"

# Create tecat output directory
mkdir -p $output_directory/$sample/tecat

# Assign tecat input and reference
tecat_input=$(find $temp_directory/$sample/input/lr/ -regex ".*\.\(fq\|fastq\|fastq\.gz\)$" | head -n 1)
tecat_ref=$(find $temp_directory/$sample/input/ref/ -regex ".*\.\(fa\|fasta\|fasta\.gz\)$" | head -n 1)

# Run tecat
./tecat_pipeline.R -i $tecat_input -o $output_directory/$sample/tecat -t $threads -p $sample -r $output_directory/$sample/tecat/tecat_results.csv -f $platform -x $tecat_ref

# TECAT end time
tc_end_time=$(date +%s)

# Echo end time
echo "TECAT end time: ${tc_end_time}"

# Calculate time taken to run tecat
tc_time_taken=$((tc_end_time-tc_start_time))

# Echo time taken
echo "TECAT time taken: ${tc_time_taken}"

# Print time taken for each tool
echo "Telseq time taken: ${ts_time_taken}"
echo "Telomerehunter time taken: ${th_time_taken}"
echo "Telogator time taken: ${tg_time_taken}"
echo "TECAT time taken: ${tc_time_taken}"

# Print the output directory
echo "Output directory: $output_directory/$sample"

# Print the end time
echo "End time: $(date)"

# Clean up the temp directory
rm -rf $temp_directory/$sample

# Echo ...
echo "JOBS DONE... :)"
