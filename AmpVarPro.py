import os
import subprocess
import sys
import argparse
from multiprocessing import Pool
import shutil
import glob

def process_barcode(barcode, args):
    try:
        barcode_str = f"{barcode:02d}"
        barcode_dir = f"barcode{barcode_str}"
        input_fastq_pattern = os.path.join(args.input_directory, barcode_dir, "*.fastq")
        input_fastq_files = glob.glob(input_fastq_pattern)

        if not input_fastq_files:
            print(f"No fastq files found for barcode {barcode_str}. Skipping.")
            return

        input_fastq = input_fastq_files[0]
        output_directory = barcode_dir
        os.makedirs(output_directory, exist_ok=True)
        output_fastq = os.path.join(output_directory, f"{barcode_str}.fastq")
        output_filt_fastq = os.path.join(output_directory, f"{barcode_str}_filt.fastq")
        output_sorted_bam = os.path.join(output_directory, f"{barcode_str}.sorted.bam")
        output_final_vcf = os.path.join(output_directory, f"{barcode_str}.final.vcf")
        output_ann_vcf = os.path.join(output_directory, f"{barcode_str}.ann.vcf")
        output_ann_final_vcf = os.path.join(output_directory, f"{barcode_str}.ann.final.vcf")

        print(f"Processing barcode {barcode_str}")

        # Concatenating fastq files
        print("Concatenating fastq files")
        subprocess.run(["cp", input_fastq, output_fastq], check=True)

        # Read filtering
        print("Read filtering")
        subprocess.run(["fastp", "-i", output_fastq, "-o", output_filt_fastq, "-q", args.phred_quality, "--length_required", str(args.min_length), "--length_limit", str(args.max_length)], check=True)

        # Mapping to reference genome
        print("Mapping to reference genome")
        minimap_cmd = f"minimap2 -ax map-ont -t {args.threads} {args.reference_genome} {output_filt_fastq} | samtools view -bS -F 4 | samtools sort -o {output_sorted_bam} -T reads.tmp"
        subprocess.run(minimap_cmd, shell=True, check=True)

        # Indexing mapped file
        print("Indexing mapped file")
        subprocess.run(["samtools", "index", output_sorted_bam], check=True)

        # Variant calling
        print("Variant calling")
        clair3_output_dir = output_directory  # No additional subfolder for Clair3 output
        os.makedirs(clair3_output_dir, exist_ok=True)
        subprocess.run(["run_clair3.sh", "-b", output_sorted_bam, "-f", args.reference_genome, "-m", args.model, "-t", str(args.threads), "-p", "ont", "-o", clair3_output_dir, "--haploid_precise", "--include_all_ctgs", "--no_phasing_for_fa"], check=True)

        # Gunzip and move VCF file
        print("Gunzipping and moving VCF file")
        subprocess.run(["gunzip", os.path.join(clair3_output_dir, "merge_output.vcf.gz")], check=True)
        subprocess.run(["mv", os.path.join(clair3_output_dir, "merge_output.vcf"), output_final_vcf], check=True)

        # Variant annotation
        print("Variant annotation")
        process = subprocess.run(["snpEff", "-v", args.snpeff_data, output_final_vcf], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        with open(output_ann_vcf, 'w') as output_file:
            output_file.write(process.stdout.decode('utf-8'))

        # Filter lines starting with '##' in the annotated VCF file
        print("Filtering VCF file")
        with open(output_ann_vcf, 'r') as input_file, open(output_ann_final_vcf, 'w') as output_file:
            for line in input_file:
                if not line.startswith('##'):
                    output_file.write(line)

        print(f"Processing barcode {barcode_str} completed")

    except subprocess.CalledProcessError as e:
        print(f"Error processing barcode {barcode_str}: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error processing barcode {barcode_str}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Barcode Processing Script")
    parser.add_argument("input_directory", type=str, help="Input directory containing barcode subdirectories")
    parser.add_argument("min_length", type=int, help="Minimum read length")
    parser.add_argument("max_length", type=int, help="Maximum read length")
    parser.add_argument("threads", type=int, help="Number of threads")
    parser.add_argument("reference_genome", type=str, help="Path to the reference genome")
    parser.add_argument("phred_quality", type=str, help="Phred quality")
    parser.add_argument("snpeff_data", type=str, help="Path to SNPeff data")
    parser.add_argument("model", type=str, help="Path to the model")
    args = parser.parse_args()

    # Define the number of processes (adjust as needed)
    num_processes = 4

    # Create a Pool for parallel processing
    with Pool(num_processes) as pool:
        # Iterate over barcodes in parallel
        pool.starmap(process_barcode, [(barcode, args) for barcode in range(1, 97)])

    print("The End")

