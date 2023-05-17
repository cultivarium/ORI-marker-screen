import argparse
import glob
import os
import subprocess
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from unidecode import unidecode


def run(fastq_directory, mapping_file):

    print("Making database...")
    barcodes = []
    f = open("barcode_references.fasta", "w+")
    d = pd.read_csv("./barcodes.csv")
    for index, row in d.iterrows():
        barcode = unidecode(row["Corresponding ORI"].replace(" ", "_"))
        f.write(">" + barcode + "\n")
        f.write(row["ORI-barcode sequence"] + "\n")
        barcodes.append(barcode)
    f.close()

    results = defaultdict(list)
    stats = []

    samples = pd.read_csv(mapping_file)
    for index, row in samples.iterrows():
        filename = row["FileName"]
        print(filename)
        fn1 = glob.glob(fastq_directory + "/" + filename + "*_R1_*.fastq.gz")

        if len(fn1) == 0:  # file not found
            print("File not found: " + row['Sample'])
            barcode_hits = defaultdict(int)
            results["Sample"].append(row["Sample"])
            for barcode in barcodes:
                results[barcode].append(barcode_hits[barcode])

            stats_sample = {
                "Sample": row["Sample"],
                "Total_Reads": 0,
                "Quality_trimmed": 0,
                "Adapter_trimmed": 0,
                "Trimmed_Reads": 0,
                "Merged": 0,
                "Matched": 0,
            }
            stats.append(stats_sample)
            continue
        fn1 = fn1[0]
        fn2 = fn1.replace("_R1_", "_R2_")
        base = "/".join(fn1.split("/")[:-1]) + "/" + fn1.split("/")[-1].split("_")[0]

        # Run BBduk
        command = (
            BBTOOLS
            + "/bbduk.sh in1={f1} in2={f2} out1={f3} out2={f4} ref=./adapters.fa ktrim=r k=21 qtrim=r trimq=20 maq=20 minlen=50 entropy=0.3 threads=12"
        )
        command = command.format(
            f1=fn1,
            f2=fn2,
            f3=base + ".1.clean.fq.gz",
            f4=base + ".2.clean.fq.gz",
            f5=base,
        )

        output = subprocess.check_output(
            command, shell=True, stderr=subprocess.STDOUT
        ).decode()
        for line in output.split("\n"):
            if line.startswith("Input:"):
                reads = int(line.split(" reads")[0].split()[-1]) / 2
                bases = int(line.split(" bases")[0].split()[-1])
            if line.startswith("QTrimmed"):
                qtrimmed = int(line.split(" reads")[0].split()[-1]) / 2
            elif line.startswith("KTrimmed"):
                ktrimmed = int(line.split(" reads")[0].split()[-1]) / 2
            elif line.startswith("Result:"):
                result = int(line.split(" reads")[0].split()[-1]) / 2

        fn1 = base + ".1.clean.fq.gz"
        fn2 = fn1.replace(".1.clean.fq.gz", ".2.clean.fq.gz")

        ## Run BBmerge
        command = BBTOOLS + "/bbmerge.sh in1={f1} in2={f2} out={merged}.fasta"
        command = command.format(f1=fn1, f2=fn2, merged=base)

        output = subprocess.check_output(
            command, shell=True, stderr=subprocess.STDOUT
        ).decode()
        for line in output.split("\n"):
            if line.startswith("Joined"):
                merged = int(line.split()[1])

        ## Run VSEARCH
        command = "vsearch --search_exact {base}.fasta --db barcode_references.fasta --blast6out {base}.blast"
        command = command.format(base=base)
        output = subprocess.check_output(
            command, shell=True, stderr=subprocess.STDOUT
        ).decode()
        for line in output.split("\n"):
            if line.startswith("Matching"):
                matched = int(line.split(":")[1].strip().split()[0])

        ## Get hits of each type
        barcode_hits = defaultdict(int)

        f = open(base + ".blast")
        for line in f.readlines():
            barcode = line.split("\t")[1]
            barcode_hits[barcode] += 1

        results["Sample"].append(row["Sample"])
        for barcode in barcodes:
            results[barcode].append(barcode_hits[barcode])

        stats_sample = {
            "Sample": row["Sample"],
            "Total_Reads": reads,
            "Quality_trimmed": qtrimmed,
            "Adapter_trimmed": ktrimmed,
            "Trimmed_Reads": result,
            "Merged": merged,
            "Matched": matched,
        }
        stats.append(stats_sample)

    results = pd.DataFrame(results)
    stats = pd.DataFrame(stats)

    stats.to_csv("barcode_stats.tsv", sep="\t", index=False)
    results.to_csv("barcode_results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Whole plasmid sequencing of a plasmid ORI pool."
    )
    parser.add_argument(
        "-d",
        "--fastq_directory",
        action="store",
        default=None,
        required=True,
        help="Directory of FASTQ files. File names must take the form: sample_*_R1_*.fastq.gz",
    )
    parser.add_argument(
        "-m",
        "--mapping_file",
        action="store",
        default=None,
        required=True,
        help="Mapping file of comma separated columns Sample,Pool.",
    )
    parser.add_argument(
        "-b",
        "--bbmap_folder",
        action="store",
        default="~/bbmap/",
        required=False,
        help="Directory containing BBTools on your system",
    )

    args = parser.parse_args()

    BBTOOLS = args.bbmap_folder.rstrip("/")
    mapping = pd.read_csv(args.mapping_file)

    run(args.fastq_directory.rstrip("/"), args.mapping_file)
