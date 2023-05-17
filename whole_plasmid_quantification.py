import argparse
import glob
import os
import subprocess
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO


## Constants
MAX_MISMATCH = 1
MIN_MAPPING_LENGTH = 100
BBTOOLS = ""


def quality_trim(fastq_directory):
    print("Quality trimming")
    for fn1 in glob.glob(fastq_directory.rstrip("/") + "/*_R1_*.fastq.gz"):
        fn2 = fn1.replace("_R1_", "_R2_")
        base = "/".join(fn1.split("/")[:-1]) + "/" + fn1.split("/")[-1].split("_")[0]

        # Run BBduk
        if not os.path.isfile(base + ".1.clean.fq.gz"):
            command = (
                BBTOOLS.rstrip("/")
                + "/bbduk.sh in1={f1} in2={f2} out1={f3} out2={f4} ref=./adapters.fa ktrim=r k=21 qtrim=r trimq=20 maq=20 minlen=50 entropy=0.3 threads=12"
            )
            command = command.format(
                f1=fn1,
                f2=fn2,
                f3=base + ".1.clean.fq.gz",
                f4=base + ".2.clean.fq.gz",
                f5=base,
            )
            process = subprocess.Popen(command, shell=True)
            process.wait()


def make_index():
    print("Indexing")
    for fasta in glob.glob("./pool_data/*.fasta"):
        name = fasta.split(".fasta")[0]
        cmd = "bowtie2-build {fasta} {fasta}"
        cmd = cmd.format(fasta=fasta)
        print(cmd)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()


def map_reads(fastq_directory, mapping):
    print("Mapping")

    if not os.path.exists('./results_bowtie2/'):
        os.makedirs(path)

    already_sequenced = "|".join(glob.glob("./results_bowtie2/*.bam"))
    for index, row in mapping.iterrows():
        ## make index file
        print("Running " + row["Sample"])
        r1 = glob.glob(
            fastq_directory.rstrip("/") + "/" + row["Sample"] + "*.1.clean.fq.gz"
        )[0]
        r2 = glob.glob(
            fastq_directory.rstrip("/") + "/" + row["Sample"] + "*.2.clean.fq.gz"
        )[0]
        cmd = """bowtie2 -x {index} -1 {r1} -2 {r2} -p 12 | samtools view -bS > {output}.bam; 
        samtools sort -@ 6 -o {output}.sort.bam {output}.bam; rm {output}.bam; samtools index {output}.sort.bam"""
        cmd = cmd.format(
            index="pool_data/" + row["Pool"] + ".fasta",
            output="results_bowtie2/" + row["Sample"],
            r1=r1,
            r2=r2,
        )
        if row["Sample"] not in already_sequenced:
            process = subprocess.Popen(cmd, shell=True)
            process.wait()


def identify_origins(mapping):
    print("Reading BAMs...")
    origins = pd.read_csv("origins.tsv", sep="\t")
    origins = origins.groupby("Construct").first().reset_index()
    origins = origins.set_index("Construct")

    ## Pysam time

    # for each file
    results = defaultdict(list)

    for index, row in mapping.iterrows():
        samfile = pysam.AlignmentFile(
            "{output}.sort.bam".format(output="results_bowtie2/" + row["Sample"]), "rb"
        )
        fasta = "pool_data/" + row["Pool"] + ".fasta"
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id in origins.index:
                print(record.id)
                start_pos = origins.loc[record.id]["start"] - 1
                stop_pos = origins.loc[record.id]["stop"] - 1
                coverages = {}
                for i in range(start_pos, stop_pos + 1):
                    coverages[i] = 0
                for pileupcolumn in samfile.pileup(
                    record.id,
                    start_pos,
                    stop_pos,
                    truncate=False,
                    min_mapping_quality=20,
                    ignore_overlap=False,
                ):
                    for pileupread in pileupcolumn.pileups:
                        if (
                            not pileupread.is_del
                            and not pileupread.is_refskip
                            and float(pileupread.alignment.get_tag("NM"))
                            <= MAX_MISMATCH
                            and len(pileupread.alignment.get_reference_positions())
                            >= MIN_MAPPING_LENGTH
                        ):
                            if pileupcolumn.pos in coverages:
                                coverages[pileupcolumn.pos] += 1

                mean_coverage = str(round(np.mean(list(coverages.values())), 2))
                breadth = len([x for x in coverages.values() if x > 0]) / float(
                    stop_pos - start_pos
                )

                results["Sample"].append(row["Sample"])
                results["Backbone"].append("_".join(record.id.split("_")[:-1]))
                results["Origin"].append(origins.loc[record.id]["Origin"])
                results["Coverage"].append(mean_coverage)
                results["Breadth"].append(breadth)
            else:
                print("ERR: Cannot find origin for " + record.id)

    results = pd.DataFrame(results)
    pd.pivot(results, index="Sample", columns="Origin", values="Coverage").to_csv(
        "plasmid_library_coverage.tsv", sep="\t"
    )
    pd.pivot(results, index="Sample", columns="Origin", values="Breadth").to_csv(
        "plasmid_library_breadth.tsv", sep="\t"
    )


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

    BBTOOLS = args.bbmap_folder
    mapping = pd.read_csv(args.mapping_file)

    quality_trim(args.fastq_directory)
    make_index()
    map_reads(args.fastq_directory, mapping)
    identify_origins(mapping)
