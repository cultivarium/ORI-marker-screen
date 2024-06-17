import argparse
import glob
import os
import subprocess
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from unidecode import unidecode


def make_plots(portal_ingest, stats, output_directory):
    for_heatmap = portal_ingest
    for_heatmap.loc[
        for_heatmap["Fold enrichment"] < for_heatmap["Cutoff"], "Fold enrichment"
    ] = None
    for_heatmap["ORI"] = for_heatmap["Plasmid (pGL2)"].str.split("_").str[-1]
    for_heatmap["Strain2"] = for_heatmap["Strain"] + "-" + for_heatmap["Sample_Name"]

    h = for_heatmap.pivot(
        index="Strain2", columns="ORI", values="Fold enrichment"
    ).sort_index()

    if not for_heatmap["Fold enrichment"].isnull().all():

        fig_size = (10, len(h.index) * 0.2)
        fig, ax = plt.subplots(figsize=fig_size)

        heatmap = sns.heatmap(
            h,
            cmap="viridis",
            linewidth=0.1,
            linecolor="grey",
            norm=matplotlib.colors.LogNorm(),
            cbar_kws={"format": "%.0f"},
        )
    else:
        heatmap = sns.heatmap(
            h.fillna(0),
            cmap="viridis",
            linewidth=0.1,
            linecolor="grey",
            cbar_kws={"format": "%.0f"},
        )

    heatmap.set_yticklabels(heatmap.get_yticklabels(), ha="right")
    plt.yticks(range(len(h.index)), h.index)

    plt.savefig(
        os.path.join(output_directory, "heatmap_result.pdf"), bbox_inches="tight"
    )

    # Sort the dataframe by "Matched reads" column
    stats["Matched_P"] = round(stats["Matched"] / stats["Total_Reads"] * 100, 1)

    s_sorted = stats.set_index("Sample").sort_values(by="Matched")

    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 12))

    ax1.axvline(1000, color="red")

    # Plot the "Matched reads" bar chart
    s_sorted["Matched"].plot.barh(ax=ax1)
    ax1.set_xlabel("Number of reads")
    ax1.set_ylabel("")
    ax1.set_xscale("log")  # Set x-axis to log scale

    ax1.set_title("Reads hitting ORI barcodes")

    s_sorted = stats.set_index("Sample").sort_values(by="Matched_P")

    # Plot the "Matched/Total_Reads" bar chart
    s_sorted["Matched_P"].plot.barh(ax=ax2)
    ax2.set_xlabel("Percentage of reads matched (%)")
    ax2.set_ylabel("")
    ax2.set_title("Percentage of total reads matching ORI barcodes")
    ax2.axvline(20, color="red")

    # Adjust the layout and display the plot
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_directory, "quality_results.pdf"), bbox_inches="tight"
    )


def run(fastq_directory, mapping_file, output_directory, unmerged_reads):
    all_lib_info = pd.read_csv("library_info.csv")

    samples = pd.read_csv(mapping_file)
    libraries = list(samples.Library.unique())
    barcodes = set()
    lib_to_cutoff = {}  # The cutoff value for each library

    for library in libraries:
        print("Making Database for {}".format(library))

        f = open("{}.fasta".format(library), "w+")
        lib_info = all_lib_info[all_lib_info.Library == library]
        lib_to_cutoff[library] = lib_info["Negative control cutoff"].to_list()[
            0
        ]  # Should all be the same, take first

        for index, row in lib_info.iterrows():
            barcode = row["ORI"]
            f.write(">" + barcode + "\n")
            f.write(row["ORI-barcode sequence"] + "\n")
            barcodes.add(barcode)
        f.close()

    results = defaultdict(list)  # Form the results table
    portal_ingest = []  # Form the final table for portal ingest
    stats = []  # Form the stats table

    for index, row in samples.iterrows():
        filename = row["FileName"]
        print(filename)
        fn1 = glob.glob(fastq_directory + "/" + filename + "*_R1_*.fastq.gz")

        if len(fn1) == 0:  # file not found
            print("File not found: " + row["Sample"])
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

        ## If file is actually found
        fn1 = fn1[0]
        fn2 = fn1.replace("_R1_", "_R2_")
        base = "/".join(fn1.split("/")[:-1]) + "/" + fn1.split("/")[-1].split("_")[0]

        ## Get library info
        lib_info = all_lib_info[all_lib_info.Library == row["Library"]]
        barcode_to_pgl0 = pd.Series(lib_info.pGL0.values, index=lib_info.ORI).to_dict()
        barcode_to_pgl2 = pd.Series(lib_info.pGL2.values, index=lib_info.ORI).to_dict()

        # Run BBduk
        command = (
            BBTOOLS
            + "bbduk.sh in1={f1} in2={f2} out1={f3} out2={f4} ref=./adapters.fa ktrim=r k=21 qtrim=r trimq=15 maq=15 minlen=30 entropy=0.3 threads=12 forcetrimleft=8 forcetrimright2=8"
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

        if unmerged_reads:
            command = (
                BBTOOLS
                + "bbmerge.sh in1={f1} in2={f2} out={merged}.fasta maxloose=t outu={merged}_u1.fasta outu2={merged}_u2.fasta"
            )
        else:
            command = (
                BBTOOLS + "bbmerge.sh in1={f1} in2={f2} out={merged}.fasta maxloose=t"
            )
        command = command.format(f1=fn1, f2=fn2, merged=base)

        output = subprocess.check_output(
            command, shell=True, stderr=subprocess.STDOUT
        ).decode()
        for line in output.split("\n"):
            if line.startswith("Joined"):
                merged = int(line.split()[1])

        # Cat R1 and merged files
        if unmerged_reads:
            with open(base + ".fasta", "r") as file1, open(
                base + "_u1.fasta", "r"
            ) as file2:
                file1_contents = file1.read()
                file2_contents = file2.read()
            with open(base + "_combined.fasta", "w") as combined_file:
                combined_file.write(file1_contents + "\n" + file2_contents)
            command = "vsearch --usearch_global {base}_combined.fasta --id 0.95 --db {library}.fasta --blast6out {base}.blast"
        else:
            command = "vsearch --usearch_global {base}.fasta --id 0.95 --db {library}.fasta --blast6out {base}.blast"

        ## Run VSEARCH
        command = command.format(
            base=base, library=row["Library"]
        )  # Map to correct library for this sample
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
            if barcode != "Dummy" and barcode in barcode_to_pgl0:
                dat = {
                    "Sample_Name": row["Sample"],
                    "Sample_ID": row["FileName"],
                    "Strain": row["Strain"],
                    "ORI (pGL0)": barcode_to_pgl0[barcode],
                    "Plasmid (pGL2)": barcode_to_pgl2[barcode],
                    "Fold enrichment": round(
                        barcode_hits[barcode] / (1 + barcode_hits["Dummy"]), 2
                    ),
                    "Cutoff": lib_to_cutoff[row["Library"]],
                }
                portal_ingest.append(dat)

        stats_sample = {
            "Sample": row["Sample"],
            "Sample_ID": row["FileName"],
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
    portal_ingest = pd.DataFrame(portal_ingest)

    stats.to_csv(
        os.path.join(output_directory, "barcode_stats.tsv"), sep="\t", index=False
    )
    results.to_csv(
        os.path.join(output_directory, "barcode_results.tsv"), sep="\t", index=False
    )
    portal_ingest.to_csv(
        os.path.join(output_directory, "portal_ingest.tsv"), sep="\t", index=False
    )

    ## Now, let's make some plots
    make_plots(portal_ingest, stats, output_directory)


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
        help="Mapping file of comma separated columns FileName,Sample,Strain,Library.",
    )

    parser.add_argument(
        "-b",
        "--bbmap_folder",
        action="store",
        default="~/bbmap/",
        required=False,
        help="Directory containing BBTools on your system",
    )
    parser.add_argument(
        "-o",
        "--output_folder",
        action="store",
        default=".",
        required=False,
        help="Directory to store output files",
    )

    parser.add_argument(
        "--unmerged_reads",
        action="store_true",
        default=False,
        required=False,
        help="Also process unmerged R1 (useful with 2x75 bp or lower quality reads)",
    )

    args = parser.parse_args()

    BBTOOLS = args.bbmap_folder.rstrip("/") + "/"
    mapping = pd.read_csv(args.mapping_file)

    run(
        args.fastq_directory.rstrip("/"),
        args.mapping_file,
        args.output_folder.rstrip("/"),
        args.unmerged_reads,
    )
