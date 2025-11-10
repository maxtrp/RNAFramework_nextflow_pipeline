#!/usr/bin/env python3

"""
Adapted from nf-core/rnaseq v3.12.0:
https://github.com/nf-core/rnaseq/blob/3.12.0/bin/check_samplesheet.py
"""

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the input samplesheet follows the following structure:

    sample,treated_fastq_1,treated_fastq_2,control_fastq_1,control_fastq_2
    sample_name,/path/to/fastqs/treated_sample_R1.fastq.gz,/path/to/fastqs/treated_sample_R2.fastq.gz,/path/to/fastqs/control_sample_R1.fastq.gz,/path/to/fastqs/control_sample_R2.fastq.gz

    Treatment and control channels are separated in the output csv:
    sample,treatment,fastq_1,fastq_2
    sample_name,treated,/path/to/fastqs/treated_sample_R1.fastq.gz,/path/to/fastqs/treated_sample_R2.fastq.gz
    sample_name,control,/path/to/fastqs/control_sample_R1.fastq.gz,/path/to/fastqs/control_sample_R2.fastq.gz
    """

    sample_mapping_list = []
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        HEADER = ["sample", "treated_fastq_1", "treated_fastq_2", "control_fastq_1", "control_fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}")
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            if line.strip():
                lspl = [x.strip().strip('"') for x in line.strip().split(",")]

                ## Check valid number of columns per row
                if len(lspl) < len(HEADER):
                    print_error(
                        f"Invalid number of columns (minimum = {len(HEADER)})!",
                        "Line",
                        line,
                    )

                ## Extract sample information
                sample, treated_fastq_1, treated_fastq_2, control_fastq_1, control_fastq_2 = lspl[: len(HEADER)]
                if sample.find(" ") != -1:
                    print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                    sample = sample.replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)

                ## Check FastQ file extensions and contents
                for fastq in [treated_fastq_1, treated_fastq_2, control_fastq_1, control_fastq_2]:
                    if fastq:
                        if fastq.find(" ") != -1:
                            print_error("FastQ file contains spaces!", "Line", line)
                        if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )

                ## Ensure both treated_fastq_1 and treated_fastq_2 are provided for treated samples
                if not treated_fastq_1 or not treated_fastq_2:
                    print_error("Both treated_fastq_1 and treated_fastq_2 must be provided for treated samples!", "Line", line)

                ## Ensure both control_fastq_1 and control_fastq_2 are provided for control samples
                if not control_fastq_1 or not control_fastq_2:
                    print_error("Both control_fastq_1 and control_fastq_2 must be provided for control samples!", "Line", line)

                ## Add sample entries to the list
                sample_mapping_list.append([sample, "treated", treated_fastq_1, treated_fastq_2])
                sample_mapping_list.append([sample, "control", control_fastq_1, control_fastq_2])

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_list) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample",  "treatment", "fastq_1", "fastq_2"]) + "\n")
            for entry in sample_mapping_list:
                fout.write(",".join(entry) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
