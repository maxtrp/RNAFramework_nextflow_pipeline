#!/usr/bin/env python3

"""
raw_counts.py

Description:
Outputs a csv file containing mutation and coverage counts for the treated and control samples.

Usage:
python raw_counts.py <treated_counts.txt> <control_counts.txt>

Author:
Max Walk
"""

import sys
import pandas as pd

def load_counts(counts_txt_file):

    num_fields = 5 # number of lines associated with a given transcript (nb last line is empty)

    with open(counts_txt_file, 'r') as infile:

        all_lines = [line.strip() for line in infile]

    num_transcripts = len(all_lines)//num_fields

    counts_data_dict = {} # for storing per-position data
    transcript_seq_dict = {} # for storing transcript sequence
    
    for entry_num in range(num_transcripts):

        transcript_id = all_lines[(num_fields*entry_num)]
        sequence_DNA = all_lines[(num_fields*entry_num) + 1]
        mutations = [int(val) for val in all_lines[(num_fields*entry_num) + 2].split(',')]
        coverage = [int(val) for val in all_lines[(num_fields*entry_num) + 3].split(',')]

        assert(len(sequence_DNA) == len(mutations)), f"Invalid entry in {counts_txt_file} at line {num_fields*entry_num}"

        counts_data_dict[transcript_id] = {
            'mutations': mutations,
            'coverage': coverage
        }

        # convert transcript sequence to RNA
        sequence_RNA = sequence_DNA.replace('T', 'U')
        transcript_seq_dict[transcript_id] = list(sequence_RNA)
    
    return counts_data_dict, transcript_seq_dict

def main():

    treated_counts_txt = sys.argv[1]
    control_counts_txt = sys.argv[2]

    treated_counts_dict, treated_transcripts_dict = load_counts(treated_counts_txt)
    control_counts_dict, control_transcripts_dict = load_counts(control_counts_txt)

    assert treated_transcripts_dict.keys() == control_transcripts_dict.keys(), "Transcripts in treated/control counts files dont match!"

    transcript_ids = [transcript for transcript in treated_counts_dict.keys()]

    for transcript in transcript_ids:

        transcript_sequence = treated_transcripts_dict[transcript]
        
        treated_data = treated_counts_dict[transcript]
        control_data = control_counts_dict[transcript]

        treated_df = pd.DataFrame.from_dict(treated_data)
        control_df = pd.DataFrame.from_dict(control_data)

        treated_df.columns = ['treated_' + col for col in treated_df.columns]
        control_df.columns = ['control_' + col for col in control_df.columns]

        combined_df = pd.concat([treated_df, control_df], axis=1)

        # add nucleotide/position columns
        combined_df.insert(0, 'nucleotide', transcript_sequence)
        combined_df.insert(0, 'position', range(1, len(combined_df) + 1))

        combined_df.to_csv(f'{transcript}_raw_counts.csv', index=False)

if __name__ == "__main__":
    main()