#!/usr/bin/env python3
"""
fix_draco_json

Description:
Adds 'filenames' field to draco json file to enable downstream processing using rf-json2rc. 
Fixes incompatibility between draco v1.3 and RNA Framework v2.9.5.

Usage:
python fix_draco_json.py <draco.json>

Authors:
Max Walk
"""

import sys
import json

def main():

    input_json = sys.argv[1]

    input_json_filename = input_json.split('/')[-1]
    output_json_filename = f'fixed_{input_json_filename}'

    with open(input_json, 'r') as f: 
        json_dict = json.load(f)

    json_dict['filenames'] = [input_json_filename]

    with open(output_json_filename, 'w') as f:
        json.dump(json_dict, f)

if __name__ == '__main__':
    main()