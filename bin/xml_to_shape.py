#!/usr/bin/env python3
"""
xml_to_shape

Description:
Converts reactivities in XML format (by RNAframework rf-norm) to shape format (https://rna.urmc.rochester.edu/Text/File_Formats.html#SHAPE)

Usage:
python xml_to_shape.py <reactivities.xml> <reactivities.shape>

Authors:
Alfredo Smart & Max Walk
"""


import xml.etree.ElementTree as ET
import sys

input_xml_file = sys.argv[1]
output_shape_file = sys.argv[2]

# parse XML
tree = ET.parse(input_xml_file)
root = tree.getroot()

transcript = root.find('.//transcript')
if transcript is None:
    raise ValueError('Could not find <transcript> element in the XML file.')

sequence = transcript.findtext('sequence', default='').strip().replace('\n', '').replace(' ', '')

# Extract reactivities
reactivities_str = transcript.findtext('reactivity', default='').strip().replace('\n', '').replace(' ', '').replace('\t','')
reactivities = reactivities_str.split(',')

# write cleaned .shape file
with open(output_shape_file, 'w') as f:
    for i, val in enumerate(reactivities, start=1):
        try:
            # Handle missing or invalid values
            if val.lower() == 'nan' or val.strip() == '':
                value = -999
            else:
                num = float(val)
                # Handle negative values
                if num < 0: 
                    value = -999
                else:
                    value = num
        except ValueError:
            print(f'Value error of {val} at position {i} in XML file, setting to -999.')
            value = -999

        f.write(f'{i}\t{value}\n')

print(f'Clean .shape file written: {output_shape_file}')
