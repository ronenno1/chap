import json
import csv
import sys
import os

def get_leaves(item, key=None):
    if isinstance(item, dict):
        leaves = []
        for i in item.keys():
            if key==None:
                leaves.extend(get_leaves(item[i], i))
            else:
                leaves.extend(get_leaves(item[i], key+'_'+i))
        return leaves
    else:
        return [(key, item)]

with open(sys.argv[1]+'.txt') as f_input, open(sys.argv[1]+'.json', 'w') as f_output:
    text = f_input.read()
    text =  text.replace('}}}','}}},')
    text = '[' + text.replace('{"category":"heartbeat","statuscode":200}','') + '{}]'
    f_output.write(text)
    
with open(sys.argv[1]+'.json') as f_input, open(sys.argv[1]+'.csv', 'w') as f_output:
    csv_output = csv.writer(f_output)
    write_header = True
    for entry in json.load(f_input):
        leaf_entries = sorted(get_leaves(entry))
        if write_header:
            csv_output.writerow([k for k, v in leaf_entries])
            write_header = False
        csv_output.writerow([v for k, v in leaf_entries])
