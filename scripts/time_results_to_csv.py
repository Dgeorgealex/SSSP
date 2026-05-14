#!/usr/bin/env python3
"""
Parse time_*.txt files from data/graphs and generate a CSV with mean times for original, v2, and v3 versions.

Usage:
    python time_results_to_csv.py [--algorithm ALGO] [--output OUTPUT.csv]
    
    --algorithm: Specify which algorithm to extract (e.g., 'GOR', 'PADSCALING'). 
                 If not specified, uses the first algorithm found in each file.
    --output: Output CSV file path. Default: time_results.csv
"""

import os
import re
import sys
import argparse
from collections import defaultdict
from pathlib import Path

def extract_mean_value(file_path, algorithm=None):
    """
    Extract mean value from a time file.
    
    Args:
        file_path: Path to the time file
        algorithm: Optional algorithm name to search for (e.g., 'GOR'). 
                   If None, returns the first mean found.
    
    Returns:
        Float mean value or -1 if parsing fails or algorithm timed out
    """
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check if the algorithm timed out
        # if algorithm:
        #     timeout_pattern = rf'TIMEOUT:.*?{algorithm}'
        #     if re.search(timeout_pattern, content, re.IGNORECASE | re.DOTALL):
        #         return -1
        
        # Pattern to match: ALGORITHM_NAME - ... mean = XXXX ms
        if algorithm:
            # Look for specific algorithm
            pattern = rf'\b{algorithm}\b.*?mean\s*=\s*(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*ms'
            match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        else:
            # Find any mean value
            pattern = r'mean\s*=\s*(\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*ms'
            match = re.search(pattern, content)
        
        if match:
            return float(match.group(1))
        else:
            return -1
    except Exception as e:
        print(f"Error parsing {file_path}: {e}", file=sys.stderr)
        return -1

def get_base_name(filename):
    """
    Extract base name from time_*.txt filename.
    
    Examples:
        time_big_aug_gor_2e6_6.txt -> big_aug_gor_2e6_6
        time_big_aug_gor_2e6_6_v2.txt -> big_aug_gor_2e6_6
        time_big_aug_gor_2e6_6_v3.txt -> big_aug_gor_2e6_6
    """
    # Remove 'time_' prefix and '.txt' suffix
    name = filename[5:-4] if filename.startswith('time_') else filename[:-4]
    
    # Remove version suffix
    name = re.sub(r'_(v2|v3)', '', name)
    
    return name

def get_version(filename):
    """
    Extract version from filename.
    Returns: 'original', 'v2', or 'v3'
    """
    if '_v2' in filename:
        return 'v2'
    elif '_v3' in filename:
        return 'v3'
    else:
        return 'original'

def main():
    parser = argparse.ArgumentParser(
        description='Parse time_*.txt files and generate CSV with results'
    )
    parser.add_argument('--algorithm', default='GOR',help='Algorithm to extract (e.g., GOR, PADSCALING)')
    parser.add_argument('--output', default='time_results.csv', help='Output CSV file')
    parser.add_argument('--dir', default='data/graphs/', help='Directory containing time_*.txt files')
    
    args = parser.parse_args()
    
    # Find all time_*.txt files
    data_dir = Path(args.dir)
    if not data_dir.exists():
        print(f"Error: Directory {args.dir} does not exist", file=sys.stderr)
        sys.exit(1)
    
    time_files = sorted(data_dir.glob('time_*.txt'))
    
    if not time_files:
        print(f"No time_*.txt files found in {args.dir}", file=sys.stderr)
        sys.exit(1)
    
    # Group files by base name
    grouped_files = defaultdict(dict)
    
    for file_path in time_files:
        filename = file_path.name
        base_name = get_base_name(filename)
        version = get_version(filename)
        grouped_files[base_name][version] = file_path
    
    # Collect results for each base_name
    results = {}
    for base_name in grouped_files:
        versions = grouped_files[base_name]
        means = {}
        for version_type in ['original', 'v2', 'v3']:
            if version_type in versions:
                mean = extract_mean_value(versions[version_type], args.algorithm)
                means[version_type] = mean
            else:
                means[version_type] = None
        valid_means = [m for m in means.values() if m is not None and m >= 0]
        average = (sum(valid_means) / len(valid_means)) if valid_means else None
        results[base_name] = {'means': means, 'average': average}

    # Desired output ordering: first entries with no suffix (nothing), then '_t', then '_0', then '_-1'
    group_order = ['nothing', 't', '0', '-1']
    edge_order = ['5e5', '2e6', '5e6', '1e7', '2e7']
    edge_index = {tok: i for i, tok in enumerate(edge_order)}

    buckets = {g: [] for g in group_order}
    others = []

    for name in results:
        if '_-1' in name:
            buckets['-1'].append(name)
        elif '_0' in name:
            buckets['0'].append(name)
        elif '_t' in name:
            buckets['t'].append(name)
        else:
            buckets['nothing'].append(name)

    def sort_key(n):
        m = re.search(r'(5e5|2e6|5e6|1e7|2e7)', n)
        if m:
            return (edge_index[m.group(1)], n)
        return (len(edge_order), n)

    # Flatten rows in the requested group sequence (no explicit group headers)
    csv_rows = []
    csv_rows.append(['file_name', 'original', 'v2', 'v3', 'average'])

    for grp in group_order:
        for base_name in sorted(buckets[grp], key=sort_key):
            means = results[base_name]['means']
            avg = results[base_name]['average']
            row = [base_name]
            for version_type in ['original', 'v2', 'v3']:
                val = means[version_type]
                row.append(str(int(val)) if val is not None else '')
            row.append(str(int(avg)) if avg is not None else '')
            csv_rows.append(row)

    # Append any names not matching the groups (should be none)
    for base_name in sorted(others, key=sort_key):
        means = results[base_name]['means']
        avg = results[base_name]['average']
        row = [base_name]
        for version_type in ['original', 'v2', 'v3']:
            val = means[version_type]
            row.append(str(int(val)) if val is not None else '')
        row.append(str(int(avg)) if avg is not None else '')
        csv_rows.append(row)
    
    # Write CSV
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        for row in csv_rows:
            f.write(','.join(row) + '\n')
    
    print(f"CSV written to {output_path}")
    print(f"Processed {len(grouped_files)} unique graph configurations")
    if args.algorithm:
        print(f"Algorithm: {args.algorithm}")

if __name__ == '__main__':
    main()
