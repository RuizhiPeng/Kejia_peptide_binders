#!/usr/bin/env python3
"""
Script to replace path prefixes in a text file containing PDB file paths.
"""

import os
import argparse

def replace_path_prefix(input_file, output_file, old_prefix, new_prefix):
    """
    Replace path prefixes in a text file.
    
    Args:
        input_file (str): Path to input file
        output_file (str): Path to output file
        old_prefix (str): Old prefix to replace
        new_prefix (str): New prefix to use
    """
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        updated_lines = []
        for line in lines:
            line = line.strip()
            if line.startswith(old_prefix):
                # Replace the old prefix with the new prefix
                relative_path = line[len(old_prefix):]
                # Remove leading slash if present to avoid double slashes
                if relative_path.startswith('/'):
                    relative_path = relative_path[1:]
                new_line = os.path.join(new_prefix, relative_path)
                updated_lines.append(new_line + '\n')
            else:
                # Keep lines that don't match the old prefix unchanged
                updated_lines.append(line + '\n')
        
        with open(output_file, 'w') as f:
            f.writelines(updated_lines)
        
        print(f"Successfully updated {len(updated_lines)} lines")
        print(f"Output written to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
    except Exception as e:
        print(f"Error: {e}")

def main():
    # Default values based on your example
    default_old_prefix = "/net/scratch/kejiawu/pp/beta"
    default_new_prefix = "/storage/d1/users/ruizhi/logos"
    
    parser = argparse.ArgumentParser(description="Replace path prefixes in PDB file path lists")
    parser.add_argument("input_file", help="Input text file containing PDB paths")
    parser.add_argument("-o", "--output", help="Output file (default: input_file with _updated suffix)")
    parser.add_argument("--old-prefix", default=default_old_prefix, 
                       help=f"Old prefix to replace (default: {default_old_prefix})")
    parser.add_argument("--new-prefix", default=default_new_prefix,
                       help=f"New prefix to use (default: {default_new_prefix})")
    
    args = parser.parse_args()
    
    # Set default output filename if not provided
    if not args.output:
        base_name = os.path.splitext(args.input_file)[0]
        extension = os.path.splitext(args.input_file)[1]
        args.output = f"{base_name}_updated{extension}"
    
    print(f"Input file: {args.input_file}")
    print(f"Output file: {args.output}")
    print(f"Replacing '{args.old_prefix}' with '{args.new_prefix}'")
    print("-" * 50)
    
    replace_path_prefix(args.input_file, args.output, args.old_prefix, args.new_prefix)

if __name__ == "__main__":
    main()
