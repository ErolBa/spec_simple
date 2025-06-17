#/usr/bin/env python

import os
import re

def split_args(text):
    """
    Split arguments inside SALLOCATE(...) respecting nested parentheses.
    Returns a list of [xxx, yyy, zzz]
    """
    args = []
    current = ''
    depth = 0
    for c in text:
        if c == ',' and depth == 0:
            args.append(current.strip())
            current = ''
        else:
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            current += c
    if current:
        args.append(current.strip())
    return args

def process_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    changed = False

    for line in lines:
        code_only = line.split('!')[0].rstrip()

        # Match SALLOCATE with outermost parenthesis
        match = re.match(r'^\s*\bSALLOCATE\s*\((.*)\)\s*$', code_only)
        if match:
            args_text = match.group(1)
            try:
                args = split_args(args_text)
                if len(args) == 3:
                    xxx, yyy, zzz = args
                    joined = f"{xxx}{yyy}".replace(' ', '')
                    replacement = (
                        f"if( allocated( {xxx} ) ) deallocate( {xxx} )\n"
                        f"allocate( {joined}, stat=astat )\n"
                        f"{joined} = {zzz}\n"
                    )
                    new_lines.append(replacement)
                    changed = True
                    continue  # skip appending original line
            except Exception as e:
                print(f"Error parsing line in {path}: {line.strip()} — {e}")
        new_lines.append(line)

    if changed:
        with open(path, 'w') as f:
            f.writelines(new_lines)
        print(f"Modified: {path}")

def find_fortran_files(root):
    extensions = ('.f90', '.f95', '.f03', '.f08', '.f', '.for')
    for dirpath, _, filenames in os.walk(root):
        for fname in filenames:
            if fname.lower().endswith(extensions):
                yield os.path.join(dirpath, fname)

if __name__ == '__main__':
    import sys
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    for file_path in find_fortran_files(directory):
        process_file(file_path)

