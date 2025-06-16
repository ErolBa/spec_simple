#/usr/bin/env python

import os
import re

def find_fortran_files(directory):
    """Return a list of Fortran source files in the directory tree."""
    return [os.path.join(root, file)
            for root, _, files in os.walk(directory)
            for file in files if file.lower().endswith(('.f90', '.f95', '.f03', '.f08', '.f', '.for'))]

def extract_subroutine_names(content):
    """Extract subroutine names using regex."""
    pattern = re.compile(r'^\s*subroutine\s+(\w+)', re.IGNORECASE | re.MULTILINE)
    return pattern.findall(content)

def extract_all_words(content):
    """Return all words from code (excluding comments)."""
    # Remove Fortran-style comments
    code_lines = [line.split('!')[0] for line in content.splitlines()]
    code = '\n'.join(code_lines)
    return re.findall(r'\b\w+\b', code, re.IGNORECASE)

def main(directory='.'):
    subroutines = {}     # subroutine name -> file
    all_words = []       # list of all words used in all files (excluding comments)

    # Pass 1: Collect definitions and all words from all files
    for file in find_fortran_files(directory):
        with open(file, 'r', errors='ignore') as f:
            content = f.read()

        for name in extract_subroutine_names(content):
            subroutines[name.lower()] = file

        words = extract_all_words(content)
        all_words.extend(word.lower() for word in words)

    # Pass 2: Count how often each subroutine name appears
    unused = [name for name in subroutines if all_words.count(name) < 3]

    # Report
    if unused:
        print("Unused subroutines (only mentioned once):")
        for name in unused:
            print(f" - {name} (defined in {subroutines[name]})")
    else:
        print("✅ No unused subroutines found.")

if __name__ == '__main__':
    import sys
    target = sys.argv[1] if len(sys.argv) > 1 else '.'
    main(target)

