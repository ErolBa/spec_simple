#/usr/bin/env python

import os
import re

def find_fortran_files(directory):
    return [os.path.join(root, f)
            for root, _, files in os.walk(directory)
            for f in files if f.lower().endswith(('.f90', '.f95', '.f03', '.f08', '.f', '.for'))]

def extract_subroutine_definitions(content):
    """Return a list of subroutine names defined in the content."""
    pattern = re.compile(r'^\s*subroutine\s+(\w+)', re.IGNORECASE | re.MULTILINE)
    return pattern.findall(content)

def extract_subroutine_calls(content):
    """Return a list of subroutine names that are called via 'call xxx' or 'WCALL(..., xxx)'."""
    calls = set()

    # Match: call xxx
    call_pattern = re.compile(r'\bcall\s+(\w+)', re.IGNORECASE)
    calls.update(call_pattern.findall(content))

    # Match: WCALL(anything, xxx) — extract second argument (ignoring whitespace)
    wcall_pattern = re.compile(r'\bW?CALL\s*\([^,]+,\s*(\w+)', re.IGNORECASE)
    calls.update(wcall_pattern.findall(content))

    return list(calls)

def remove_comments(content):
    """Remove Fortran-style comments (everything after a !)"""
    lines = content.splitlines()
    return '\n'.join(line.split('!')[0] for line in lines)

def main(directory='.'):
    subroutine_defs = {}   # name -> file
    called_subroutines = set()

    for filepath in find_fortran_files(directory):
        with open(filepath, 'r', errors='ignore') as f:
            raw_content = f.read()
            content = remove_comments(raw_content)

        # Definitions
        for name in extract_subroutine_definitions(content):
            subroutine_defs[name.lower()] = filepath

        # Calls
        for call in extract_subroutine_calls(content):
            called_subroutines.add(call.lower())

    # Compare
    unused_subroutines = [name for name in subroutine_defs if name not in called_subroutines]

    # Report
    if unused_subroutines:
        print("Unused subroutines:")
        for name in unused_subroutines:
            print(f" - {name} (defined in {subroutine_defs[name]})")
    else:
        print("✅ No unused subroutines found.")

if __name__ == '__main__':
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else '.'
    main(path)

