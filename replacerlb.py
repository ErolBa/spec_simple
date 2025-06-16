import re
import os

def split_args(text):
    """
    Split arguments inside parentheses, respecting nested parentheses.
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
        stripped = line.strip()
        match = re.match(r'^(\s*)RlBCAST\s*\((.*)\)\s*$', stripped)
        if match:
            indent = match.group(1)
            args_text = match.group(2)

            try:
                args = split_args(args_text)
                if len(args) == 3:
                    a1, a2, a3 = args
                    new_line = f"{indent}call MPI_BCAST({a1}, {a2}, MPI_DOUBLE_PRECISION, {a3}, MPI_COMM_SPEC, ierr)\n"
                    new_lines.append(new_line)
                    changed = True
                    continue
            except Exception as e:
                print(f"⚠️ Skipping malformed RlBCAST in {path}: {line.strip()}")

        new_lines.append(line)

    if changed:
        with open(path, 'w') as f:
            f.writelines(new_lines)
        print(f"✅ Modified: {path}")

def find_fortran_files(directory):
    extensions = ('.f90', '.f95', '.f03', '.f08', '.f', '.for')
    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith(extensions):
                yield os.path.join(root, file)

if __name__ == '__main__':
    import sys
    directory = sys.argv[1] if len(sys.argv) > 1 else '.'
    for filepath in find_fortran_files(directory):
        process_file(filepath)

