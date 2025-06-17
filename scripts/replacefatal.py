import re
import os

def split_args(text):
    """
    Split arguments in FATAL(...) handling nested parentheses and commas.
    Returns a list of 3 elements: [msg, condition, extra]
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

def ensure_string_literal(s):
    """
    Ensure a string is quoted (for Fortran string literals).
    """
    s = s.strip()
    return s if (s.startswith('"') and s.endswith('"')) else f'"{s}"'

def process_file(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    changed = False

    for line in lines:
        # Remove trailing comment and whitespace
        code_line = line.split('!')[0].rstrip()

        match = re.match(r'^\s*\bFATAL\s*\((.*)\)\s*$', code_line)
        if match:
            try:
                args = split_args(match.group(1))
                if len(args) == 3:
                    msg, cond, extra = args

                    msg_q = ensure_string_literal(msg)
                    extra_q = ensure_string_literal(extra)

                    # Build replacement block
                    block = (
                        f"if( {cond} ) then\n"
                        f"     write(6,'({msg_q} :      fatal : myid=\",i3,\" ; {cond} ; {extra};\")') myid\n"
                        f"     call MPI_ABORT( MPI_COMM_SPEC, 1, ierr )\n"
                        f"     stop {msg_q} // ' : ' // '{cond}' // ' : ' // {extra_q} // ' ;'\n"
                        f"   endif\n"
                    )
                    new_lines.append(block)
                    changed = True
                    continue
            except Exception as e:
                print(f"⚠️ Error parsing FATAL in {path}: {line.strip()} — {e}")
        new_lines.append(line)

    if changed:
        with open(path, 'w') as f:
            f.writelines(new_lines)
        print(f"✅ Modified: {path}")

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

