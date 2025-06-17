#!/bin/bash

# Loop over Fortran files and apply substitution
find . -type f \( -iname "*.f90" -o -iname "*.f95" -o -iname "*.f03" -o -iname "*.f08" -o -iname "*.f" -o -iname "*.for" \) | while read -r file; do
    echo "Processing $file"
    sed -E -i.bak 's/\bWCALL\s*\(\s*[^,]+,\s*([^,]+),\s*([^)]+)\s*\)/call \1\2/g' "$file"
done
