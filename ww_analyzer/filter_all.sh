#!/bin/bash

# Checks for an index file
if [ "$#" -ne 1 ]; then
  echo "Error: No index file specified"
  echo "Usage: $0 <path/to/index_file.txt>"
  exit 1
fi

# Creates the output file name
IFS='/' read -ra INPUT <<<"$1"
IFS='_' read -ra NAME <<<${INPUT[-1]}
OUTPUT_NAME="${NAME[0]}_${NAME[1]}_${NAME[2]}_${NAME[6]}.root"

# Checks if directory exists. If not, creates it
DIR="rootfiles_${NAME[0]}_${NAME[1]}_${NAME[2]}_${NAME[6]}"
if [ ! -d $DIR ]; then
  echo "Creating output directory as $DIR"
  mkdir $DIR
  cd $DIR
else
  echo "Output directory already exists"
  cd $DIR
fi

i=0
while IFS= read -r line; do
  OUTPUT_NAME_PARTIAL="${NAME[0]}_${NAME[1]}_${NAME[2]}_${NAME[6]}_p$i.root"
  macro_command="../filter.C( \"${line}\", \"${OUTPUT_NAME_PARTIAL}\")"
  root -l -q "${macro_command}"
  ((i = i + 1))
done <"../${INPUT[-1]}"

hadd -f "../$OUTPUT_NAME" *.root
