#!/bin/env bash

export MCDBIN="${HOME}/programs.aa/mcd-molcas/mcd-molcas"

# Path to the options file
OPTIONS_TEMPLATE="options.template"
OPTIONS_FILE="options.dat"

# Define the range and step size for the temperature (in Kelvin). Needed for C-term spectra.
START_TEMP=0
END_TEMP=50
STEP=5

# Ensure output directory exists
OUTPUT_DIR="output"
rm -fr "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

for TEMP in $(seq $START_TEMP $STEP $END_TEMP); do
  # Format the temperature with leading zeros
  FORMATTED_TEMP=$(printf "%02d" $TEMP)

  echo "Running calculation at ${FORMATTED_TEMP}K..."

  # Update the options.dat file with the current temperature
  sed "s/temp = TEMP/temp = $TEMP/g" "$OPTIONS_TEMPLATE" > "$OPTIONS_FILE"

  # Run the calculation
  ${MCDBIN}/mcd-c-molcas
  ${MCDBIN}/plot-mcdspectrum < mcd-c-spectrum-0

  # Rename and move output files with formatted temperature
  mv graph.dat "${OUTPUT_DIR}/graph-c-${FORMATTED_TEMP}K.dat"
  mv impulses.dat "${OUTPUT_DIR}/impulses-c-${FORMATTED_TEMP}K.dat"

done

echo "All calculations are completed."


