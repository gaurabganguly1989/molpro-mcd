#!/bin/env bash

export MCDBIN="${HOME}/programs.aa/mcd-molcas/mcd-molcas"

${MCDBIN}/mcd-a-molcas 
${MCDBIN}/mcd-b-molcas 
${MCDBIN}/mcd-c-molcas 

for x in "a" "b" "c"
do
  ${MCDBIN}/plot-mcdspectrum < mcd-"$x"-spectrum-0
  mv graph.dat graph-$x.dat
  mv impulses.dat impulses-$x.dat
done
