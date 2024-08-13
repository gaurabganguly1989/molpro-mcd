#!/usr/bin/env bash

# export Project="nip"
# file=$Project.cas_16e_11o.out

########################
# DATA post processing #
########################

if [[ -z $1 ]];
then
    echo "no input passed."
else
    echo "input passed = $1"
fi

file=$1

# create energy files
grep "\!MCSCF.*1 Energy" $file | awk '{print $(NF)}' > MCSCF_Energy_Ag
grep "\!MCSCF.*2 Energy" $file | awk '{print $(NF)}' > MCSCF_Energy_B3u
grep "\!MCSCF.*3 Energy" $file | awk '{print $(NF)}' > MCSCF_Energy_B2u

MCSCF_Ag=`awk '{print $1}' MCSCF_Energy_Ag | tail -1`
awk -v var="$MCSCF_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' MCSCF_Energy_B3u > MCSCF_Energy_B3u_eV
awk -v var="$MCSCF_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' MCSCF_Energy_B2u > MCSCF_Energy_B2u_eV

lines=`wc -l MCSCF_Energy_B2u_eV | awk '{print $1}'`

grep "\!MCSCF trans.*1|DMX|.*.2>" $file | awk '{print $4}' | head -$lines > TRDX_B3u
grep "\!MCSCF trans.*1|DMX|.*.3>" $file | awk '{print $4}' | head -$lines > TRDX_B2u
grep "\!MCSCF trans.*1|DMY|.*.2>" $file | awk '{print $4}' | head -$lines > TRDY_B3u
grep "\!MCSCF trans.*1|DMY|.*.3>" $file | awk '{print $4}' | head -$lines > TRDY_B2u

paste MCSCF_Energy_B3u_eV MCSCF_Energy_B2u_eV TRDX_B3u TRDY_B2u > temp

awk '{print NR "," $0}' temp > temp1

echo "# st  Ag->B3u  Ag->B2u  <Ag|MuX|B3u>  <Ag|MuY|B2u> " > 2A1g-to-2Eu

# MCSCF dipole thr = 0.4 au.
awk '{if (sqrt($4*$4)>=0.4) printf " %4s  %8.3f %8.3f %12.5f %12.5f %8s\n", $1, $2, $3, $4, $5, "<-bright";
      else printf " %4s  %8.3f %8.3f %12.5f %12.5f\n", $1, $2, $3, $4, $5}' temp1 >> 2A1g-to-2Eu

# awk '{printf " %4s  %8.2f %8.2f %12.5f %12.5f\n", $1, $2, $3, $4, $5}' temp1 >> 2A1g-to-2Eu
# grep "\!NEVPT2 STATE" $file | grep "1 Energy" | awk '{print $(NF)}' > NEVPT2_Ag
# grep "\!NEVPT2 STATE" $file | grep "2 Energy" | awk '{print $(NF)}' > NEVPT2_B3u
# grep "\!NEVPT2 STATE" $file | grep "3 Energy" | awk '{print $(NF)}' > NEVPT2_B2u

grep '\!NEVPT2 STATE.*1 Energy' $file | awk '{print $(NF)}' > pc_NEVPT2_Ag
grep '\!NEVPT2 STATE.*2 Energy' $file | awk '{print $(NF)}' > pc_NEVPT2_B3u
grep '\!NEVPT2 STATE.*3 Energy' $file | awk '{print $(NF)}' > pc_NEVPT2_B2u

grep -A1 '\!NEVPT2 STATE.*1 Energy' $file | grep -v '\!NEVPT2 STATE.*1 Energy' | grep "Strongly contracted energy" | awk '{print $(NF)}' > sc_NEVPT2_Ag
grep -A1 '\!NEVPT2 STATE.*2 Energy' $file | grep -v '\!NEVPT2 STATE.*2 Energy' | grep "Strongly contracted energy" | awk '{print $(NF)}' > sc_NEVPT2_B3u
grep -A1 '\!NEVPT2 STATE.*3 Energy' $file | grep -v '\!NEVPT2 STATE.*3 Energy' | grep "Strongly contracted energy" | awk '{print $(NF)}' > sc_NEVPT2_B2u

pc_NEV_Ag=`awk '{print $1}' pc_NEVPT2_Ag | tail -1`
sc_NEV_Ag=`awk '{print $1}' sc_NEVPT2_Ag | tail -1`

awk -v var="$pc_NEV_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' pc_NEVPT2_B3u > pc_NEVPT2_B3u_eV
awk -v var="$sc_NEV_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' sc_NEVPT2_B3u > sc_NEVPT2_B3u_eV
awk -v var="$pc_NEV_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' pc_NEVPT2_B2u > pc_NEVPT2_B2u_eV
awk -v var="$sc_NEV_Ag" '{printf " %8.3f\n", ($1-var)*27.2114}' sc_NEVPT2_B2u > sc_NEVPT2_B2u_eV

# echo "pc-NEVPT2 sc-NEVPT2" > NEVPT2_B3u_eV
# echo "pc-NEVPT2 sc-NEVPT2" > NEVPT2_B2u_eV

paste  pc_NEVPT2_B3u_eV sc_NEVPT2_B3u_eV >  NEVPT2_B3u_eV
paste  pc_NEVPT2_B2u_eV sc_NEVPT2_B2u_eV >  NEVPT2_B2u_eV

cp NEVPT2_B3u_eV NEVPT2_B2u_eV

echo "#                   Ag->B3u                              Ag->B2u                                                " >  Result
echo "# st    E(CASSCF) E(pc-NEVPT) E(sc-NEVPT)    E(CASSCF) E(pc-NEVPT) E(sc-NEVPT)    <Ag|MuX|B3u>    <Ag|MuY|B2u>  " >> Result

grep "<-bright" 2A1g-to-2Eu > bright
# grep " 4," 2A1g-to-2Eu > bright
# grep "12," 2A1g-to-2Eu >> bright
paste bright NEVPT2_B3u_eV NEVPT2_B2u_eV > temp_bright

awk '{printf " %4s  %6.3f %10.3f %10.3f %14.3f %10.3f %10.3f %18.6f %15.6f\n", $1, $2, $7, $8, $3, $9, $10, $4, $5}' temp_bright >> Result

echo " " >> Result
echo " " >> Result

echo "# st    E(CASSCF) E(pc-NEVPT) E(sc-NEVPT)    E(CASSCF) E(pc-NEVPT) E(sc-NEVPT)    <Ag|MuX|B3u>    <Ag|MuY|B2u>  pc-NEVPT     sc-NEVPT " >> Result

awk '{printf " %4s  %8'\''d %10'\''d %10'\''d %14'\''d %10'\''d %10'\''d %16.6f %15.6f %12.6E %12.6E\n", $1, $2*8065.54, $7*8065.54, $8*8065.54, $3*8065.54, $9*8065.54, $10*8065.54, $4, $5, 2*(2/3)*$7/27.2114*$4*$4, 2*(2/3)*$8/27.2114*$4*$4}' temp_bright >> Result

rm -f Energy_* TRDX_* TRDY_* temp* pc_NEVPT2* sc_NEVPT2* NEVPT2_* MCSCF_* MRCI_* bright

# THE END #





