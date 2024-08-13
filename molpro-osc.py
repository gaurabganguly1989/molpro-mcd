#!/usr/bin/env python

"""
Copyright (c) [2024] 
Dr. Gaurab Ganguly
Email: gaurabganguly1989@gmail.com
University of Vienna

"""

import numpy as np
import pandas as pd
import math
import os

###########################################################################
# PARAMETERS FOR DATA
###########################################################################

# Irrep = number of symmetry considered in calculation
nIrrep = 4

# SFS = Dictionary of spin-free states (SFS) per irrep
SFS = {1:  1,
       2:  1,
       3: 60,
       4: 60}

# total number of spin-free states (SFS)
nSFS = sum(SFS.values())

# MULT = multiplicity of each spin-free state [ 4 = quartet ]
MULT = {i: 4 for i in range(1, nIrrep + 1)}

# SOS = total number of spin-orbit states (SOS)
nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

# degeneracy of the SF and SO ground states
sfdgen = 2
sodgen = sfdgen * MULT[1]

# MOLPRO file.out
fout = "cop.nevpt-so.out"

###########################################################################
#  SF-ENERGIES
###########################################################################
def parse_sf_energy(fout, SFS, nSFS):
    """
    Parses spin-free energies from a MOLPRO log file.

    Parameters:
        fout (str): Path to the MOLPRO log file.
        SFS (dict): Number of spin-free states per symmetry in the calculation.
        nSFS (int): Total number of spin-free states.

    Returns:
        numpy.ndarray: Array containing spin-free energies.
    """
    try:
        with open(fout, 'r') as f:
            lines = f.readlines()
            sf_energy = np.zeros(nSFS, dtype=float)

            i = 0
            iSFS = 0
            while i < len(lines):
                if f"SETTING HLSDIAG({iSFS + 1}" in lines[i]:
                    sf_energy[iSFS] = float(lines[i].split()[-1])
                    iSFS += 1
                i += 1 

        return sf_energy
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred while parsing spin-free energies: {str(e)}")
        return None

###########################################################################
#  SO-ENERGIES
###########################################################################
def parse_so_energy(fout, SFS):
    """
    Parses spin-orbit energies from a MOLPRO log file.

    Parameters:
        fout (str): Path to the MOLPRO log file.
        SFS (dict): Number of spin-free states per symmetry in the calculation.

    Returns:
        numpy.ndarray: Array containing spin-orbit energies.
    """
    try:
        with open(fout, 'r') as f:
            lines = f.readlines()
            nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))
            so_energy = np.zeros(nSOS, dtype=float)

            for i, line in enumerate(lines):
                if "Spin-orbit eigenstates   (energies)" in line:
                    for j in range(nSOS):
                        so_energy[j] = float(lines[i + 5 + j].split()[1])
                    break

        return so_energy
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred while parsing spin-orbit energies: {str(e)}")
        return None

###########################################################################
#  PARSE SPIN-FREE PROPERTIES
###########################################################################

def parse_sf_prop(fout, sf_prop_type, nSFS):
    """
    Parses the spin-free properties matrix from a MOLPRO log file.

    Parameters:
        fout (str): The path to the log file containing the spin-free properties data.
        sf_prop_type (str): The type of property (e.g., 'Dipole', 'Quadrupole') to be parsed.
        nSOS (int): The number of spin-free states in the calculation.

    Returns:
        tuple: A tuple containing two NumPy arrays representing the real and imaginary parts
        of the parsed spin-free properties matrix. If the file is not found, returns (None, None).

    Raises:
        FileNotFoundError: If the specified log file is not found.
        Exception: If there is an error during parsing.
    """
    try:
        with open(fout, 'r') as ofile:
            lines = ofile.readlines()
            i, k = 0, 0

            # Initialize NumPy arrays to store real and imaginary parts of spin-free properties matrix
            sf_prop_real = np.zeros((nSFS, nSFS), dtype=float)
            sf_prop_imag = np.zeros((nSFS, nSFS), dtype=float)

            # Loop through the lines in the file to find and parse spin-free properties matrix
            while i < len(lines):
                if f"Property matrix for the {sf_prop_type} operator" in lines[i]:
                    i += 3
                    blocks = math.ceil(nSFS / 8)
                    for _ in range(blocks):
                        # Check for the presence of relevant data in the current line
                        if "Nr  Nr" in lines[i]:
                            nCols = int(len(lines[i].split())) - 4
                            for j in range(nSFS):
                                # Extract values from the line and store them in the arrays
                                line = lines[(i + 1) + j].strip()
                                line_split = line.split()
                                line_values = [float(x) for x in line_split[-nCols:]]
                                for m in range(len(line_values)):
                                    sf_prop_real[j, k + m] = line_values[m]
                                    sf_prop_imag[j, k + m] = 0.0  # Imaginary part is always 0 for spin-free properties

                            k += nCols
                            i += nSFS + 2
                        else:
                            break  # Exit the loop if "Nr  Nr" is not found
                else:
                    i += 1

        # Return the parsed spin-free properties matrix as NumPy arrays
        return sf_prop_real, sf_prop_imag

    # Handle the case where the specified file is not found
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None, None  # Return None if file not found

    # Handle any other exceptions that might occur during parsing
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is a parsing error

###########################################################################
#  PARSE SPIN-ORBIT PROPERTIES
###########################################################################
def parse_so_prop(fout, so_prop_type, SFS):
    """
    Parses the spin-orbit properties matrix from a MOLPRO log file.

    Parameters:
        fout (str): The path to the log file containing the spin-orbit properties data.
        so_prop_type (str): The type of property (e.g., 'Dipole', 'Quadrupole') to be parsed.
        nSOS (int): The number of spin-orbit states in the calculation.

    Returns:
        tuple: A tuple containing two NumPy arrays representing the real and imaginary parts
        of the parsed spin-orbit properties matrix. If the file is not found, returns (None, None).

    Raises:
        FileNotFoundError: If the specified log file is not found.
        Exception: If there is an error during parsing.
    """
    try:
        with open(fout, 'r') as ofile:
            lines = ofile.readlines()
           #i, k = 0, 0

            # total number of spin-free states (SFS)
            nSFS = sum(SFS.values())

            # SOS = total number of spin-orbit states (SOS)
            nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

            # Initialize NumPy arrays to store real and imaginary parts of spin-orbit properties matrix
            so_prop_real = np.zeros((nSOS, nSOS), dtype=float)
            so_prop_imag = np.zeros((nSOS, nSOS), dtype=float)

            blocks = math.ceil(nSOS / 8)
            # Loop through the lines in the file to find and parse spin-orbit properties matrix
            i = 0
            while i < len(lines):

                if f"{so_prop_type} (TRANSFORMED, REAL)" in lines[i]:
                    i += 1
                    k = 0
                    for _ in range(blocks):
                        nCols = int(len(lines[i].split()))
                        for j in range(nSOS):
                            real_line = lines[(i + 1) + j].strip().split()
                           #print(real_line)
                            real_vals = [float(x) for x in real_line[-nCols:]]
                            for m in range(len(real_vals)):
                                so_prop_real[j, k + m] = real_vals[m]
                        k += nCols
                        i += nSOS + 1

                elif f"{so_prop_type} (TRANSFORMED, IMAG)" in lines[i]:
                    i += 1
                    k = 0
                    for _ in range(blocks):
                        nCols = int(len(lines[i].split()))
                        for j in range(nSOS):
                            imag_line = lines[(i + 1) + j].strip().split()
                           #print(real_line)
                            imag_vals = [float(x) for x in imag_line[-nCols:]]
                            for m in range(len(imag_vals)):
                                so_prop_imag[j, k + m] = imag_vals[m]
                        k += nCols
                        i += nSOS + 1
                else:
                    i += 1

        # Return the parsed spin-orbit properties matrix as NumPy arrays
        return so_prop_real, so_prop_imag

    # Handle the case where the specified file is not found
    except FileNotFoundError:
        print(f"Error: File '{log_file}' not found.")
        return None, None  # Return None if file not found

    # Handle any other exceptions that might occur during parsing
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is a parsing error

###########################################################################
#  PARSE SPIN-ORBIT EIGENVECTORS
###########################################################################

def parse_so_eigvecs(fout, SFS):
    """
    Parses eigenvectors from a MOLPRO file.out generated by a QM calculation tool.
    
    Parameters:
        fout (str): The path to the MOLPRO file.out containing the eigenvector data.
        SFS (dictionary): The number of spin-free states (SFS) per symm in the calculation.
        
    Returns:
        tuple: A tuple containing two NumPy arrays representing the real and imaginary parts
        of the parsed eigenvectors. If the file is not found, returns (None, None).
        
    Raises:
        FileNotFoundError: If the specified log file is not found.
        Exception: If there is an error during parsing.
    """
    try:
        # Attempt to open the specified log file in read mode
        with open(fout, 'r') as ofile:
            # Read all lines from the file
            lines = ofile.readlines()
            i, k = 0, 0

            # total number of spin-free states (SFS)
            nSFS = sum(SFS.values())

            # SOS = total number of spin-orbit states (SOS)
            nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

            # Initialize NumPy arrays to store real and imag parts of eigenvectors
            eigvecs_real = np.zeros((nSOS, nSOS), dtype=float)
            eigvecs_imag = np.zeros((nSOS, nSOS), dtype=float)

            # Calculate the no. of blocks in MOLPRO file.out, based on nSOS.
            # Typical printing pattern of MOLPRO: 8 columns are printed.
            blocks = math.ceil(nSOS / 8)
            
            # Loop through the lines in the file to find and parse eigenvectors
            while i < len(lines):
                if "Eigenvectors of spin-orbit matrix" in lines[i]:
                    i += 5
                    for _ in range(blocks):
                        # Check for the presence of relevant data in the current line
                        if "Nr   State  S   Sz" in lines[i]:
                            nCols = int(len(lines[i].split())) - 4
                            for j in range(nSOS):
                                # Extract real and imag parts of eigenvectors from the lines
                                real_line = lines[(i + 2) + 3 * j + 0].strip().split()
                                imag_line = lines[(i + 2) + 3 * j + 1].strip().split()

                                # Convert extracted values to float and store them in the arrays
                                real_vals = [float(x) for x in real_line[-nCols:]]
                                imag_vals = [float(x) for x in imag_line[-nCols:]]
                                for m in range(len(real_vals)):
                                    eigvecs_real[j, k + m] = real_vals[m]
                                for m in range(len(imag_vals)):
                                    eigvecs_imag[j, k + m] = imag_vals[m]

                            k += nCols
                            i += 3 * nSOS + 3
                        else:
                            break  # Exit the loop if "Nr Sym" is not found
                else:
                    i += 1

            #####################################################
            # The way SO eigenvectors are parsed in MOLPRO format
            #####################################################
            molpro_fmt_eigvecs_real = np.zeros((nSOS, nSOS), dtype=float)
            molpro_fmt_eigvecs_imag = np.zeros((nSOS, nSOS), dtype=float)
            molpro_fmt_eigvecs_real = eigvecs_real.copy()
            molpro_fmt_eigvecs_imag = eigvecs_imag.copy()
            
            #########################################################
            # Rearrange basis of the SO eigenvectors in MOLCAS format
            #########################################################
            molcas_fmt_eigvecs_real = np.zeros((nSOS, nSOS), dtype=float)
            molcas_fmt_eigvecs_imag = np.zeros((nSOS, nSOS), dtype=float)
            
            # Example for: 2 Irreps, 1 state in Irrep(1) and 2 states in Irrep(2), and Mult = 3
            # =====================    ====================
            #       MOLPRO(EigVec)          MOLCAS(EigVec)
            # =====================    ====================
            #  i    Irrep, Ms, SFS      j   Irrep, SFS, Ms
            # =====================    ====================
            #  0   (    1, +1,   1)     0  (    1,   1, -1)
            #  1   (    1,  0,   1)     1  (    1,   1,  0)
            #  2   (    1, -1,   1)     2  (    1,   1, +1)
            #  3   (    2, +1,   1)     3  (    2,   1, -1)
            #  4   (    2, +1,   2)     4  (    2,   1,  0)
            #  5   (    2,  0,   1)     5  (    2,   1, +1)
            #  6   (    2,  0,   2)     6  (    2,   2, -1)
            #  7   (    2, -1,   1)     7  (    2,   2,  0)
            #  8   (    2, -1,   2)     8  (    2,   2, +1)
            # =====================    ====================

            i = 0 # MOLPRO idx
            for iIrrep in range(1, nIrrep + 1): # MOLPRO 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
                iS = (MULT[iIrrep] - 1) / 2
                for iSz in np.arange(iS, -(iS + 1), -1): # MOLPRO 2nd loop: +Ms, ..., 0, ..., -Ms
                    for iSFS in range(1, (SFS[iIrrep] + 1)): # MOLPRO 3rd loop: SFS(1), SFS(2), ..., SFS(n)

                        j = 0 # MOLCAS idx
                        for jIrrep in range(1, nIrrep + 1): # MOLCAS 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
                            for jSFS in range(1, (SFS[jIrrep] + 1)): # MOLCAS 2nd loop: SFS(1), SFS(2), ..., SFS(n)
                                jS = (MULT[iIrrep] - 1) / 2
                                for jSz in np.arange(-jS, (jS + 1), 1): # MOLCAS 3rd loop: -Ms, ..., 0, ..., +Ms

                                    if iIrrep == jIrrep and iSFS == jSFS and iSz == jSz:
                                        molcas_fmt_eigvecs_real[j,:] = molpro_fmt_eigvecs_real[i,:]
                                        molcas_fmt_eigvecs_imag[j,:] = molpro_fmt_eigvecs_imag[i,:]

                                    j += 1
                        i += 1 

            # Return the parsed eigenvectors as NumPy arrays
            return molcas_fmt_eigvecs_real, molcas_fmt_eigvecs_imag
           #return molpro_fmt_eigvecs_real, molpro_fmt_eigvecs_imag

    # Handle the case where the specified file is not found
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None, None  # Return None if file not found

    # Handle any other exceptions that might occur during parsing
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is a parsing error

###########################################################################
#  SFS-to-SOS MATRIX EXPANSION -> FOLLOWED BY UNIRATY TRANSFORMATION
###########################################################################
def sf2so_prop_trafo(U_real, U_imag, sf_prop, nIrrep, SFS, MULT):
    """
    Transform single-particle (sf) properties to single-occupation (so) basis
    using provided transformation matrices and spin multiplicity.

    Parameters:
    - U_real (numpy.ndarray): Real part of the transformation matrix from sf to so basis.
    - U_imag (numpy.ndarray): Imaginary part of the transformation matrix from sf to so basis.
    - sf_prop (numpy.ndarray): 2D array representing sf properties in sf basis.
    - sf_states (int): Number of single-particle states.
    - spin_mult (int): Spin multiplicity of the system.

    Returns:
    - sf2so_prop_real (numpy.ndarray): Real part of the transformed so properties.
    - sf2so_prop_imag (numpy.ndarray): Imaginary part of the transformed so properties.

    If there is an error during transformation, the function prints an error message and returns None for both matrices.
    """
    try:
        # total number of spin-free states (SFS)
        nSFS = sum(SFS.values())
        
        # SOS = total number of spin-orbit states (SOS)
        nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

        # Create an empty expanded matrix filled with zeros
        sf_prop_expand = np.zeros((nSOS, nSOS), dtype=float)

       ############################################################
       ## MOLPRO-style expansion of Spin-free to spin-orbit matrix
       ############################################################
       #i = 0
       #for iIrrep in range(1, nIrrep + 1): # MOLPRO 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
       #    iS = (MULT[iIrrep] - 1) / 2             
       #    for iSz in np.arange(iS, -(iS + 1), -1): # MOLPRO 2nd loop: +Ms, ..., 0, ..., -Ms
       #        j = 0
       #        for jIrrep in range(1, nIrrep + 1): # MOLPRO 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
       #            jS = (MULT[iIrrep] - 1) / 2
       #            for jSz in np.arange(jS, -(jS + 1), -1): # MOLPRO 2nd loop: +Ms, ..., 0, ..., -Ms
       #                row_start = i * SFS[iIrrep]
       #                col_start = j * SFS[jIrrep]
       #                m = 0
       #                for iSFS in range(1, (SFS[iIrrep] + 1)): # MOLPRO 3rd loop: Irrep(1), Irrep(2), ..., Ireep(n)
       #                    n = 0
       #                    for jSFS in range(1, (SFS[jIrrep] + 1)): # MOLPRO 3rd loop: Irrep(1), Irrep(2), ..., Ireep(n)
       #                        if iSz == jSz:
       #                            sf_prop_expand[row_start + m, col_start + n] = sf_prop[i, j]
       #                        n += 1
       #                    m += 1
       #                j += 1
       #        i += 1
       ############################################################

        ###########################################################
        # MOLCAS-style expansion of Spin-free to spin-orbit matrix
        ###########################################################
        i = 0 
        for iIrrep in range(1, nIrrep + 1): # MOLCAS 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
            for iSFS in range(1, SFS[iIrrep] + 1): # MOLCAS 2nd loop: SFS(1), SFS(2), ..., SFS(n)
                j = 0
                for jIrrep in range(1, nIrrep + 1): # MOLCAS 1st loop: Irrep(1), Irrep(2), ..., Ireep(n)
                    for jSFS in range(1, SFS[jIrrep] + 1): # MOLCAS 2nd loop: SFS(1), SFS(2), ..., SFS(n)
                        row_start = i * MULT[iIrrep]
                        col_start = j * MULT[jIrrep]
                        iS = (MULT[iIrrep] - 1) / 2
                        jS = (MULT[jIrrep] - 1) / 2
                        m = 0
                        for iSz in np.arange(-iS, (iS + 1),  1): # MOLCAS 3rd loop: -Ms, ..., 0, ..., +Ms
                            n = 0
                            for jSz in np.arange(-jS, (jS + 1),  1): # MOLCAS 3rd loop: -Ms, ..., 0, ..., +Ms
                                if iSz == jSz:
                                    sf_prop_expand[row_start + m, col_start + n] = sf_prop[i, j]
                                n += 1
                            m += 1
                        j += 1
                i += 1
        ###########################################################

        interim_real = np.zeros((nSOS, nSOS), dtype=float)
        interim_imag = np.zeros((nSOS, nSOS), dtype=float)

        # Calculate intermediate matrices
        for iSOS in range(nSOS):
            for jSOS in range(nSOS):
                for kSOS in range(nSOS):
                    interim_real[iSOS][jSOS] += U_real[kSOS][iSOS] * sf_prop_expand[kSOS][jSOS]
                    interim_imag[iSOS][jSOS] -= U_imag[kSOS][iSOS] * sf_prop_expand[kSOS][jSOS]

        sf2so_prop_real = np.zeros((nSOS, nSOS), dtype=float)
        sf2so_prop_imag = np.zeros((nSOS, nSOS), dtype=float)

        # Calculate final so matrix
        for iSOS in range(nSOS):
            for jSOS in range(nSOS):
                for kSOS in range(nSOS):
                    sf2so_prop_real[iSOS][jSOS] += (interim_real[iSOS][kSOS] * U_real[kSOS][jSOS] - interim_imag[iSOS][kSOS] * U_imag[kSOS][jSOS])
                    sf2so_prop_imag[iSOS][jSOS] += (interim_real[iSOS][kSOS] * U_imag[kSOS][jSOS] + interim_imag[iSOS][kSOS] * U_real[kSOS][jSOS])

       #return sf_prop_expand, sf2so_prop_real, sf2so_prop_imag
        return sf2so_prop_real, sf2so_prop_imag

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is an error during transformation


############################################################################
# Boltzman distribution of Ground States split under SO coupling
############################################################################
def boltzmann_distro(rel_energy_au, st_idx, kB = 3.1668105e-6, T = 300):
    """
    Calculate the Boltzmann factor for a specific state index within sodgen states.

    Parameters:
        rel_energy_sodgen_au (numpy.ndarray): Array containing the energy differences of each spin-orbit state (in atomic units).
        state_index (int): The index of the state to calculate the Boltzmann factor for.
        T (float, optional): The temperature in Kelvin. Default is 300 K.

    Returns:
        float: The Boltzmann factor for the specified state index.
    """
    # Calculate beta
    beta = 1 / (kB * T)

    # Calculate Boltzmann factor for the specified state
    boltz_fac = np.exp(-beta * rel_energy_au[st_idx])

    # Normalize Boltzmann factor
    boltz_fac /= np.sum(np.exp(-beta * rel_energy_au))

    return boltz_fac

############################################################################
# Group the nearly degenerate states
############################################################################
def group_and_sum(df):
    """
    Group values in Column1 if they differ by less than 5 and sum corresponding values in Column2.

    Args:
        df (DataFrame): Input DataFrame with two columns - Column1 and Column2.

    Returns:
        DataFrame: DataFrame with grouped values in Column1 and summed values in Column2.
    """
    result = []
    current_group = []
    current_sum = 0

    for i, row in df.iterrows():
        if not current_group or abs(current_group[-1] - row['Column1']) <= 20:
            current_group.append(row['Column1'])
            current_sum += row['Column2']
        else:
            result.append({'Column1': sum(current_group) / len(current_group), 'Column2': current_sum})
            current_group = [row['Column1']]
            current_sum = row['Column2']

    if current_group:
        result.append({'Column1': sum(current_group) / len(current_group), 'Column2': current_sum})

    return pd.DataFrame(result)



############################################################################
############################################################################
############################################################################


############################################################################
# THE MAIN PROGRAM
############################################################################
def main():
    """
    Main function to parse and process spin-orbit properties matrices from MOLPRO output.
    Writes the processed data into output files for further analysis.
    """
    try:
        sf_energy = parse_sf_energy(fout, SFS, nSFS)
        so_energy = parse_so_energy(fout, SFS)
        sf_dmx_real, sf_dmx_imag = parse_sf_prop(fout, "DMX", nSFS)
        sf_dmy_real, sf_dmy_imag = parse_sf_prop(fout, "DMY", nSFS)
        sf_dmz_real, sf_dmz_imag = parse_sf_prop(fout, "DMZ", nSFS)
        so_dmx_real, so_dmx_imag = parse_so_prop(fout, "DMX", SFS)
        so_dmy_real, so_dmy_imag = parse_so_prop(fout, "DMY", SFS)
        so_dmz_real, so_dmz_imag = parse_so_prop(fout, "DMZ", SFS)
        eigvecs_real, eigvecs_imag = parse_so_eigvecs(fout, SFS)
        sf2so_dmx_real, sf2so_dmx_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmx_real, nIrrep, SFS, MULT)
        sf2so_dmy_real, sf2so_dmy_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmy_real, nIrrep, SFS, MULT)
        sf2so_dmz_real, sf2so_dmz_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmz_real, nIrrep, SFS, MULT)

#      ###########################################################################
#      ## DEBUG: PRINT SF-DMX, SF-DMY, SF-DMZ 
#      ###########################################################################
#      #with open('sf-dmx.txt', 'w') as sfdmx, \
#      #     open('sf-dmy.txt', 'w') as sfdmy, \
#      #     open('sf-dmz.txt', 'w') as sfdmz:
#      #     for iSFS in range(nSFS):
#      #        for jSFS in range(nSFS):
#      #            sfdmx.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmx_real[iSFS, jSFS]:16.8E}  {sf_dmx_imag[iSFS, jSFS]:16.8E}\n") 
#      #            sfdmy.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmy_real[iSFS, jSFS]:16.8E}  {sf_dmy_imag[iSFS, jSFS]:16.8E}\n") 
#      #            sfdmz.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmz_real[iSFS, jSFS]:16.8E}  {sf_dmz_imag[iSFS, jSFS]:16.8E}\n") 
#      ###########################################################################
#      ## END DEBUG
#      ###########################################################################

#      ###########################################################################
#      ## DEBUG: PRINT SO-DMX, SO-DMY, SO-DMZ and SF2SO-DMX, SF2SO-DMY, SF2SO-DMZ
#      ###########################################################################
#      #with open('so-dmx.txt', 'w') as sodmx, \
#      #     open('so-dmy.txt', 'w') as sodmy, \
#      #     open('so-dmz.txt', 'w') as sodmz, \
#      #     open('sf2so-dmx.txt', 'w') as sf2sodmx, \
#      #     open('sf2so-dmy.txt', 'w') as sf2sodmy, \
#      #     open('sf2so-dmz.txt', 'w') as sf2sodmz:
#      #     for iSOS in range(nSOS):
#      #        for jSOS in range(nSOS):
#      #            sodmx.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmx_real[iSOS, jSOS]:16.8E}  {so_dmx_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sodmy.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmy_real[iSOS, jSOS]:16.8E}  {so_dmy_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sodmz.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmz_real[iSOS, jSOS]:16.8E}  {so_dmz_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmx.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmx_real[iSOS, jSOS]:16.8E}  {sf2so_dmx_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmy.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmy_real[iSOS, jSOS]:16.8E}  {sf2so_dmy_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmz.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmz_real[iSOS, jSOS]:16.8E}  {sf2so_dmz_imag[iSOS, jSOS]:16.8E}\n") 
#      ###########################################################################
#      ## END DEBUG
#      ###########################################################################

        #############################################
        # SF-ABSRPTION SPECTRA (DIRECT FROM SF-TDMs)
        #############################################
        rel_energy_au = np.zeros(nSFS, dtype=float)
        rel_energy_ev = np.zeros(nSFS, dtype=float)
        fosc          = np.zeros(nSFS, dtype=float)

        for jSFS in range(sfdgen):
            for iSFS in range(nSFS):
                rel_energy_au[iSFS] = sf_energy[iSFS] - sf_energy[0]
                rel_energy_ev[iSFS] = rel_energy_au[iSFS] * 27.2114
                fosc[iSFS] = (1/sfdgen) * (2/3) * (rel_energy_au[iSFS]) * \
                             (sf_dmx_real[iSFS,jSFS]**2 + sf_dmx_imag[iSFS,jSFS]**2 + \
                              sf_dmy_real[iSFS,jSFS]**2 + sf_dmy_imag[iSFS,jSFS]**2 + \
                              sf_dmz_real[iSFS,jSFS]**2 + sf_dmz_imag[iSFS,jSFS]**2)

        with open('sf-osc-ev.stk', 'w') as f:
            for jSFS in range(sfdgen):
                for iSFS in range(nSFS):
                    f.write(f"{rel_energy_ev[iSFS]:16.4f} {fosc[iSFS]:16.2E}\n")

        ##############################################
        # BOLTZMANN POPULATIONS OF THE SO GROUND STATES
        ##############################################
        rel_energy_sodgen_au = np.zeros(sodgen, dtype=float)
        for jSOS in range(sodgen):
            rel_energy_sodgen_au[jSOS] = so_energy[jSOS] - so_energy[0]    

        with open('boltzmann.txt', 'w') as f:
            for jSOS in range(sodgen):
                boltz_fac = boltzmann_distro(rel_energy_sodgen_au, jSOS)
                f.write(f"{jSOS : 8d} {boltz_fac : 16.4f}\n")
       
        ##############################################
        # SO-ABSRPTION SPECTRA (DIRECT FROM SO-TDMs)
        ##############################################
        rel_energy_au = np.zeros(nSOS, dtype=float)
        rel_energy_ev = np.zeros(nSOS, dtype=float)
        fosc          = np.zeros(nSOS, dtype=float)

        for jSOS in range(sodgen):
            boltz_fac = boltzmann_distro(rel_energy_sodgen_au, jSOS)
            for iSOS in range(nSOS):
                rel_energy_au[iSOS] = so_energy[iSOS] - so_energy[0]
                rel_energy_ev[iSOS] = rel_energy_au[iSOS] * 27.2114
                fosc[iSOS] = (boltz_fac) * (2/3) * (rel_energy_au[iSOS]) * \
                             (so_dmx_real[iSOS,jSOS]**2 + so_dmx_imag[iSOS,jSOS]**2 + \
                              so_dmy_real[iSOS,jSOS]**2 + so_dmy_imag[iSOS,jSOS]**2 + \
                              so_dmz_real[iSOS,jSOS]**2 + so_dmz_imag[iSOS,jSOS]**2)

        with open('so-osc-ev.stk', 'w') as f:
            for jSOS in range(sodgen):
                for iSOS in range(nSOS):
                    f.write(f"{rel_energy_ev[iSOS]:16.4f} {fosc[iSOS]:16.2E}\n")

        ##############################################
        # SO-ABSRPTION SPECTRA (FROM SF-2-SO-TDMs)
        ##############################################
        rel_energy_au = np.zeros(nSOS, dtype=float)
        rel_energy_ev = np.zeros(nSOS, dtype=float)
        fosc          = np.zeros(nSOS, dtype=float)

        for jSOS in range(sodgen):
            boltz_fac = boltzmann_distro(rel_energy_sodgen_au, jSOS)
            for iSOS in range(nSOS):
                rel_energy_au[iSOS] = so_energy[iSOS] - so_energy[0]
                rel_energy_ev[iSOS] = rel_energy_au[iSOS] * 27.2114
                fosc[iSOS] = (boltz_fac) * (2/3) * (rel_energy_au[iSOS]) * \
                             (sf2so_dmx_real[iSOS,jSOS]**2 + sf2so_dmx_imag[iSOS,jSOS]**2 + \
                              sf2so_dmy_real[iSOS,jSOS]**2 + sf2so_dmy_imag[iSOS,jSOS]**2 + \
                              sf2so_dmz_real[iSOS,jSOS]**2 + sf2so_dmz_imag[iSOS,jSOS]**2)

        with open('sf2so-osc-ev.stk', 'w') as f:
            for jSOS in range(sodgen):
                for iSOS in range(nSOS):
                    f.write(f"{rel_energy_ev[iSOS]:16.4f} {fosc[iSOS]:16.2E}\n")

        ##################################################################
        # Print the bright transitions only SF
        ##################################################################
        data_sf = []
        with open('sf-osc-ev.stk', 'r') as infile:
            for line in infile:
                columns = line.split()
                data_sf.append((float(columns[0]), float(columns[1])))
        
        sorted_data_sf = sorted(data_sf, key=lambda x: x[0])
        
        grouped_data_sf = {}
        for energy_ev, oscillator in sorted_data_sf:
            if energy_ev in grouped_data_sf:
                grouped_data_sf[energy_ev] += oscillator
            else:
                grouped_data_sf[energy_ev] = oscillator
        
       #with open('sf-osc.stk', 'w') as ofile_rcm, open('sf-osc-ev.stk', 'w') as ofile_ev:
        with open('sf-osc-rcm.stk', 'w') as ofile_rcm:
            for energy_ev, oscillator in grouped_data_sf.items():
                if oscillator > 1e-3:  # Check if the oscillator value is greater than 10^-5
                    ofile_rcm.write(f"{(round(int(energy_ev * 8065.54) / 10) * 10) : 8d}  {oscillator : 16.2E}\n")
                   #ofile_ev.write(f"{int(row['Column1']):16.3f}  {row['Column2']:16.2E}\n")
        ##################################################################

        ##################################################################
        # Print the bright transitions only SO (grouped)
        ##################################################################
        data_so = []
        with open('so-osc-ev.stk', 'r') as infile:
            for line in infile:
                columns = line.split()
                data_so.append((float(columns[0]), float(columns[1])))
        
        sorted_data_so = sorted(data_so, key=lambda x: x[0])
        
        grouped_data_so = {}
        for energy_ev, oscillator in sorted_data_so:
            if energy_ev in grouped_data_so:
                grouped_data_so[energy_ev] += oscillator
            else:
                grouped_data_so[energy_ev] = oscillator
        
        with open('tmp_so.stk', 'w') as outfile:
            for energy_ev, oscillator in grouped_data_so.items():
                if oscillator > 1e-3:  # Check if the oscillator value is greater than 10^-5
                    outfile.write(f"{int(energy_ev*8065.54):16d}  {oscillator:16.2E}\n")

        # Almost degenerate states are grouped
        df = pd.read_csv('tmp_so.stk', delim_whitespace=True, header=None, names=['Column1', 'Column2'])
        result_df = group_and_sum(df)

        with open('so-osc-rcm.stk', 'w') as ofile_rcm:
            for _, row in result_df.iterrows():
                ofile_rcm.write(f"{(round(int(row['Column1']) / 10) * 10) : 8d}  {row['Column2'] : 16.2E}\n")
        os.remove('tmp_so.stk')
        ##################################################################


        ##################################################################
        # Print the bright transitions only SF2SO (grouped)
        ##################################################################
        data_sf2so = []
        with open('sf2so-osc-ev.stk', 'r') as infile:
            for line in infile:
                columns = line.split()
                data_sf2so.append((float(columns[0]), float(columns[1])))
        
        sorted_data_sf2so = sorted(data_sf2so, key=lambda x: x[0])
        
        grouped_data_sf2so = {}
        for energy_ev, oscillator in sorted_data_sf2so:
            if energy_ev in grouped_data_sf2so:
                grouped_data_sf2so[energy_ev] += oscillator
            else:
                grouped_data_sf2so[energy_ev] = oscillator
        
        with open('tmp_sf2so.stk', 'w') as outfile:
            for energy_ev, oscillator in grouped_data_sf2so.items():
                if oscillator > 1e-3:  # Check if the oscillator value is greater than 10^-5
                    outfile.write(f"{int(energy_ev*8065.54):16d}  {oscillator:16.2E}\n")

        # Almost degenerate states are grouped
        df = pd.read_csv('tmp_sf2so.stk', delim_whitespace=True, header=None, names=['Column1', 'Column2'])
        result_df = group_and_sum(df)

        with open('sf2so-osc-rcm.stk', 'w') as ofile_rcm:
            for _, row in result_df.iterrows():
                ofile_rcm.write(f"{(round(int(row['Column1']) / 10) * 10) : 8d}  {row['Column2'] : 16.2E}\n")
        os.remove('tmp_sf2so.stk')

#      ###########################################################################
#      ## DEBUG: PRINT SO-DMX, SO-DMY, SO-DMZ and SF2SO-DMX, SF2SO-DMY, SF2SO-DMZ
#      ###########################################################################
#      #with open('so-dmx.txt', 'w') as sodmx, \
#      #     open('so-dmy.txt', 'w') as sodmy, \
#      #     open('so-dmz.txt', 'w') as sodmz, \
#      #     open('sf2so-dmx.txt', 'w') as sf2sodmx, \
#      #     open('sf2so-dmy.txt', 'w') as sf2sodmy, \
#      #     open('sf2so-dmz.txt', 'w') as sf2sodmz:
#      #     for iSOS in range(nSOS):
#      #        for jSOS in range(nSOS):
#      #            sodmx.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmx_real[iSOS, jSOS]:16.8E}  {so_dmx_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sodmy.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmy_real[iSOS, jSOS]:16.8E}  {so_dmy_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sodmz.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmz_real[iSOS, jSOS]:16.8E}  {so_dmz_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmx.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmx_real[iSOS, jSOS]:16.8E}  {sf2so_dmx_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmy.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmy_real[iSOS, jSOS]:16.8E}  {sf2so_dmy_imag[iSOS, jSOS]:16.8E}\n") 
#      #            sf2sodmz.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmz_real[iSOS, jSOS]:16.8E}  {sf2so_dmz_imag[iSOS, jSOS]:16.8E}\n") 
#      ###########################################################################
#      ## END DEBUG
       ###########################################################################
        
    except Exception as e:
        print(f"An error occurred while processing spin-free energies: {str(e)}")

if __name__ == '__main__':
    main()



