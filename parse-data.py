#!/usr/bin/env python

import os
import re
import numpy as np
import math

###########################################################################
# PARAMETERS FOR DATA
###########################################################################

# Irrep = number of symmetry considered in calculation
nIrrep = 4

# SFS = Dictionary of spin-free states (SFS) per irrep
SFS = {1: 1,
       2: 1,
       3:60,
       4:60}

# total number of spin-free states (SFS)
nSFS = sum(SFS.values())

# MULT = multiplicity of each spin-free state [ 2 = triplet ]
MULT = {i: 4 for i in range(1, nIrrep + 1)}

# SOS = total number of spin-orbit states (SOS)
nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

# MOLPRO file.out
fout = "cop.nevpt-so.out"

###########################################################################
#  SPIN-ORBIT ENERGIES
###########################################################################
def parse_so_energy(fout, SFS):
    """
    Parses the spin-free properties matrix from a MOLPRO log file.

    Parameters:
        fout (str): The path to the log file containing the spin-free properties data.
        SFS (dictionary): The number of spin-free states per symm in the calculation.

    Returns:
        tuple: A tuple containing two NumPy arrays representing the SO energies.
        If the file is not found, returns (None, None).

    Raises:
        FileNotFoundError: If the specified log file is not found.
        Exception: If there is an error during parsing.
    """
    try:
        with open(fout, 'r') as ofile:
            lines = ofile.readlines()
            i = 0
       
            nSFS = sum(SFS.values())
            nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))
 
            # Initialize NumPy arrays to store real and imaginary parts of spin-free properties matrix
            so_energy = np.zeros(nSOS, dtype=float)
        
            # Loop through the lines in the file to find and parse spin-free properties matrix
            while i < len(lines):
                if "Spin-orbit eigenstates   (energies)" in lines[i]:
                    i += 4
                    for iSOS in range(nSOS):
                        so_energy[iSOS] = float(lines[(i + 1) + iSOS].strip().split()[1])
                    else:
                        break  # Exit the loop if "Nr  Nr" is not found
                else:
                    i += 1

        # Return the parsed spin-free properties matrix as NumPy arrays
        return so_energy

    # Handle the case where the specified file is not found
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None, None  # Return None if file not found

    # Handle any other exceptions that might occur during parsing
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is a parsing error

###########################################################################
#  CONSTRUCT SPIN MATRICES
###########################################################################
def make_spin_matrix(SFS, MULT):
    """
    Generate spin matrices (Sx, Sy, Sz) for a set of irreducible representations.

    Parameters:
    - SFS (list): Number of states for each irreducible representation.
    - MULT (list): Multiplicity of each irreducible representation.

    Returns:
    Tuple of NumPy arrays representing the real and imaginary parts of Sx, Sy, and Sz matrices.

    Each matrix element is computed based on quantum numbers (iIrrep, iSFS, iSz) and (jIrrep, jSFS, jSz).

    Example:
    sx_real, sx_imag, sy_real, sy_imag, sz_real, sz_imag = make_spin_matrix([2, 1], [2, 1])
    """

    nSFS = sum(SFS.values())
    nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

    sx_real = np.zeros((nSOS, nSOS), dtype=float)
    sx_imag = np.zeros((nSOS, nSOS), dtype=float)
    sy_real = np.zeros((nSOS, nSOS), dtype=float)
    sy_imag = np.zeros((nSOS, nSOS), dtype=float)
    sz_real = np.zeros((nSOS, nSOS), dtype=float)
    sz_imag = np.zeros((nSOS, nSOS), dtype=float)

    # Sx (REAL) MATRIX
    i = 0
    for iIrrep in range(1, nIrrep + 1):
        for iSFS in range(1, SFS[iIrrep] + 1):
            iS = (MULT[iIrrep] - 1) / 2
            for iSz in np.arange(-iS, (iS + 1), 1):
                j = 0
                for jIrrep in range(1, nIrrep + 1):
                    for jSFS in range(1, SFS[jIrrep] + 1):  # Fixed this line
                        jS = (MULT[jIrrep] - 1) / 2
                        for jSz in np.arange(-jS, (jS + 1), 1):

                            if iIrrep == jIrrep and iSFS == jSFS and iSz == (jSz + 1):
                                sx_real[i, j] = 0.5 * np.sqrt(jS * (jS + 1) - jSz * (jSz + 1))
                            elif iIrrep == jIrrep and iSFS == jSFS and iSz == (jSz - 1):
                                sx_real[i, j] = 0.5 * np.sqrt(jS * (jS + 1) - jSz * (jSz - 1))
                            j += 1
                i += 1

    # Sy (IMAG) MATRIX
    i = 0
    for iIrrep in range(1, nIrrep + 1):
        for iSFS in range(1, SFS[iIrrep] + 1):
            iS = (MULT[iIrrep] - 1) / 2
            for iSz in np.arange(-iS, (iS + 1), 1):
                j = 0
                for jIrrep in range(1, nIrrep + 1):
                    for jSFS in range(1, SFS[jIrrep] + 1):  # Fixed this line
                        jS = (MULT[jIrrep] - 1) / 2
                        for jSz in np.arange(-jS, (jS + 1), 1):

                            if iIrrep == jIrrep and iSFS == jSFS and iSz == (jSz + 1):
                                sy_imag[i, j] =-0.5 * np.sqrt(jS * (jS + 1) - jSz * (jSz + 1))
                            elif iIrrep == jIrrep and iSFS == jSFS and iSz == (jSz - 1):
                                sy_imag[i, j] = 0.5 * np.sqrt(jS * (jS + 1) - jSz * (jSz - 1))
                            j += 1
                i += 1

    # Sz (REAL) MATRIX
    i = 0
    for iIrrep in range(1, nIrrep + 1):
        for iSFS in range(1, SFS[iIrrep] + 1):
            iS = (MULT[iIrrep] - 1) / 2
            for iSz in np.arange(-iS, (iS + 1), 1):
                j = 0
                for jIrrep in range(1, nIrrep + 1):
                    for jSFS in range(1, SFS[jIrrep] + 1):  # Fixed this line
                        jS = (MULT[jIrrep] - 1) / 2
                        for jSz in np.arange(-jS, (jS + 1), 1):

                            if iIrrep == jIrrep and iSFS == jSFS and iSz == jSz:
                                sz_real[i, j] = jSz
                            j += 1
                i += 1

    return sx_real, sx_imag, sy_real, sy_imag, sz_real, sz_imag

###########################################################################
#  SF2SO SPIN MATRIX TRANSFORMATION
###########################################################################
def sf2so_spin_mat_trafo(U_real, U_imag, spin_mat, SFS, MULT):
    """
    Transform a spin matrix from the spin-free basis to the spin-orbit basis.

    Parameters:
    - U_real (numpy.ndarray): Real part of the unitary transformation matrix.
    - U_imag (numpy.ndarray): Imaginary part of the unitary transformation matrix.
    - spin_mat (numpy.ndarray): Spin matrix in the spin-free basis.
    - SFS (dict): Dictionary representing the number of states for each irreducible representation.
    - MULT (list): List representing the multiplicity of each irreducible representation.

    Returns:
    Tuple of NumPy arrays representing the real and imaginary parts of the spin matrix
    transformed to the spin-orbit basis.

    The transformation is performed using the provided unitary transformation matrices (U_real and U_imag).

    Example:
    U_real = np.array([[1, 0], [0, 1]])
    U_imag = np.array([[0, 0], [0, 0]])
    spin_mat = np.array([[1, 0], [0, -1]])
    SFS = {1: 2}
    MULT = [2]
    sf2so_real, sf2so_imag = sf2so_spin_mat_trafo(U_real, U_imag, spin_mat, SFS, MULT)
    """
    try:
        # total number of spin-free states (SFS)
        nSFS = sum(SFS.values())

        # SOS = total number of spin-orbit states (SOS)
        nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))

        interim_real = np.zeros((nSOS, nSOS), dtype=float)
        interim_imag = np.zeros((nSOS, nSOS), dtype=float)

        # Calculate intermediate matrices
        for iSOS in range(nSOS):
            for jSOS in range(nSOS):
                for kSOS in range(nSOS):
                    interim_real[iSOS][jSOS] += U_real[kSOS][iSOS] * spin_mat[kSOS][jSOS]
                    interim_imag[iSOS][jSOS] -= U_imag[kSOS][iSOS] * spin_mat[kSOS][jSOS]

        sf2so_spin_mat_real = np.zeros((nSOS, nSOS), dtype=float)
        sf2so_spin_mat_imag = np.zeros((nSOS, nSOS), dtype=float)

        # Calculate final so matrix
        for iSOS in range(nSOS):
            for jSOS in range(nSOS):
                for kSOS in range(nSOS):
                    sf2so_spin_mat_real[iSOS][jSOS] += (interim_real[iSOS][kSOS] * U_real[kSOS][jSOS] - interim_imag[iSOS][kSOS] * U_imag[kSOS][jSOS])
                    sf2so_spin_mat_imag[iSOS][jSOS] += (interim_real[iSOS][kSOS] * U_imag[kSOS][jSOS] + interim_imag[iSOS][kSOS] * U_real[kSOS][jSOS])

        return sf2so_spin_mat_real, sf2so_spin_mat_imag

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is an error during transformation

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
        
        # MOLCAS-style expansion
        i = 0
        for iIrrep in range(1, nIrrep + 1):
            for iSFS in range(1, SFS[iIrrep] + 1):
                j = 0
                for jIrrep in range(1, nIrrep + 1):
                    for jSFS in range(1, SFS[jIrrep] + 1):
                       #print("i = ", i, "j = ", j)
                        row_start = i * MULT[iIrrep]
                        col_start = j * MULT[jIrrep]
                        iSz = (MULT[iIrrep] - 1) / 2
                        jSz = (MULT[jIrrep] - 1) / 2
                        ikk = 0
                        for ik in np.arange(-iSz, (iSz + 1),  1):
                            jkk = 0
                            for jk in np.arange(-jSz, (jSz + 1),  1):
                               #if ikk == jkk:
                                if iSz == jSz:
                                    sf_prop_expand[row_start + ikk, col_start + jkk] = sf_prop[i, j]
                                jkk += 1
                            ikk += 1
                        j += 1
                i += 1

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

            # MOLPRO and MOLCAS eigenvectors
           #molpro_fmt_eigvecs_real = np.zeros((nSOS, nSOS), dtype=float)
           #molpro_fmt_eigvecs_imag = np.zeros((nSOS, nSOS), dtype=float)

            molcas_fmt_eigvecs_real = np.zeros((nSOS, nSOS), dtype=float)
            molcas_fmt_eigvecs_imag = np.zeros((nSOS, nSOS), dtype=float)

           #molpro_fmt_eigvecs_real = eigvecs_real.copy()
           #molpro_fmt_eigvecs_imag = eigvecs_imag.copy()
            
            i = 0
            for iIrrep in range(1, nIrrep + 1):
                iS = (MULT[iIrrep] - 1) / 2
                for iSz in np.arange(iS, -(iS + 1), -1):
                    for iSFS in range(1, (SFS[iIrrep] + 1)):
                        j = 0
                        for jIrrep in range(1, nIrrep + 1):
                            for jSFS in range(1, (SFS[jIrrep] + 1)):
                                jS = (MULT[iIrrep] - 1) / 2
                                for jSz in np.arange(-jS, (jS + 1), 1):
                                    if iIrrep == jIrrep and iSFS == jSFS and iSz == jSz:
                                        molcas_fmt_eigvecs_real[j,:] = eigvecs_real[i,:]
                                        molcas_fmt_eigvecs_imag[j,:] = eigvecs_imag[i,:]
                                    j += 1
                        i += 1 

            # Return the parsed eigenvectors as NumPy arrays
           #return molpro_fmt_eigvecs_real, molpro_fmt_eigvecs_imag, molcas_fmt_eigvecs_real, molcas_fmt_eigvecs_imag
            return molcas_fmt_eigvecs_real, molcas_fmt_eigvecs_imag

    # Handle the case where the specified file is not found
    except FileNotFoundError:
        print(f"Error: File '{fout}' not found.")
        return None, None  # Return None if file not found

    # Handle any other exceptions that might occur during parsing
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None, None  # Return None if there is a parsing error

###########################################################################
#  ABSRPTION SPECTRA
###########################################################################
#   def absorption(dgen, so_energy, sf2so_dmx_real, sf2so_dmx_imag, sf2so_dmy_real, sf2so_dmy_imag, sf2so_dmz_real, sf2so_dmz_imag):
#       try:
#           nSFS = sum(SFS.values())
#           nSOS = sum(SFS[i] * MULT[i] for i in range(1, nIrrep + 1))
#   
#           rel_energy_au = np.zeros(nSOS, dtype=float)
#           rel_energy_ev = np.zeros(nSOS, dtype=float)
#           fosc          = np.zeros(nSOS, dtype=float)
#           
#           for i in range(nSOS):
#               rel_energy_au[i] = so_energy[i] - so_energy[0]
#               rel_energy_ev[i] = rel_energy_au[i] * 27.2114
#               fosc[i] = 2/3 * rel_energy_au[i] * (sf2so_dmx_real[i,j]**2 + sf2so_dmx_imag[i,j]**2 + \
#                                                           sf2so_dyx_real[i,j]**2 + sf2so_dyx_imag[i,j]**2 + \
#                                                           sf2so_dzx_real[i,j]**2 + sf2so_dzx_imag[i,j]**2)
#   
#           return rel_energy_ev, fosc
#   
#       # Handle the case where the specified file is not found
#       except FileNotFoundError:
#           print(f"Error: File '{fout}' not found.")
#           return None, None  # Return None if file not found
#   
#       # Handle any other exceptions that might occur during parsing
#       except Exception as e:
#           print(f"An error occurred: {str(e)}")
#           return None, None  # Return None if there is a parsing error








###########################################################################
# THE MAIN PROGRAM, CALLING ALL FUNCTIONS
###########################################################################

def main():
    """
    Main function to parse and process spin-free and spin-orbit properties matrices from MOLPRO output.
    Writes the processed data into output files for further analysis.

    Raises:
        Exception: If there is an error during the parsing or processing of the data.
    """
    try:
        #################################################
        # Parse SO energies
        #################################################
        so_energy = parse_so_energy(fout, SFS)

        # Write spin-free properties to output files
        with open('energies.txt', 'w') as f:
            # Write header to output files
            f.write(f"{nSOS} (atomic units)\n")
            # Write data to output files
            for iSOS in range(nSOS):
                f.write(f"{so_energy[iSOS]:14.8f}\n")

        #################################################
        # Create Spin matrices
        #################################################
        sf_sx_real, sf_sx_imag, sf_sy_real, sf_sy_imag, sf_sz_real, sf_sz_imag = make_spin_matrix(SFS, MULT)

       ## Write spin-free properties to output files
       #with open('sf-spin-1.txt', 'w') as f1, \
       #     open('sf-spin-2.txt', 'w') as f2, \
       #     open('sf-spin-3.txt', 'w') as f3:
       #    # Write header to output files
       #    f1.write(" #NROW NCOL REAL IMAG\n")
       #    f2.write(" #NROW NCOL REAL IMAG\n")
       #    f3.write(" #NROW NCOL REAL IMAG\n")
       #    # Write data to output files
       #    for iSOS in range(nSOS):
       #        for jSOS in range(nSOS):
       #            f1.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf_sx_real[iSOS, jSOS]:16.8E}  {sf_sx_imag[iSOS, jSOS]:16.8E}\n")
       #            f2.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf_sy_real[iSOS, jSOS]:16.8E}  {sf_sy_imag[iSOS, jSOS]:16.8E}\n")
       #            f3.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf_sz_real[iSOS, jSOS]:16.8E}  {sf_sz_imag[iSOS, jSOS]:16.8E}\n")

        #################################################
        # Parse spin-free properties
        #################################################
        sf_dmx_real, sf_dmx_imag = parse_sf_prop(fout, "DMX", nSFS)
        sf_dmy_real, sf_dmy_imag = parse_sf_prop(fout, "DMY", nSFS)
        sf_dmz_real, sf_dmz_imag = parse_sf_prop(fout, "DMZ", nSFS)
        sf_lx_real, sf_lx_imag = parse_sf_prop(fout, "LX", nSFS)
        sf_ly_real, sf_ly_imag = parse_sf_prop(fout, "LY", nSFS)
        sf_lz_real, sf_lz_imag = parse_sf_prop(fout, "LZ", nSFS)

       ## Write spin-free properties to output files
       #with open('sf-dipole-1.txt', 'w') as f1, \
       #     open('sf-dipole-2.txt', 'w') as f2, \
       #     open('sf-dipole-3.txt', 'w') as f3, \
       #     open('sf-angmom-1.txt', 'w') as f4, \
       #     open('sf-angmom-2.txt', 'w') as f5, \
       #     open('sf-angmom-3.txt', 'w') as f6:
       #    # Write header to output files
       #    f1.write(" #NROW NCOL REAL IMAG\n")
       #    f2.write(" #NROW NCOL REAL IMAG\n")
       #    f3.write(" #NROW NCOL REAL IMAG\n")
       #    f4.write(" #NROW NCOL REAL IMAG\n")
       #    f5.write(" #NROW NCOL REAL IMAG\n")
       #    f6.write(" #NROW NCOL REAL IMAG\n")
       #    # Write data to output files
       #    for iSFS in range(nSFS):
       #        for jSFS in range(nSFS):
       #            f1.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmx_real[iSFS, jSFS]:16.8E}  {sf_dmx_imag[iSFS, jSFS]:16.8E}\n")
       #            f2.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmy_real[iSFS, jSFS]:16.8E}  {sf_dmy_imag[iSFS, jSFS]:16.8E}\n")
       #            f3.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_dmz_real[iSFS, jSFS]:16.8E}  {sf_dmz_imag[iSFS, jSFS]:16.8E}\n")
       #            f4.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_lx_real[iSFS, jSFS]:16.8E}  {sf_lx_imag[iSFS, jSFS]:16.8E}\n")
       #            f5.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_ly_real[iSFS, jSFS]:16.8E}  {sf_ly_imag[iSFS, jSFS]:16.8E}\n")
       #            f6.write(f"{iSFS + 1:4d}  {jSFS + 1:2d}  {sf_lz_real[iSFS, jSFS]:16.8E}  {sf_lz_imag[iSFS, jSFS]:16.8E}\n")

        #################################################
        # Parse spin-orbit properties
        #################################################
       #so_dmx_real, so_dmx_imag = parse_so_prop(fout, "DMX", SFS)
       #so_dmy_real, so_dmy_imag = parse_so_prop(fout, "DMY", SFS)
       #so_dmz_real, so_dmz_imag = parse_so_prop(fout, "DMZ", SFS)
       #so_lx_real, so_lx_imag = parse_so_prop(fout, "LX", SFS)
       #so_ly_real, so_ly_imag = parse_so_prop(fout, "LY", SFS)
       #so_lz_real, so_lz_imag = parse_so_prop(fout, "LZ", SFS)

       ## Write spin-orbit properties to output files
       #with open('so-dipole-1.txt', 'w') as f1, \
       #     open('so-dipole-2.txt', 'w') as f2, \
       #     open('so-dipole-3.txt', 'w') as f3, \
       #     open('so-angmom-1.txt', 'w') as f4, \
       #     open('so-angmom-2.txt', 'w') as f5, \
       #     open('so-angmom-3.txt', 'w') as f6:
       #    # Write header to output files
       #    f1.write(" #NROW NCOL REAL IMAG\n")
       #    f2.write(" #NROW NCOL REAL IMAG\n")
       #    f3.write(" #NROW NCOL REAL IMAG\n")
       #    f4.write(" #NROW NCOL REAL IMAG\n")
       #    f5.write(" #NROW NCOL REAL IMAG\n")
       #    f6.write(" #NROW NCOL REAL IMAG\n")
       #    # Write data to output files
       #    for iSOS in range(nSOS):
       #        for jSOS in range(nSOS):
       #            f1.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmx_real[iSOS, jSOS]:16.8E}  {so_dmx_imag[iSOS, jSOS]:16.8E}\n")
       #            f2.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmy_real[iSOS, jSOS]:16.8E}  {so_dmy_imag[iSOS, jSOS]:16.8E}\n")
       #            f3.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_dmz_real[iSOS, jSOS]:16.8E}  {so_dmz_imag[iSOS, jSOS]:16.8E}\n")
       #            f4.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_lx_real[iSOS, jSOS]:16.8E}  {so_lx_imag[iSOS, jSOS]:16.8E}\n")
       #            f5.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_ly_real[iSOS, jSOS]:16.8E}  {so_ly_imag[iSOS, jSOS]:16.8E}\n")
       #            f6.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {so_lz_real[iSOS, jSOS]:16.8E}  {so_lz_imag[iSOS, jSOS]:16.8E}\n")

        #################################################
        # Parse SO eigenvectors
        #################################################
        eigvecs_real, eigvecs_imag = parse_so_eigvecs(fout, SFS)

        # Write MOLPRO eigenvalues and eigenvectors to output file
        with open('eigvectors.txt', 'w') as f:
            # Write header to output file
            f.write(" #NROW NCOL REAL IMAG\n")
            # Write data to output file
            for iSOS in range(nSOS):
                for jSOS in range(nSOS):
                    f.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {eigvecs_real[iSOS, jSOS]:16.8E}  {eigvecs_imag[iSOS, jSOS]:16.8E}\n")

        #################################################
        # Convert SF to SO property
        #################################################
        sf2so_dmx_real, sf2so_dmx_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmx_real, nIrrep, SFS, MULT)
        sf2so_dmy_real, sf2so_dmy_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmy_real, nIrrep, SFS, MULT)
        sf2so_dmz_real, sf2so_dmz_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_dmz_real, nIrrep, SFS, MULT)

        sf2so_lx_real, sf2so_lx_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_lx_real, nIrrep, SFS, MULT)
        sf2so_ly_real, sf2so_ly_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_ly_real, nIrrep, SFS, MULT)
        sf2so_lz_real, sf2so_lz_imag = sf2so_prop_trafo(eigvecs_real, eigvecs_imag, sf_lz_real, nIrrep, SFS, MULT)

        sf2so_sx_real, sf2so_sx_imag = sf2so_spin_mat_trafo(eigvecs_real, eigvecs_imag, sf_sx_real, SFS, MULT)
        sf2so_sy_imag, sf2so_sy_real = sf2so_spin_mat_trafo(eigvecs_real, eigvecs_imag, sf_sy_imag, SFS, MULT)
        sf2so_sz_real, sf2so_sz_imag = sf2so_spin_mat_trafo(eigvecs_real, eigvecs_imag, sf_sz_real, SFS, MULT)

        # Write spin-orbit properties to output files
        with open('dipole-1.txt', 'w') as f1, \
             open('dipole-2.txt', 'w') as f2, \
             open('dipole-3.txt', 'w') as f3, \
             open('angmom-1.txt', 'w') as f4, \
             open('angmom-2.txt', 'w') as f5, \
             open('angmom-3.txt', 'w') as f6, \
             open('spin-1.txt', 'w') as f7, \
             open('spin-2.txt', 'w') as f8, \
             open('spin-3.txt', 'w') as f9:
            # Write header to output files
            f1.write(" #NROW NCOL REAL IMAG\n")
            f2.write(" #NROW NCOL REAL IMAG\n")
            f3.write(" #NROW NCOL REAL IMAG\n")
            f4.write(" #NROW NCOL REAL IMAG\n")
            f5.write(" #NROW NCOL REAL IMAG\n")
            f6.write(" #NROW NCOL REAL IMAG\n")
            f7.write(" #NROW NCOL REAL IMAG\n")
            f8.write(" #NROW NCOL REAL IMAG\n")
            f9.write(" #NROW NCOL REAL IMAG\n")
            # Write data to output files
            for iSOS in range(nSOS):
                for jSOS in range(nSOS):
                    f1.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmx_real[iSOS, jSOS]:16.8E}  {sf2so_dmx_imag[iSOS, jSOS]:16.8E}\n")
                    f2.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmy_real[iSOS, jSOS]:16.8E}  {sf2so_dmy_imag[iSOS, jSOS]:16.8E}\n")
                    f3.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_dmz_real[iSOS, jSOS]:16.8E}  {sf2so_dmz_imag[iSOS, jSOS]:16.8E}\n")
                    f4.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_lx_real[iSOS, jSOS]:16.8E}  {sf2so_lx_imag[iSOS, jSOS]:16.8E}\n")
                    f5.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_ly_real[iSOS, jSOS]:16.8E}  {sf2so_ly_imag[iSOS, jSOS]:16.8E}\n")
                    f6.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_lz_real[iSOS, jSOS]:16.8E}  {sf2so_lz_imag[iSOS, jSOS]:16.8E}\n")
                    f7.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_sx_real[iSOS, jSOS]:16.8E}  {sf2so_sx_imag[iSOS, jSOS]:16.8E}\n")
                    f8.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_sy_real[iSOS, jSOS]:16.8E}  {sf2so_sy_imag[iSOS, jSOS]:16.8E}\n")
                    f9.write(f"{iSOS + 1:4d}  {jSOS + 1:2d}  {sf2so_sz_real[iSOS, jSOS]:16.8E}  {sf2so_sz_imag[iSOS, jSOS]:16.8E}\n")

    # Handle exceptions
    except Exception as e:
        print(f"An error occurred: {str(e)}")


if __name__ == '__main__':
    main()


