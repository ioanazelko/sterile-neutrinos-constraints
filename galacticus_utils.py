### I. Zelko,  December 2021

import matplotlib.pyplot as plt
import numpy as np



def get_transfer_function(file, output):
    """
    Returns the wavenumber k and the transfer function for the given file and output group
    """
    transfer_function = file['Outputs'][output]['transferFunction'][()]
    wavenumber = file['Outputs'][output]['wavenumber'][()]
    return wavenumber, transfer_function

def read_halo_mass_function(file, output):
    
    """
    input: 
        file: hdf5 object
        output: output group from the hdf5 object
    output:
        halo_mass: field halo mass in units of Solar Masses
        halo_mass_f: for the given output, return specifically the number of halos per unit mass, per natural logarithm of halo mass; in units of Mpc-3
    """
    f = h5py.File(file,'r')
    
    halo_mass=f['Outputs'][output]['haloMass'][()] #(in units of Solar masses
    halo_mass_f=f['Outputs'][output]['haloMassFunctionLnM'][()] #specifically the number of halos per unit mass, per natural logarithm of halo mass; in units of Mpc-3
    f.close()
    return halo_mass, halo_mass_f