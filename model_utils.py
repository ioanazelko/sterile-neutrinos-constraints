### I. Zelko,  January 2021
### Functions to define the cosmological models used in the analysis

import numpy as np
from scipy import interpolate



##### Thermal Relic Warm Dark Mattter Functions
def alpha(Omega_X, h, m_X, model_type):
    """
    Specifies the alphas used in Bode et al. 2001 and Viel et al. 2005, which were
    also used in Schneider et al. 2012, Gilman et al. 2020, Abazajian and Kusenko 2019
    
    Note: for the Viel case, they say this function is valid for k< 5h Mpc-1
    in: 
        Omega_X - omega dark matter
        h - Hubble parameter  (H/(100 *(km/s)/Mpc))
        m_X - dark matter particle mass [keV]
        model_type  - either Bode2001 or Viel2005
    out: alpha
    """
    
    if model_type == "bode":
        ### assumes m_X is in keV
        g_X=1.5 #By convention, g is taken to be 1.5 for the warm particle, based on the equivalent contribution of a light neutrino species (according to Bode et al. 2001)
        alph = 0.048*(Omega_X/0.4)**0.15*(h/0.65)**1.3*((1/m_X)**1.15)*(1.5/g_X)**0.29
    elif model_type == "viel":
        
        alph = 0.049*(1/m_X)**1.11 *(Omega_X/0.25)**0.11*(h/0.7)**1.22
        
    else:
        raise ValueError("Specify if the model is Bode or Viel")
    return alph

def get_mu(model_type):
    """
    Constant that comes from the fits performed in Bode2001 and Viel2005
    """
    if model_type == "bode":
        return 1.2
    elif model_type == "viel":
        return 1.12
    else:
        raise ValueError("Specify if the model is Bode or Viel")
def T_WDM(Omega_X,h,m_X,k, model_type):
    """
    Function that implements A9, the transfer function for WDM, of Bode et al. 2001
    in: k: assumed in units of [h/Mpc]
        model_type: choses between the constants used in Bode et al. 2001, or Viel et al. 2005
    out: Transfer function relative to the CDM transfer function, as a function of wavenumber k
    """
    ### Get the value for the nu constant, depending on the fit used.
    mu=get_mu(model_type=model_type)
    ### Get the value of alpha, which depends on the thermal relic mass and the cosmology, as well as the fit used
    alph = alpha(Omega_X,h,m_X,model_type=model_type)
    
    T_K = (1+ (alph*k)**(2*mu) )**(-5/mu)
    return T_K
    
def get_k_hm_for_transfer_function(wavenumber, transfer_f):
    """
    For a given transfer function sampled at the wavenumber array, find the wavenumber where transfer function
    becomes 0.5 (also called the 'half-mode')
    """
    k_vs_T = interpolate.interp1d(transfer_f,wavenumber)
    k_hm= k_vs_T(0.5)
    return k_hm

def m_WDM_from_k_hm(k_hm, Omega_X, h, model_type):
    """
    For the case of thermal relic warm dark matter, there is an analytical expression between the 'half-mode' 
    scale at which the thermal relic WDM transfer function decreases to 1/2.
    in:
        k_hm - the wavenumber at which the transfer function decreases to 1/2. Assumed in units of [h/Mpc]
        Omega_X -  the abundance of the dark matter particle (ususally = Omega_Matter-Omega_Baryon)
        h - the dimmensionless Hubble constant = H/100 km/s/Mpc
        model_type - chooses between Bode2001 and Viel2005
    out:
        m_WDM - the thWDM mass which gives the transfer function that becomes 1/2 at the given k_hm. 
        Assumed in units of [keV]
        
    """
    mu =get_mu(model_type)
    if model_type == "viel":
        m_WDM = k_hm**(1/1.11)*(0.049*(2**(mu/5)-1)**(-1/2/mu) *  (Omega_X/0.25)**0.11*(h/0.7)**1.22 )**(1/1.11)
    elif model_type == "bode":
        raise ValueError("implement this")
    else:
        raise ValueError("wrong model")
    return m_WDM

def get_mWDM_for_transfer_function(wavenumber, transfer_f,Omega_X, h,th_WDM_model_type):
    """
    For a given transfer function, finds the corresponding thWDM particle mass that would give a transfer function
    with the same half-mode k_hm

    in:
        wavenumber - the wavenumber at which the transfer function decreases to 1/2. Assumed in units of [h/Mpc]
        transfer_f - the transfer function to be matched
        Omega_X -  the abundance of the dark matter particle (ususally = Omega_Matter-Omega_Baryon)
        h - the dimmensionless Hubble constant = H/100 km/s/Mpc
        th_WDM_model_type - chooses between Bode2001 and Viel2005
    out:
        m_WDM - the thWDM mass which gives the transfer function that becomes 1/2 at the given k_hm. 
        Assumed in units of [keV]
        
    """
    #### Find the half-mode
    k_hm = get_k_hm_for_transfer_function(wavenumber, transfer_f)
    #### Find the wdm mass which gives the wdm transfer function with the same half-mode
    m_WDM = m_WDM_from_k_hm(k_hm, Omega_X, h,th_WDM_model_type)
    return m_WDMs
### Fit functions

def polynomial_fit(x,y,deg):
    """
    General function to perform a polynomial fit of degree deg
    """
    poly = np.polyfit(x, y, deg=deg)
    return poly
    
def power_law(x,a,b):
    """
    General function that returns a power law function
    """
    return a*(x**b)



#### Cosmological Parameters 
def cosmo_param_Planck_2018():
    """
    Returns cosmological parameters to the values from Plank 2018
    """
    H_0 = 67.36000  #Planck 2018; arXiv:1807.06211 
    Omega_Matter  = 0.31530 #Planck 2018; arXiv:1807.06211 
    Omega_Lambda = 0.68470  #Planck 2018; arXiv:1807.06211 
    Omega_Baryon   = 0.04930 #Planck 2018; arXiv:1807.06211 
    Omega_Neutrinos = 0.0014197466735303636 #Planck 2018; arXiv:1807.06211 
    T_cmb  =2.72548 #Fixsen 2009 https://arxiv.org/abs/0911.1955 
    h = H_0/100 # where it was assumed H0 is in km s-1 Mpc-1; this is called the Hubble parameter

    return H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, Omega_Neutrinos, T_cmb, h


def cosmo_param_Kev_2021():
    """
    Returns the cosmological parameters used by Kev in his (old) transfer function calculations from 2021
    """

    H_0 = 67.81 #km/s/Mpc WHY?
    Omega_Matter = 0.30704
    Omega_Lambda = 0.69296
    Omega_Baryon = 0.04868
    T_cmb = 2.72548 #K
    h = H_0/100 # where it was assumed H0 is in km s-1 Mpc-1; this is called the Hubble parameter

    return H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, T_cmb, h

