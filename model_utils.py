### I. Zelko,  January 2021
### Functions to define the cosmological models used in the analysis

import numpy as np


def cosmo_param_Kev_2021():

    H_0 = 67.81 #km/s/Mpc WHY?
    Omega_Matter = 0.30704
    Omega_Lambda = 0.69296
    Omega_Baryon = 0.04868
    T_cmb = 2.72548 #K
    h = H_0/100 # where it was assumed H0 is in km s-1 Mpc-1; this is called the Hubble parameter

    return H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, T_cmb, h


def cosmo_param_Planck_2018():
    H_0 = 67.36000  #Planck 2018; arXiv:1807.06211 
    Omega_Matter  = 0.31530 #Planck 2018; arXiv:1807.06211 
    Omega_Lambda = 0.68470  #Planck 2018; arXiv:1807.06211 
    Omega_Baryon   = 0.04930 #Planck 2018; arXiv:1807.06211 
    Omega_Neutrinos = 0.0014197466735303636 #Planck 2018; arXiv:1807.06211 
    T_cmb  =2.72548 #Fixsen 2009 https://arxiv.org/abs/0911.1955 
    h = H_0/100 # where it was assumed H0 is in km s-1 Mpc-1; this is called the Hubble parameter

    return H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, Omega_Neutrinos, T_cmb, h


def alpha(Omega_X, h, m_X, model_type):
    """
    Specifies the alphas used in Bode et al. 2001 and Viel et al. 2001, which were
    also used in Schneider et al. 2012, Gilman et al. 2020, Abazajian and Kusenko 2019
    
    Note: for the Viel case, they say this function is valid for k< 5h Mpc-1
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

def T_WDM(Omega_X,h,m_X,k, model_type):
    """
    Function that implements A9, the transfer function for WDM, of Bode et al. 2001
    in: k: assumed in units of [h/Mpc]
        model_type: choses between the constants used in Bode et al. 2001, or Viel et al. 2005
    """
    
    nu=mu(model_type=model_type)
    alph = alpha(Omega_X,h,m_X,model_type=model_type)
    
    T_K = (1+ (alph*k)**(2*nu) )**(-5/nu)
    return T_K

def mu(model_type):
    if model_type == "bode":
        return 1.2
    elif model_type == "viel":
        return 1.12
    else:
        raise ValueError("Specify if the model is Bode or Viel")

def polynomial_fit(x,y,deg):

    poly = np.polyfit(x, y, deg=deg)
    return poly
    
def power_law(x,a,b):
    return a*(x**b)