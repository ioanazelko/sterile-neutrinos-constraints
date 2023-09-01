### Created March 2022, I.Zelko
import os
import sys


import matplotlib.pyplot as plt
import numpy as np 


DARK_MATTER_PAPER_LOCATION=os.environ["DARK_MATTER_PAPER_LOCATION"]
DARK_MATTER_CODE_LOCATION = os.environ["DARK_MATTER_CODE_LOCATION"]

sys.path.insert(0, DARK_MATTER_CODE_LOCATION)


import distribution_functions as dist
import model_utils
import thermal_relic as th
import plot_utils



def m_sn_from_m_thWDM(poly_dict_thWDM_to_sn, m_thWDM, degree):
    mass_coefficients = poly_dict_thWDM_to_sn[degree]
    m_sn = np.zeros(len(m_thWDM))
    for i in range(degree+1):
        m_sn+=mass_coefficients[i]*(m_thWDM**(degree-i))
    return m_sn
def p_sn_from_p_tWDM(poly_dict_thWDM_to_sn, m_thWDM,p_thWDM, degree):
    
    mass_coefficients = poly_dict_thWDM_to_sn[degree]
    d_m_sn_d_m_thWDM = np.zeros(len(p_thWDM))
    for i in range(degree+1):
        d_m_sn_d_m_thWDM += mass_coefficients[i]*(degree-i)*(m_thWDM**(degree-1-i))
    p_sn = p_thWDM /np.abs(d_m_sn_d_m_thWDM)
    return p_sn
def sn_from_thWDM(poly_dict_thWDM_to_sn, m_thWDM,p_thWDM, degree):
    print("Calculate sn for degree ",degree)
    m_sn = m_sn_from_m_thWDM(poly_dict_thWDM_to_sn, m_thWDM, degree)
    p_sn = p_sn_from_p_tWDM(poly_dict_thWDM_to_sn, m_thWDM,p_thWDM, degree)
    m_sn_lower_limit, m_sn_upper_limit, p_sn_norm = dist.norm_and_limit(m_sn,p_sn)
    return m_sn, p_sn_norm, m_sn_lower_limit
def plot_degrees(poly_dict_thWDM_to_sn,m_thWDM,p_thWDM):
    fig, ax  =  plt.subplots(figsize=(3.6,3))
    m_sn, p_sn_n, m_sn_lower_limit1 = sn_from_thWDM(poly_dict_thWDM_to_sn,m_thWDM,p_thWDM,  degree=1)
    print(p_sn_n)
    print(m_sn)
    ax.plot(m_thWDM,m_sn, label='1')
    m_sn, p_sn_n, m_sn_lower_limit2 = sn_from_thWDM(poly_dict_thWDM_to_sn,m_thWDM,p_thWDM,  degree=2)
    print(p_sn_n)
    print(m_sn)
    ax.plot(m_thWDM,m_sn,label='2')

    m_sn, p_sn_n, m_sn_lower_limit3 = sn_from_thWDM(poly_dict_thWDM_to_sn,m_thWDM,p_thWDM,  degree=3)
    print(p_sn_n)
    print(m_sn)
    ax.plot(m_thWDM,m_sn,label=3)
    ax.legend()
    ax.set_xlabel("Thermal WDM Mass [keV]")
    ax.set_ylabel("Sterile Neutrino Mass [keV]")
def run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn, degree=2):
    case_list = ["include_baryons","no_baryons" ]
    label_list= ["Case I","Case II"]
    m_sn_list = []
    p_sn_list = []
    sn_lower_limit_list = []
    for case_index in range(len(case_list)):
        case  = case_list[case_index]

        label=label_list[case_index]
        print("Calculating ",label," ", case)

        print("Calculating the thWDM posterior")
        m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM = th.calculate_thWDM(m_hm,p_m_hm,case)
        print("Calculating the sn posterior")
        m_sn, p_sn_n, m_sn_lower_limit = sn_from_thWDM(poly_dict_thWDM_to_sn,m_thWDM,p_thWDM,  degree=degree)
        m_sn_list.append(m_sn)
        p_sn_list.append(p_sn_n)
        sn_lower_limit_list.append(m_sn_lower_limit)
    return m_sn_list, p_sn_list,sn_lower_limit_list

def plot_sn_posteriors(m_hm,p_m_hm,t1,t2):
    """
    Makes the plots of the posteriors of the PL, KTY
    in: m_hm: half-mode mass array
        p_m_hm: posterior of the half-mode mass array
        t1, t2 - objects of the transfer function class, for PK and KTY
    """
  
    
    poly_dict_thWDM_to_sn1, poly_dict_sn_to_thWDM1 =t1.fit_polynomials()
    poly_dict_thWDM_to_sn2, poly_dict_sn_to_thWDM1 =t2.fit_polynomials()

    set_labels = ["PK","KTY"]
    m_sn_list1, p_sn_list1, sn_lower_limit_list1 = run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn1)
    m_sn_list2, p_sn_list2, sn_lower_limit_list2= run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn2)
    
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    color_dict =plot_utils.colorblind_color_dict_15()
    fig, ax = plt.subplots(figsize=(3.6, 3.))


    ax.plot(m_sn_list2[0], p_sn_list2[0],color=color_dict["cb_magenta"],label="KTY - I")
    ax.axvline(x=sn_lower_limit_list2[0], color=color_dict["cb_magenta"],linestyle='--')
    ax.plot(m_sn_list2[1], p_sn_list2[1],color=color_dict["cb_light_pink"],label="KTY - II")
    ax.axvline(x=sn_lower_limit_list2[1], color=color_dict["cb_light_pink"],linestyle='--')
    ax.plot(m_sn_list1[0], p_sn_list1[0],color=color_dict["cb_blue"], label="PK - I")
    ax.axvline(x=sn_lower_limit_list1[0], color=color_dict["cb_blue"],linestyle='--')
    ax.plot(m_sn_list1[1], p_sn_list1[1],color=color_dict["cb_light_blue"],label="PK - II")
    ax.axvline(x=sn_lower_limit_list1[1], color=color_dict["cb_light_blue"],linestyle='--')
    ax.set_xlabel(r"$m_{\textrm{sn}}$[keV]")
    ax.set_title(r"$m_{\textrm{sn}}$ Normalized Posterior",fontsize=fontsize)
    ax.set_xscale("log")
    ax.legend(loc="upper left",handlelength =1,fontsize=fontsize)
    plt.tight_layout()

def m_thWDM_to_msnDW(m_thWDM, Omega_DM):
    """
    Connection between thermal relic WDM and Dodelson-Widrow sterile neutrino 
    Eq. 5 from Viel et al 2005
    """
    msn = 4.43 *m_thWDM**(4/3) *(0.25/Omega_DM)**(1/3)
    return msn
def p_thWDM_to_p_snDW(m_thWDM, p_thWDM, Omega_DM):
    m_sn =  m_thWDM_to_msnDW(m_thWDM, Omega_DM)
    deriv_msn_to_m_thWDM = 4.43*(4/3)*(m_thWDM)**(1/3) *(0.25/Omega_DM)**(1/3)
    p_snDW = p_thWDM/np.abs(deriv_msn_to_m_thWDM)
    lower_95_limit, upper_95_limit, p_norm = dist.norm_and_limit(m_sn, p_snDW)
    return   lower_95_limit, upper_95_limit, m_sn, p_norm

def calculate_DW(m_hm,p_m_hm):
    """
    Plots the posterior of the Dodelson-Widrow model given a posterior of the half mode mass
    """
    H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, Omega_Neutrinos, T_cmb, h = model_utils.cosmo_param_Planck_2018()
    Omega_DM = Omega_Matter-Omega_Baryon-Omega_Neutrinos
    color_dict =plot_utils.colorblind_color_dict_15()

    
    case_list = ["include_baryons","no_baryons" ]
    label_list= ["Case I","Case II"]
    color_list = [color_dict["cb_dark_green"], color_dict["cb_blue_green"]]
    
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    fig, ax = plt.subplots(figsize=(3.6, 3.))
    
    for i in range(2):
        case = case_list[i]
        label=label_list[i]
        color = color_list[i]
        
        m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM = th.calculate_thWDM(m_hm,p_m_hm,case)
        lower_95_limit, upper_95_limit, m_sn, p_norm  = p_thWDM_to_p_snDW(m_thWDM, p_thWDM, Omega_DM)
        
   
        ax.plot(m_sn, p_norm ,label=label, color=color)
        ax.axvline(x=lower_95_limit,linestyle='--', color=color)
        
        
        
    
    ax.set_xlabel(r"$m_{\textrm{sn}}$[keV]")
    ax.set_title(r"$m_{\textrm{sn}}$ Normalized Posterior",fontsize=fontsize)
    ax.set_xscale("log")
    ax.legend(loc="upper left",handlelength =1)
    plt.tight_layout()
    return 


def plot_all_sn(m_hm,p_m_hm,t1,t2, paper_plot=False):

    """
    Makes the plots of the posteriors of the PL, KTY, and DW 
    """
     
    H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, Omega_Neutrinos, T_cmb, h = model_utils.cosmo_param_Planck_2018()
    Omega_DM = Omega_Matter-Omega_Baryon-Omega_Neutrinos
    color_dict =plot_utils.colorblind_color_dict_15()

    
    case_list = ["include_baryons","no_baryons" ]
    label_list= ["DW-I","DW-II"]
    color_list = [color_dict["cb_dark_green"], color_dict["cb_blue_green"]]
    
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    
    
    poly_dict_thWDM_to_sn1, poly_dict_sn_to_thWDM1 =t1.fit_polynomials()
    poly_dict_thWDM_to_sn2, poly_dict_sn_to_thWDM1 =t2.fit_polynomials()

    set_labels = ["PK","KTY"]
    m_sn_list1, p_sn_list1, sn_lower_limit_list1 = run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn1)
    m_sn_list2, p_sn_list2, sn_lower_limit_list2= run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn2)
    
    ### Plot model KTY and PK

    ax.plot(m_sn_list2[0], p_sn_list2[0],color=color_dict["cb_magenta"],label="KTY-I")
    ax.axvline(x=sn_lower_limit_list2[0], color=color_dict["cb_magenta"],linestyle='--')
    ax.plot(m_sn_list2[1], p_sn_list2[1],color=color_dict["cb_light_pink"],label="KTY-II")
    ax.axvline(x=sn_lower_limit_list2[1], color=color_dict["cb_light_pink"],linestyle='--')
    ax.plot(m_sn_list1[0], p_sn_list1[0],color=color_dict["cb_blue"], label="PK-I")
    ax.axvline(x=sn_lower_limit_list1[0], color=color_dict["cb_blue"],linestyle='--')
    ax.plot(m_sn_list1[1], p_sn_list1[1],color=color_dict["cb_light_blue"],label="PK-II")
    ax.axvline(x=sn_lower_limit_list1[1], color=color_dict["cb_light_blue"],linestyle='--')

   ### Plot the DW model
    for i in range(2):
        case = case_list[i]
        label=label_list[i]
        color = color_list[i]
        
        m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM = th.calculate_thWDM(m_hm,p_m_hm,case)
        lower_95_limit, upper_95_limit, m_sn, p_norm  = p_thWDM_to_p_snDW(m_thWDM, p_thWDM, Omega_DM)
        
   
        ax.plot(m_sn, p_norm ,label=label, color=color)
        ax.axvline(x=lower_95_limit,linestyle='--', color=color)
    ax.set_xlabel(r"$m_{\textrm{sn}}$[keV]")
    ax.set_ylabel(r"$m_{\textrm{sn}}$ Normalized Posterior Distribution",fontsize=fontsize)
    ax.set_xscale("log")
    ax.set_xlim([1E-2,1E3])
    ax.legend(loc="upper left",handlelength =1,fontsize=fontsize)
    plt.tight_layout()
    if paper_plot == True:
        fig.savefig(DARK_MATTER_PAPER_LOCATION+"/m_sn_p_norm_log_single.pdf",bbox_inches = 'tight',pad_inches = 0.01)
    return


def plot_all_sn_panels(m_hm,p_m_hm,t1,t2, paper_plot=False):

    """
    Makes the plots of the posteriors of the PL, KTY, and DW for the paper, in 4 separate panels
    """
     
    H_0, Omega_Baryon, Omega_Matter, Omega_Lambda, Omega_Neutrinos, T_cmb, h = model_utils.cosmo_param_Planck_2018()
    Omega_DM = Omega_Matter-Omega_Baryon-Omega_Neutrinos
    color_dict =plot_utils.colorblind_color_dict_15()

    
    case_list = ["include_baryons","no_baryons" ]
    label_list= ["DW-I","DW-II"]
    color_list = [color_dict["cb_dark_green"], color_dict["cb_blue_green"]]
    
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    
    fig, ax = plt.subplots( 2,2,figsize=(7.2, 6.))
  
    poly_dict_thWDM_to_sn1, poly_dict_sn_to_thWDM1 =t1.fit_polynomials()
    poly_dict_thWDM_to_sn2, poly_dict_sn_to_thWDM1 =t2.fit_polynomials()

    set_labels = ["PK","KTY"]
    print("----------")
    print("PK")
    m_sn_list1, p_sn_list1, sn_lower_limit_list1 = run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn1)
    print("----------")
    print("KTY")
    m_sn_list2, p_sn_list2, sn_lower_limit_list2= run_sn_all_cases(m_hm,p_m_hm,poly_dict_thWDM_to_sn2)
    
    
    ax[0,0].plot(m_sn_list1[0], p_sn_list1[0],color=color_dict["cb_blue"], label="PK-I")
    ax[0,0].axvline(x=sn_lower_limit_list1[0], color=color_dict["cb_blue"],linestyle='--')
    ax[0,0].plot(m_sn_list1[1], p_sn_list1[1],color=color_dict["cb_light_blue"],label="PK-II")
    ax[0,0].axvline(x=sn_lower_limit_list1[1], color=color_dict["cb_light_blue"],linestyle='--')
    ax[0,0].plot(m_sn_list2[0], p_sn_list2[0],color=color_dict["cb_magenta"],label="KTY-I")
    ax[0,0].axvline(x=sn_lower_limit_list2[0], color=color_dict["cb_magenta"],linestyle='--')
    ax[0,0].plot(m_sn_list2[1], p_sn_list2[1],color=color_dict["cb_light_pink"],label="KTY-II")
    ax[0,0].axvline(x=sn_lower_limit_list2[1], color=color_dict["cb_light_pink"],linestyle='--')


    ax[0,1].plot(m_sn_list1[0], p_sn_list1[0],color=color_dict["cb_blue"], label="PK-I")
    ax[0,1].axvline(x=sn_lower_limit_list1[0], color=color_dict["cb_blue"],linestyle='--')
    ax[0,1].plot(m_sn_list1[1], p_sn_list1[1],color=color_dict["cb_light_blue"],label="PK-II")
    ax[0,1].axvline(x=sn_lower_limit_list1[1], color=color_dict["cb_light_blue"],linestyle='--')
    ax[1,0].plot(m_sn_list2[0], p_sn_list2[0],color=color_dict["cb_magenta"],label="KTY-I")
    ax[1,0].axvline(x=sn_lower_limit_list2[0], color=color_dict["cb_magenta"],linestyle='--')
    ax[1,0].plot(m_sn_list2[1], p_sn_list2[1],color=color_dict["cb_light_pink"],label="KTY-II")
    ax[1,0].axvline(x=sn_lower_limit_list2[1], color=color_dict["cb_light_pink"],linestyle='--')

    print("----------")
    print("DW")
  
    for i in range(2):
        case = case_list[i]
        label=label_list[i]
        print("Calculating ",label," ", case)

        color = color_list[i]
        
        m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM = th.calculate_thWDM(m_hm,p_m_hm,case)
        lower_95_limit, upper_95_limit, m_sn, p_norm  = p_thWDM_to_p_snDW(m_thWDM, p_thWDM, Omega_DM)
        
        ax[0,0].plot(m_sn, p_norm ,label=label, color=color)
        ax[0,0].axvline(x=lower_95_limit,linestyle='--', color=color)
        ax[1,1].plot(m_sn, p_norm ,label=label, color=color)
        ax[1,1].axvline(x=lower_95_limit,linestyle='--', color=color)
    
    
    ax[0,0].set_xlabel(r"$m_{\textrm{sn}}$[keV]")
    ax[0,0].set_ylabel(r"$m_{\textrm{sn}}$ Normalized Posterior Distribution",fontsize=fontsize)
    ax[0,0].set_xscale("log")
    
    ax[0,0].set_xlim([1E-1,1E3])
    ax[0,0].legend(loc="upper left",handlelength =1,fontsize=fontsize)
    
    for (i,j) in [(0,1),(1,0),(1,1)]:
        ax[i,j].set_xlabel(r"$m_{\textrm{sn}}$[keV]")
        ax[i,j].set_ylabel(r"$m_{\textrm{sn}}$ Normalized Posterior Distribution",fontsize=fontsize)
        ax[i,j].legend(loc="upper right",handlelength =1,fontsize=fontsize)

    plt.tight_layout()
    if paper_plot == True:
        fig.savefig(DARK_MATTER_PAPER_LOCATION+"/m_sn_p_norm_log.pdf",bbox_inches = 'tight',pad_inches = 0.01)
    return