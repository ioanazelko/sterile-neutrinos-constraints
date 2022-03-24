### Created March 2022, I.Zelko
import os
import sys


import matplotlib.pyplot as plt
import numpy as np 


DARK_MATTER_DATA_LOCATION=os.environ["DARK_MATTER_DATA_LOCATION"]
DARK_MATTER_PAPER_LOCATION=os.environ["DARK_MATTER_PAPER_LOCATION"]
DARK_MATTER_PLOTS_LOCATION = DARK_MATTER_DATA_LOCATION+ "/plots"
DARK_MATTER_CODE_LOCATION = os.environ["DARK_MATTER_CODE_LOCATION"]

sys.path.insert(0, DARK_MATTER_CODE_LOCATION)


import distribution_functions as dist
import model_utils
import plot_utils


## Relation between m_WDM and m_hm:
## m_hm = 3E8 (m_WDM/3.3keV)**(-3.33) M_Sun
## so:
##m_WDM [keV] = 3.3*((m_hm [M_Sub]/3E8)**(-1/3.33) )
## converting from M_sun to keV
#M_sun = 1.989E30 #kg
#keV=1.782662E-36 #kg

#m_thWDM_d = 3.3*(m_hm/3E8)**(-1/3.33) #d stands for decreasing, because the array is decreasing in m_WDM right now
## p_thWDM = p_hm * dm_hm/d m_WDM
#p_thWDM_d = p_m_hm*3E8/(3.3/3.33*(m_hm/3E8)**(-1/3.33-1))
def const_thWDM_hm(case):
    if case=="no_baryons":
        const = 3.01639E8
    elif case == "include_baryons":
        const = 3.80877E8
    else:
        raise ValueError("Please specify correct case for the constant")
    return const
def m_hm_to_m_thWDM(m_hm,case):
    const = const_thWDM_hm(case)
    m_thWDM_d = 3.3*(m_hm/const)**(-1/3.33) #d stands for decreasing, because the array is decreasing in m_WDM right now
    return m_thWDM_d
def m_thWDM_to_m_hm(m_thWDM,case):
    const = const_thWDM_hm(case)
    m_hm_d = const*(m_thWDM/3.3)**(-3.33)  #d stands for decreasing because the array is decreasing right now 
    return m_hm_d
def deriv_m_hm_to_m_thWDM(m_hm,case):
    const = const_thWDM_hm(case)
    deriv=const*(-3.33)/3.3*(m_hm/const)**(4.33/3.33)
    return deriv
def calculate_thWDM(m_hm,p_m_hm,case):
    m_thWDM_d = m_hm_to_m_thWDM(m_hm, case=case)
    p_thWDM_d = p_m_hm*np.abs(deriv_m_hm_to_m_thWDM(m_hm,case))
    ## Order the arrays to be increasing in the m_WDM so that I don't mess with the integration calculations
    p_thWDM = np.flip(p_thWDM_d)
    m_thWDM = np.flip(m_thWDM_d)
    ## Normalize the arrays
    m_thWDM_lower_limit, m_thWDM_upper_limit, p_norm_thWDM = dist.norm_and_limit(m_thWDM,p_thWDM)
    return m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM

def plot_thWDM(m_hm,p_m_hm):
    fontsize=8
    plot_utils.set_plot_options(fontsize)
    fig, ax = plt.subplots(figsize=(3.6, 3.))
    color_dict =plot_utils.colorblind_color_dict_15()
    color_list = [color_dict["cb_blue"],color_dict["cb_light_blue"]]

    case_list = ["include_baryons","no_baryons" ]
    label_list= ["Case I","Case II"]

    for case_index in range(len(case_list)):
        case  = case_list[case_index]
        color = color_list[case_index]
        label=label_list[case_index]
        print("Thermal relic WDM, for case ", label," ", case)
        
        m_thWDM_lower_limit, m_thWDM_upper_limit, m_thWDM , p_thWDM = calculate_thWDM(m_hm,p_m_hm,case)
        ax.plot(m_thWDM, p_thWDM, color=color, label=label)
        ax.axvline(x=m_thWDM_lower_limit, color=color,linestyle='--')

    ax.set_xlabel(r"$m_{\textrm{thWDM}}$[keV]",fontsize=fontsize)
    ax.set_xscale("log")
    ax.set_ylabel(r"$m_{\textrm{thWDM}} $ Normalized Posterior Distribution",fontsize=fontsize)
    ax.legend(loc="upper left",handlelength =1,fontsize=fontsize)
    plt.tight_layout()
    paper_plot = True
    if paper_plot==True:
        fig.savefig(DARK_MATTER_PAPER_LOCATION+"/m_thWDM_p_norm.pdf",bbox_inches = 'tight',pad_inches = 0.01)

