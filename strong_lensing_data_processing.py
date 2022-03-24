### Created March 2022, I.Zelko
import os
import sys

import matplotlib.pyplot as plt
import numpy as np 
import pickle

DARK_MATTER_CODE_LOCATION = os.environ["DARK_MATTER_CODE_LOCATION"]
DARK_MATTER_PAPER_LOCATION=os.environ["DARK_MATTER_PAPER_LOCATION"]


sys.path.insert(0, DARK_MATTER_CODE_LOCATION)

import distribution_functions as dist
import plot_utils


def load_posterior(paper):
    """
    Load the likelihood, which is provided as a 2D python function that has been pickled

    in: paper: Gilman2020 loads the posterior from Gilman et al. 2020 https://arxiv.org/abs/1908.06983
               Nadler2021 load the strong lensing posterior from Nadler et al. 2021  https//arxiv.org/pdf/2101.07810.pdf
    

    """
    objects = []
    if paper == "Nadler2021":
        filename = DARK_MATTER_CODE_LOCATION+"/strong_lensing_data/nadler2021_lensing_likelihood"
    elif paper == "Gilman2020":
        filename = DARK_MATTER_CODE_LOCATION+"/strong_lensing_data/gilman2020_likelihood_full"
    else:
        raise ValueError("specify the origina of the strong lensing data")
    with (open(filename, "rb")) as openfile:
        while True:
            try:
                objects.append(pickle.load(openfile))
            except EOFError:
                break
    f=objects[0] # returns the posterior function for a combination of m_hm  and \Sigma_sub  
    return f



def marginalizing_over_sigma(f, samples=10000):
    """
    Marginalizing over Sigma_sub

    in: 2D python function that gives the likelihood from all quad system conbined, as a function of log10(mhm) and Sigma_sub

    """

    
    marginalization_samples = samples
    sigma = np.linspace(0, 10.0, marginalization_samples)*1E-2 #\Sigma_sub [kpc-2]
    log10_m_hm = np.linspace(4.8,10.,marginalization_samples*2)  #log10(m_hm) [M_Sun]
    sigma_dif = (sigma[1]-sigma[0])*np.ones(marginalization_samples,dtype=np.float64)
    sigmav, log10_m_hmv = np.meshgrid(sigma,log10_m_hm)
    joint_flat = f(list(zip(log10_m_hmv.flatten(),sigmav.flatten())))
    joint_dist_array = joint_flat.reshape(marginalization_samples*2, marginalization_samples)
    #### Sigma upper level is 10, I should find out why, it is probably explained in the paper.
    #### Lower level is 0
    #### log10_mhm lower level is 5, upper level is 10
    fontsize=10
    plot_utils.set_plot_options(fontsize=fontsize)
    color_dict =plot_utils.colorblind_color_dict_15()
    fig, ax = plt.subplots(figsize=(10,10))
    x=sigma
    y=log10_m_hm
    color_map=ax.imshow(joint_dist_array,origin='lower',interpolation='none',extent=[x.min(),x.max(),y.min(),y.max()])
    color_map.set_cmap("Blues_r")
    ax.set_xlabel(r"$\Sigma_{\textrm{sub}}$[$\textrm{kpc}^{-2}$]")
    ax.set_ylabel(r"$\log_{10}(m_{\textrm{hm}})[M_{\odot}]$")
    ax.set_aspect(0.01)
    #     fig.colorbar(color_map,ax=ax, anchor=(0, 0.5), shrink=0.4)
    multiplied_dx = joint_dist_array*sigma_dif
    marginalized_dist = np.sum(multiplied_dx,axis=1)
    p_log10_mhm = dist.normalize_the_dist(log10_m_hm,marginalized_dist)
    return log10_m_hm,  p_log10_mhm



def plot_m_hm_posterior(log10_m_hm,p_log10_mhm):
    lower_95_limit, upper_95_limit, p_norm = dist.norm_and_limit(log10_m_hm, p_log10_mhm)
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    color_dict =plot_utils.colorblind_color_dict_15()

    fig2, ax2 = plt.subplots(figsize=(3.6, 3.))
    color=color_dict["cb_blue"]
    ax2.plot(log10_m_hm, p_norm,color=color)
    ax2.set_xlabel(r"$\log_{10}(m_{\textrm{hm}}[M_{\odot}])$")
    ax2.axvline(x=upper_95_limit, color=color,linestyle='--')

    #ax2.set_title(r"$m_{\textrm{hm}} $ Marginalized Posterior Distribution",fontsize=fontsize,pad=fontsize)
    ax2.set_ylabel(r"$\log_{10}(m_{\textrm{hm}}) $ Marginalized Posterior Dist. $\log_{10}(m_{\textrm{hm}}) $",fontsize=fontsize)
    fig2.savefig(DARK_MATTER_PAPER_LOCATION+"/m_hm_posterior_log.pdf",bbox_inches = 'tight',pad_inches = 0.01)

def mhm_log_to_linear(log10_m_hm,p_log10_mhm):
    """
    Moving from log10(mhm) posterior to linear posterior

    """

    m_hm = 10**log10_m_hm
    p_m_hm = p_log10_mhm/m_hm/np.log(10.)
    dist.norm_and_limit(m_hm, p_m_hm)
    fontsize=8
    plot_utils.set_plot_options(fontsize=fontsize)
    color_dict =plot_utils.colorblind_color_dict_15()

    fig, ax = plt.subplots(figsize=(3.6, 3.))

    ax.plot(m_hm, p_m_hm, color=color_dict["cb_blue"])
    ax.set_xlabel(r"$m_{\textrm{hm}}[M_{\odot}]$")
    ax.set_xscale("log")
    ax.set_title(r"$m_{\textrm{hm}} $ Marginalized Posterior Distribution",fontsize=fontsize,pad=fontsize*2)

    return m_hm, p_m_hm


def sample_the_posterior(f):
    """
    This is just for fast testing and visualisation set up
    """
    ### Define the  sigma and m_hm sampling space
    n_axis_samples = 30
    x = np.linspace(0.2, 9.8, n_axis_samples)*1E-2 #\Sigma_sub [kpc-2]
    y = np.linspace(5,10,n_axis_samples)  #log10(m_hm) [M_Sun]
    ### Make the meshgrid
    xv,yv = np.meshgrid(x,y)
    log_flat = f(list(zip(yv.flatten(),xv.flatten())))
    log_array = log_flat.reshape(n_axis_samples, n_axis_samples)
    
    ### Make the plot 
    fontsize=10
    plot_utils.set_plot_options(fontsize=fontsize)
    fig, ax = plt.subplots(figsize=(10,10))

    color_map=ax.imshow(log_array,origin='lower',interpolation='none',extent=[x.min(),x.max(),y.min(),y.max()])
    color_map.set_cmap("Blues_r")
    ax.set_xlabel(r"$\Sigma_{\textrm{sub}}$[$\textrm{kpc}^{-2}$]")
    ax.set_ylabel(r"$\log_{10}(m_{\textrm{hm}})[M_{\odot}]$")
    ax.set_aspect(0.01)
    fig.colorbar(color_map,ax=ax, anchor=(0, 0.5), shrink=0.4)
    #color_map.colorbar()
    fig2, ax2 = plt.subplots(figsize=(6,4))
    h = ax2.contourf(x, y, log_array)
    #plt.axis('scaled')
    fig2.colorbar(h)
    #plt.show()
    return log_array


# ### Check an interpolation method, it gives similar answers
# b = log10_m_hm
# cdf = np.cumsum(p_norm)
# cdf = cdf / float(cdf[-1])
# plt.plot(b, cdf)
# plt.xlim(4.8, 10)

# ### INVERT THE CDF, DETERMINE THE MASS M 95% OF SAMPLES ARE SMALLER ###
# cdf_inverse = interp1d(cdf, b)
# mhm_95 = cdf_inverse(0.954)
# print(mhm_95)