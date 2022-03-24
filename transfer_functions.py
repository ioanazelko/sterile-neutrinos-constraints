from __future__ import division,print_function

import os
import sys

from configparser import ConfigParser                  ### For parsing the configuration file
import csv
import h5py
import healpy as hp
import pandas as pd
import matplotlib as mpl             ### For plotting options
import matplotlib.pyplot as plt
import numpy as np
import pprint
from scipy.optimize import curve_fit
import time

DARK_MATTER_PAPER_LOCATION=os.environ["DARK_MATTER_PAPER_LOCATION"]
DARK_MATTER_CODE_LOCATION = os.environ["DARK_MATTER_CODE_LOCATION"]
sys.path.insert(0, DARK_MATTER_CODE_LOCATION)

import model_utils
import plot_utils




class TransferFunctions():
    def __init__(self, cosmo_constants_type, transfer_function_type,th_WDM_model_type):
        
        ### Choosing between Bode and Viel for the thermal WDM formula fit
        self.th_WDM_model_type=th_WDM_model_type
        ### Setting up the cosmological model parameters according to the choice of models
        if cosmo_constants_type == "Planck2018":
            self.H_0, self.Omega_baryon, self.Omega_Matter, self.Omega_Lambda, self.Omega_Neutrinos,self.T_cmb, self.h = model_utils.cosmo_param_Planck_2018()
        elif cosmo_constants_type == "Kev_2021":
            self.H_0, self.Omega_baryon, self.Omega_Matter, self.Omega_Lambda, self.T_cmb, self.h = model_utils.cosmo_param_Kev_2021()
        else:
            raise ValueError("The cosmological constants type is not valid")
        ################################################################################
        ### Loading up the transfer functions data according to the the choice of models
        ################################################################################
  
        if transfer_function_type == "PK":
            energy_list = ['1','2','3','7','14','15','30','60']
            file_list = ["PKTk"+ x + "keV.txt" for x in energy_list ]
            self.energy_str_list = [x + "keV" for x in energy_list ]
            self.energies = np.array(energy_list,dtype="float64")
            self.transfer_functions={}
            folder_location=DARK_MATTER_CODE_LOCATION+"/transfer_functions"

            for file_index in range(len(file_list)):
                file_name = file_list[file_index]
                energy = self.energies[file_index]
                loaded = np.loadtxt(folder_location+'/'+file_name)
                wavenumber= loaded.T[0]#[h/Mpc], the units the paper uses
                wavenumber_g = wavenumber*self.h #[1/Mpc], the units Galacticus uses
                transfer_f = loaded.T[1]
                self.transfer_functions[energy] = {}
                self.transfer_functions[energy]['wavenumber_g'] = wavenumber_g
                self.transfer_functions[energy]['wavenumber'] = wavenumber
                self.transfer_functions[energy]['transfer_f'] = transfer_f

        elif transfer_function_type == "KTY":
            energy_list = ['1','2','3','7','14']
            file_list = ["KTYTk"+ x + "keV.txt" for x in energy_list ]
            self.energy_str_list = [x + "keV" for x in energy_list ]
            self.energies = np.array(energy_list,dtype="float64")
            self.transfer_functions={}
            folder_location=DARK_MATTER_CODE_LOCATION+"/transfer_functions"

            for file_index in range(len(file_list)):
                file_name = file_list[file_index]
                energy = self.energies[file_index]
                loaded = np.loadtxt(folder_location+'/'+file_name)
                wavenumber= loaded.T[0]#[h/Mpc], the units the paper uses
                wavenumber_g = wavenumber*self.h #[1/Mpc], the units Galacticus uses
                transfer_f = loaded.T[1]
                self.transfer_functions[energy] = {}
                self.transfer_functions[energy]['wavenumber_g'] = wavenumber_g
                self.transfer_functions[energy]['wavenumber'] = wavenumber
                self.transfer_functions[energy]['transfer_f'] = transfer_f

        else:
            raise ValueError("The transfer function type is not valid")

    def plot_transfer_functions(self,galacticus_units=True,paper_plot=False):
        """
        Plotting the transfer function
        """
        ### Set the fontsize for the plot and the general properties
        plot_utils.set_plot_options(fontsize=17)
        ### Get the colorblind colors
        color_list =plot_utils.colorblind_color_list_15()
        full_fig, full_ax = plt.subplots(figsize=(7.4,3.8))

        for energy_index in range(len(self.energies)):
            energy = self.energies[energy_index]
            energy_label = self.energy_str_list[energy_index]
            color = color_list[energy_index]

            if galacticus_units == True:
                full_ax.plot(self.transfer_functions[energy]['wavenumber_g'], self.transfer_functions[energy]['transfer_f'],label=energy_label,color=color)
            else:
                ## to match with the units in the Abazajian et al. 2019 paper [h/Mpc]
                full_ax.plot(self.transfer_functions[energy]['wavenumber'] , self.transfer_functions[energy]['transfer_f'],label=energy_label,color=color)


        full_ax.set_xscale('log')

        full_ax.set_yscale('log')
        full_ax.legend(loc='lower left')
        if galacticus_units == True:
            full_ax.set_xlabel(r"k[$\textrm{Mpc}^{-1}$]")
        else:
            full_ax.set_xlabel(r"k[h $ \textrm{Mpc}^{-1}$]")
        full_ax.set_ylabel(r"$T_s(k)$")
        plt.tight_layout()
        full_ax.set_xlim([1E-1,1E3])
        full_ax.set_ylim([0.05,1])
        if paper_plot == True:
            
            full_fig.savefig(DARK_MATTER_PAPER_LOCATION+"/sn_transfer_functions.pdf",bbox_inches = 'tight',pad_inches = 0.01)
        #full_ax.set_xlim([None,1])
        #full_ax.set_ylim([0.99,1])
        #plt.tight_layout()
        #full_fig.savefig(DARK_MATTER_PAPER_LOCATION+"/transfer_functions_zoom.pdf",bbox_inches = 'tight',pad_inches = 0.01)
        #full_fig.clear()
        #plt.clf()
        #plt.close()
    
    
    def objective(self,k,m_X):
        Omega_X =self.Omega_Matter-self.Omega_baryon
        return model_utils.T_WDM(Omega_X,self.h,m_X,k,self.th_WDM_model_type)
    
    def fit_sn_with_WDM(self,paper_plot=False):
        """
        Fitting the sterile neutrinos transfer functions from 
        with the WDM transfer functions analytical formulas from Bode et al. (2001), eq. A9 
        """
        ### Set the fontsize for the plot and the general properties
        plot_utils.set_plot_options(fontsize=17)
        full_fig, full_ax = plt.subplots(figsize=(14, 7))
        ### Get the colorblind colors
        #color_list =plot_utils.colorblind_color_list()[0]
        color_list =plot_utils.colorblind_color_list_15()


        ### Index over energy and fit for the WDM mass
        sn_en= []
        wdm_en = []
        for energy_index in range(len(self.energies)):
            energy = self.energies[energy_index]
            energy_label = self.energy_str_list[energy_index]
            color = color_list[energy_index]
            wavenumber = self.transfer_functions[energy]['wavenumber']#[h/Mpc], the units the Bode et al. paper uses
            transfer_f = self.transfer_functions[energy]['transfer_f']
            full_ax.plot(wavenumber,transfer_f,label="sterile neutrino "+energy_label,c=color)

            # choose the input and output variables
            x, y = wavenumber, transfer_f
            # curve fit
            popt = curve_fit(self.objective, x, y)
            # summarize the parameter values
            m_WDM= popt[0][0]
            error = popt[1][0]
            T_wdm = self.objective(wavenumber, m_WDM)
            print('m_thWDM = %.5f keV +/- %.5e keV' % (m_WDM, error))
            #print("{:.2f}".format(m_WDM))
            print(m_WDM)
            sn_en.append(energy)
            wdm_en.append(m_WDM)
            full_ax.plot(wavenumber,T_wdm,label='thWDM '+"{:.2f}".format(m_WDM)+' keV',c=color, linestyle='--')


        full_ax.set_xscale('log')
        full_ax.set_yscale('log')
        full_ax.legend(loc='lower left')
        full_ax.set_ylabel(r"$T_s(k)$")
        full_ax.set_xlabel(r"k[h $ \textrm{Mpc}^{-1}$]")
        full_ax.set_xlim([1E-1,1E3])
        full_ax.set_ylim([0.05,1])    
        plt.tight_layout()

        if paper_plot == True:
            full_fig.savefig(DARK_MATTER_PAPER_LOCATION+"/transfer_functions_WDM_fit.pdf",bbox_inches = 'tight',pad_inches = 0.01)
            

        self.sn_en = sn_en
        self.wdm_en = wdm_en
        return np.array(sn_en), np.array(wdm_en)
    
    

    def fit_polynomials(self):
        degree_list = [1,2,3]

        poly_dict_thWDM_to_sn={}
        poly_dict_sn_to_thWDM = {}

        for degree_index in range(len(degree_list)):
            degree = degree_list[degree_index]
            poly_coeff_thWDM_to_sn =  model_utils.polynomial_fit(self.wdm_en,self.sn_en,deg=degree)
            poly_coeff_sn_to_thWDM =  model_utils.polynomial_fit(self.sn_en,self.wdm_en,deg=degree)
            poly_dict_thWDM_to_sn[degree] = poly_coeff_thWDM_to_sn
            poly_dict_sn_to_thWDM[degree] = poly_coeff_sn_to_thWDM 
        return poly_dict_thWDM_to_sn, poly_dict_sn_to_thWDM
    

    def fit_power_law(self):
        power_law_coeff, cov = curve_fit(f=model_utils.power_law, xdata=self.wdm_en, ydata=self.sn_en, p0=[0, 0], bounds=(-np.inf, np.inf))
        return power_law_coeff





def make_table(t1,t2):
    """
    Makes the table with the m_sn - m_thWDM fit coefficients, for the first 2 polynomial fits and for the power law fit
    """
    ### Load the coeficients from the transfer function objects
    
    poly_dict_thWDM_to_sn1, poly_dict_sn_to_thWDM1 =t1.fit_polynomials()
    poly_dict_thWDM_to_sn2, poly_dict_sn_to_thWDM1 =t2.fit_polynomials()
    power_law1 = t1.fit_power_law()
    power_law2 = t2.fit_power_law()
    count_list = ["&   1st &", "& 2nd &"]
    model_types = ["PK","KTY"]
    output = "\\begin{table}[t]\n"
    output+= r"\renewcommand{\arraystretch}{1.3}% for the vertical padding"+ "\n"
    output+= r"\setlength{\tabcolsep}{3pt}"+ "\n"
    output+="\\centering\n"
    output+="\\begin{tabular}{lllll}\n"
    output+="\\hline\n"
    output+="\\hline\n"
    output+=" & deg ~  & $a_0$~ & $a_1$~ & $a_2$~ \\\\ \n"
    output+="\\hline\n"
    for i in range(2):
        poly_dict_thWDM_to_sn = [poly_dict_thWDM_to_sn1,poly_dict_thWDM_to_sn2][i]
        model = model_types[i]
        output+="\\multirow{2}{*}{"+model+"} "
        for degree in range(1,3):
            coeffs= np.flip(poly_dict_thWDM_to_sn[degree])
            count = count_list[degree-1]
            output+= count+'&'.join(['{:.2e}'.format(coeff) for coeff in coeffs]) +"\\\\ \n"
        output+="\\hline\n"
    output += "power & law ~  & $a$~ & $b$~ \\\\ \n"
    output+="\\hline\n"
    
    for i in range(2):
        power_coeffs  = [power_law1, power_law2][i]
        model = model_types[i]
        output += model + "& &" + '&'.join(['{:.2e}'.format(coeff) for coeff in power_coeffs]) +"\\\\ \n"
        output+="\\hline\n"

    output += "\end{tabular}\n"
    output +="\caption{Coefficients for the fits for the relation between the $\\msn$ and $\\mth$, for the cases of the Higgs production mechanism (PK) and GUT scale (KTY). The power law fit, as well as the first 2 degree polynomials are shown. \label{table:fit_coeffs}}\n\end{table}"
    
    filename = DARK_MATTER_PAPER_LOCATION+"/fit_coeffs.tex"
    f = open(filename,'w')
    f.write(output)
    f.close()


def plot_sn_vs_wdm_overlay(t1, t2,paper_plot=False):
    
    """
    Plot the first 2 degree polynomials and the power law
    """
    ### Set the fontsize for the plot and the general properties
    plot_utils.set_plot_options(fontsize=10)
    ### Get the colorblind colors
    color_dict =plot_utils.colorblind_color_dict_15()
    
    ### Add the  polynomial fits
    degree_list = [1,2]
    fontsize=8
    color_list = [color_dict["cb_light_blue"],color_dict["cb_bright_pink"],color_dict["cb_orange"] ]
    fig, ax= plt.subplots(figsize=(3.6, 2.44))
    ax.scatter(t1.wdm_en, t1.sn_en,color = color_dict["cb_blue"],label="PK")
    ax.scatter(t2.wdm_en, t2.sn_en,color = color_dict["cb_red"],label="KTY")
    for degree_index in range(len(degree_list)):
        degree = degree_list[degree_index]
        color = color_list[degree_index]
        poly_coeff1 =  model_utils.polynomial_fit(t1.wdm_en,t1.sn_en,deg=degree)
        poly_coeff2 =  model_utils.polynomial_fit(t2.wdm_en,t2.sn_en,deg=degree)



        #ax.plot( wdm_en, sn_en, color = color_dict["cb_blue"],linestyle="None")
        ax.plot(t1.wdm_en,np.polyval(poly_coeff1, t1.wdm_en), label="Poly. deg "+str(degree),color=color )
        ax.plot(t2.wdm_en,np.polyval(poly_coeff2, t2.wdm_en),color=color )
    #### Add the power law fits
    power_law1 = t1.fit_power_law()
    power_law2 = t2.fit_power_law()
    color = color_list[2]
    
    ax.plot(t1.wdm_en,model_utils.power_law(t1.wdm_en, power_law1[0], power_law1[1]), label="Power Law",color=color )
    ax.plot(t2.wdm_en,model_utils.power_law(t2.wdm_en, power_law2[0], power_law2[1]),color=color)
    
    #a,b=find_relation_between_sn_and_wdm(wdm_en,sn_en)
    #ax.plot(wdm_en,a*wdm_en+b, label="linear fit", color =color_dict["cb_bright_pink"])
    ax.set_xlabel("Thermal Relic WDM Mass [keV]")
    ax.set_ylabel("Sterile Neutrino Mass [keV]")
    ax.legend(loc='upper left',fontsize=fontsize)
    plt.tight_layout()
    if paper_plot==True:
        fig.savefig(DARK_MATTER_PAPER_LOCATION+"/m_sn_vs_m_thWDM.pdf",bbox_inches = 'tight',pad_inches = 0.01)

def plot_sn_vs_wdm_polynomials(t1, t2,paper_plot=False):
    
    """
    Plots the polynomial fits up to 3rd degree
    """
    ### Set the fontsize for the plot and the general properties
    plot_utils.set_plot_options(fontsize=10)
    ### Get the colorblind colors
    color_dict =plot_utils.colorblind_color_dict_15()
    
    ### Add the fit
    degree_list = [1,2,3]
    fontsize=8
    color_list = [color_dict["cb_light_blue"],color_dict["cb_bright_pink"],color_dict["cb_orange"] ]
    fig, ax= plt.subplots(figsize=(3.6, 2.44))
    ax.scatter(t1.wdm_en, t1.sn_en,color = color_dict["cb_blue"],label="PK")
    ax.scatter(t2.wdm_en, t2.sn_en,color = color_dict["cb_red"],label="KTY")
    for degree_index in range(len(degree_list)):
        degree = degree_list[degree_index]
        color = color_list[degree_index]
        poly_coeff1 =  model_utils.polynomial_fit(t1.wdm_en,t1.sn_en,deg=degree)
        poly_coeff2 =  model_utils.polynomial_fit(t2.wdm_en,t2.sn_en,deg=degree)



        #ax.plot( wdm_en, sn_en, color = color_dict["cb_blue"],linestyle="None")
        ax.set_xlabel("Thermal Relic WDM Mass [keV]")
        ax.set_ylabel("Sterile Neutrino Mass [keV]")
        ax.plot(t1.wdm_en,np.polyval(poly_coeff1, t1.wdm_en), label="Poly. deg "+str(degree),color=color )
        ax.plot(t2.wdm_en,np.polyval(poly_coeff2, t2.wdm_en),color=color )

    #a,b=find_relation_between_sn_and_wdm(wdm_en,sn_en)
    #ax.plot(wdm_en,a*wdm_en+b, label="linear fit", color =color_dict["cb_bright_pink"])
    ax.legend(loc='upper left',fontsize=fontsize)
    plt.tight_layout()
    


def make_table_3rd_degree(t1,t2):
    """
    Makes the table with the m_sn - m_thWDM fit coefficients, up to the 3rd degree
    """
    ### Load the coeficients from the transfer function objects
    
    poly_dict_thWDM_to_sn1, poly_dict_sn_to_thWDM1 =t1.fit_polynomials()
    poly_dict_thWDM_to_sn2, poly_dict_sn_to_thWDM1 =t2.fit_polynomials()
    count_list = ["&   1st &", "& 2nd &","& 3rd &" ]
    model_types = ["PK","KTY"]
    output = "\\begin{table}[t]\n"
    output += r"\setlength{\tabcolsep}{3pt}"+ "\n"
    output+="\\centering\n"
    output+="\\begin{tabular}{llllll}\n"
    output+="\\hline\n"
    output+="\\hline\n"
    output+=" & deg ~  & $a_0$~ & $a_1$~ & $a_2$~ & $a_3$~ \\\\ \n"
    output+="\\hline\n"
    for i in range(2):
        poly_dict_thWDM_to_sn = [poly_dict_thWDM_to_sn1,poly_dict_thWDM_to_sn2][i]
        model = model_types[i]
        output+="\\multirow{3}{*}{"+model+"} "
        for degree in range(1,4):
            coeffs= np.flip(poly_dict_thWDM_to_sn[degree])
            count = count_list[degree-1]
            output+= count+'&'.join(['{:.2e}'.format(coeff) for coeff in coeffs]) +"\\\\ \n"
        output+="\\hline\n"
    output += "\end{tabular}\n"
    output +="\caption{Coefficients for the polynomials fits for the relation between the $\\msn$ and $\\mth$, for the cases of the kev miracle model and GUT scale. \label{table:fit_coeffs}}\n\end{table}"
    
    filename = DARK_MATTER_PAPER_LOCATION+"/fit_coeffs_3rd_degree.tex"
    f = open(filename,'w')
    f.write(output)
    f.close()