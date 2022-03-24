### I. Zelko,  October 2021

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np




def colorblind_color_list_15():
    """
    https://somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
    http://mkweb.bcgsc.ca/biovis2012/

    These values are RGB
    """
    cb_black = (0,0,0)
    cb_dark_green=(0,73/255,73/255)
    cb_blue_green=(0,146/255,146/255)
    cb_blue=(0,109/255,219/255)
    cb_medium_blue=(109/255,182/255,255/255)
    cb_light_blue=(182/255,219/255,255/255)
    cb_bright_pink=(255/255,109/255,182/255)
    cb_light_pink=(255/255,182/255,119/255)
    cb_magenta=(182/255,109/255,255/255)
    cb_purple=(73/255,0,146/255)
    cb_red=(146/255,0,0)
    cb_brown=(146/255,73/255,0)
    cb_orange=(219/255,209/255,0)
    cb_bright_green=(36/255,255/255,36/255)
    cb_yellow=(255/255,255/255,109/255)
    color_list = [cb_black,cb_dark_green,cb_blue_green,cb_blue,cb_medium_blue,
                  cb_light_blue,cb_bright_pink,cb_light_pink,cb_magenta,cb_purple,
                  cb_red,cb_brown,cb_orange,cb_bright_green,cb_yellow]
    return color_list
    
def colorblind_color_dict_15():
    color_list = colorblind_color_list_15()
    color_name_dict = dict(
        cb_black= color_list[0],
        cb_dark_green= color_list[1],
        cb_blue_green= color_list[2],
        cb_blue= color_list[3],
        cb_medium_blue= color_list[4],
        cb_light_blue= color_list[5],
        cb_bright_pink= color_list[6],
        cb_light_pink= color_list[7],
        cb_magenta= color_list[8],
        cb_purple= color_list[9],
        cb_red= color_list[10],
        cb_brown= color_list[11],
        cb_orange= color_list[12],
        cb_bright_green= color_list[13],
        cb_yellow= color_list[14])
    return color_name_dict





def colorblind_color_list():
    blue = '#377eb8'
    orange = '#ff7f00'
    green = '#4daf4a'
    pink = '#f781bf'
    brown = '#a65628'
    purple ='#984ea3'
    grey = '#999999'
    red = '#e41a1c'
    yellow = '#dede00'
    color_list = [brown,red,orange,yellow,green,blue,purple,pink,grey]
    color_dict = {"cb_blue":'#377eb8',"cb_orange":'#ff7f00',"cb_green" : '#4daf4a',"cb_pink" : '#f781bf',"cb_brown" : '#a65628',"cb_purple":'#984ea3',"cb_grey": '#999999',"cb_red" : '#e41a1c',"cb_yellow": '#dede00'}
    return color_list, color_dict

def set_plot_options(fontsize=10):
    #sns.set_style("white")
    #sns.set_style("ticks", {"xtick.major.size": 4, "ytick.major.size": 4,})
    mpl.rcParams['axes.linewidth'] = 1
    mpl.rc('text', usetex=True)
    mpl.rc('font', size=fontsize)
    mpl.rc('font', family='serif')