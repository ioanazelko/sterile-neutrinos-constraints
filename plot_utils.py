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




def plot_colordict():
    """
    Function to be used to visualize the choice of colors
    that work for all colorblind types
    """
    c_dict = colorblind_color_dict_15()
    fig, ax =plt.subplots(figsize=(12,8))
    linewidth=8
    ax.plot(np.array(range(10)),c=c_dict["cb_black"], label="black", linewidth=linewidth)
    ax.plot(np.array(range(10))+1,c=c_dict["cb_dark_green"], label="dark_green", linewidth=linewidth)
    ax.plot(np.array(range(10))+2,c=c_dict["cb_blue_green"], label="blue_green", linewidth=linewidth)
    ax.plot(np.array(range(10))+3,c=c_dict["cb_blue"], label="blue", linewidth=linewidth)
    ax.plot(np.array(range(10))+4,c=c_dict["cb_medium_blue"], label="medium_blue", linewidth=linewidth)
    ax.plot(np.array(range(10))+5,c=c_dict["cb_light_blue"], label="light_blue", linewidth=linewidth)
    ax.plot(np.array(range(10))+6,c=c_dict["cb_bright_pink"], label="bright_pink", linewidth=linewidth)
    ax.plot(np.array(range(10))+7,c=c_dict["cb_light_pink"], label="light_pink", linewidth=linewidth)
    ax.plot(np.array(range(10))+8,c=c_dict["cb_magenta"], label="magenta", linewidth=linewidth)
    ax.plot(np.array(range(10))+9,c=c_dict["cb_purple"], label="purple", linewidth=linewidth)
    ax.plot(np.array(range(10))+10,c=c_dict["cb_red"], label="red", linewidth=linewidth)
    ax.plot(np.array(range(10))+11,c=c_dict["cb_brown"], label="brown", linewidth=linewidth)
    ax.plot(np.array(range(10))+12,c=c_dict["cb_orange"], label="orange", linewidth=linewidth)
    ax.plot(np.array(range(10))+13,c=c_dict["cb_bright_green"], label="bright_green", linewidth=linewidth)
    ax.plot(np.array(range(10))+14,c=c_dict["cb_yellow"], label="yellow", linewidth=linewidth)
    ax.legend(loc="lower right", fontsize=15)

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
