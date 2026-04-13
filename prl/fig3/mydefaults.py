import numpy as np
from matplotlib import colors
import datetime
import os

today_ = datetime.datetime.today()
today  = today_.strftime("%Y_%m_%d")

plt_rcparams = {'text.usetex' : True,
                'font.size' : 8,
                'font.family' : 'serif',
                'text.latex.preamble' : r"\usepackage{lmodern} \usepackage{amstext}",
                'figure.figsize' : [3.4,3.4*0.7],
                'figure.dpi': 200
                }

EPS = 1e-10

def fint(num):
    return int(num+1e-10)
    
def mkdir(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)

def sci(n,dec=0):
    if n == 0:
        return "0"
    exponent = np.log10(np.abs(n))
    if exponent>=0:
        exponent = int(np.floor(exponent)+1e-10)
    else:
        exponent = int(np.floor(exponent)-1e-10)
    mantissa = n / 10**exponent
    if mantissa == 1:
        return f"$10^{{{exponent}}}$"
    else:
        return f"${mantissa:.{dec}f} \\times 10^{{{exponent}}}$"


_cs1 = '#FFFFFF #FFF5C6 #F9EE8D #CBEF96 #A1DFEF #8FD8DB #9BB7F9 #96A5FA #A48AF9 #B36EF8 #C062A5 #C55B4F #D54A0E #627207 #277B08 #0C645D #094F67 #0D3B5C #00296D #031B42 #08053E #000000'
_cs1 = np.array([colors.to_rgb(x) for x in _cs1.split(' ')])
cmWarm = colors.LinearSegmentedColormap.from_list('cmWarm', _cs1)

_cs2 = '#ffffff,#ffefe8,#fedfcf,#fbd2b0,#eeca92,#dfc276,#ccbc59,#aeb757,#8faf6d,#6fa589,#509ba4,#348fc2,#1886d0,#1372de,#145ce3,#2c3be3,#3f1fd6,#4f07bc,#450090,#320054,#160516,#000000'
_cs2 = np.array([colors.to_rgb(x) for x in _cs2.split(',')])
cmCool = colors.LinearSegmentedColormap.from_list('cmCool', _cs2)

_cs3 = ['#27186e', '#5d3992', '#875ec4', '#a487fe', '#c0bcff', '#daf6e9', '#80ef99', '#71c24c', '#a56741', '#a0273f', '#760113']
_cs3 = np.array([colors.to_rgb(x) for x in _cs3])
cmRainbow = colors.LinearSegmentedColormap.from_list('cmRainbow', _cs3)

_cs4 = ['#004409', '#226d52', '#4893a5', '#86b2d8', '#c5d8e9', '#ffffff', '#e9cbcb', '#f1889a', '#9e65b4', '#7b2eaf', '#3c099f']
_cs4 = np.array([colors.to_rgb(x) for x in _cs4])
cmBluePink = colors.LinearSegmentedColormap.from_list('cmBluePink', _cs4)