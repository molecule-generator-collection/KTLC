from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import sys

## UV data should be 10^5 L mol^-1 cm^-1
##peak_hight_para and d_forpeakfinder should be adjusted, depending on each system


ac2os = 4.32*(10**-4)

peak_hight_para = 8
d_for_peakfinder = 50
background = 0.01

def peak_find(x, y):

    max_y = max(y)
    min_y = min(y)

    height_sh = min_y+(max_y-min_y)/peak_hight_para

    peak_Index, _ = find_peaks(y, height=height_sh, distance=d_for_peakfinder)
    
    peak_WL=np.array(x[peak_Index])
    peak_IN=np.array(y[peak_Index])
    
    return peak_WL, peak_IN

def func(x, *params):

    num_func = int(len(params)/3)

    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3*i,3*(i+1),1))
        amp = params[int(param_range[0])]
        ctr = params[int(param_range[1])]
        wid = params[int(param_range[2])]
#        y = y + amp * np.exp( -((x - ctr)/wid)**2)
        y = y + amp/(np.pi*wid*(1+np.power((x - ctr)/wid,2))) 
        y_list.append(y)

    y_sum = np.zeros_like(x)
    for i in y_list:
        y_sum = y_sum + i

    y_sum = y_sum + params[-1]

    return y_sum

def fit_plot(x, *params):
    num_func = int(len(params)/3)
    y_list = []
    for i in range(num_func):
        y = np.zeros_like(x)
        param_range = list(range(3*i,3*(i+1),1))
        amp = abs(params[int(param_range[0])])
        ctr = abs(params[int(param_range[1])])
        wid = abs(params[int(param_range[2])])
#        y = y + amp * np.exp( -((x - ctr)/wid)**2) + params[-1]
        y = y + amp/(np.pi*wid*(1+np.power((x - ctr)/wid,2)))
        y_list.append(y)

    return y_list

def nu_integral(x, y):

    integral = 0.0
    for i in range(len(x)-1):
        integral = integral + (y[i]+y[i+1])*abs(x[i+1]-x[i])/2

    return  integral

usage = 'Usage; %s UV data' % sys.argv[0]

try:
    infilename = sys.argv[1]
except:
    print (usage); sys.exit()

df_UV_DCM = pd.read_table(infilename,sep='\s+', names = ['WL','Int'])

x = df_UV_DCM['WL']
omega = (1/df_UV_DCM['WL'])*(10**7)
y = df_UV_DCM['Int']

peak_center, peak_hight = peak_find(x, y)

print(peak_center)

#[amp,ctr,wid]
guess = []

for i in range(len(peak_center)):
#	peak_paraset = [5.0, peak_center[i], peak_hight[i]]
	peak_paraset = [peak_hight[i], peak_center[i], 2.0]
	guess.append(peak_paraset)


guess_total = []
for i in guess:
    guess_total.extend(i)
guess_total.append(background)

popt, pcov = curve_fit(func, x, y, p0=guess_total,  maxfev=5000)

#print(popt)

fit = func(x, *popt)
plt.scatter(x, y, s=20)
plt.plot(x, fit , ls='-', c='black', lw=1)

y_list = fit_plot(x, *popt)
#baseline = np.zeros_like(x) + popt[-1]

os = []
opt_peak_info = []
for n,i in enumerate(y_list):
    os.append(ac2os*nu_integral(omega, y_list[n]))
    opt_peakp, opt_peakh = peak_find(x, y_list[n])
    if len(opt_peakp) != 0:
        opt_peak_nm = opt_peakp[0]
    else:
        opt_peak_nm = '-'

    opt_peak_info.append(opt_peak_nm)

print ('After peak fitting, peak positions and their oscillator strengths')
print (f'Parameters background = {background}, peak_hight_para={peak_hight_para}, d_for_peakfinder={d_for_peakfinder}')
print (opt_peak_info)
print (os)

for n,i in enumerate(y_list):
    plt.fill_between(x, i, 0, facecolor=cm.rainbow(n/len(y_list)), alpha=0.6)

plt.show()
