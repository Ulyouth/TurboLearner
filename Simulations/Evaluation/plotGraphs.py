#!/usr/bin/env python
"""
NAME:   plotGraphs.py
AUTHOR: Ulyouth
DATE:   05.06.2020 
DESC:   A script to plot the results of scripts "genSamples.sh", 
        "optzSamples.py" and "evalSamples.py".
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plotCoeffErrGraph(case, x_var, other_vars, err_data, save_path):
    plt.figure(figsize=(6, 6))
    plt.plot(err_data[x_var], err_data['Error'])
    plt.axis([min(err_data[x_var]), max(err_data[x_var]), 
              0, max(err_data['Error']) * 1.1])
    
    title = x_var + " x Error [Case " + case + "]\n" + x_var + "="\
        + str(min(err_data[x_var])) + "~" + str(max(err_data[x_var])) + "  "
    
    for x in range(0, len(other_vars), 2):  
        if other_vars[x] == x_var: continue
        title += (other_vars[x] + "=" + str(other_vars[x+1]) + "  ")
    
    plt.title(title, fontsize=14)
    plt.xlabel(x_var, fontsize=14)
    plt.ylabel("Error [%]", fontsize=14)
    plt.grid(True)
    plt.savefig(save_path + "//error.png")
    
def plotCoeffErrGraphEx(data, x_var, f_var, f, x0, xf, y0, yf):
    colors = ['r', 'g', 'b', 'r--', 'g--', 'b--', 
              'm', 'y', 'c', 'm--', 'y--', 'c--']
    labels = []
    title = x_var + " x Error [" + f_var + "=" + str(f) + "]"
    i = 0
    
    plt.figure(figsize=(6, 6))
    plt.axis([x0, xf, y0, yf])
    
    for x in range(0, len(data), 2):
        try: err_data = pd.read_csv(data[x+1]) 
        except: continue
        
        err_data = err_data[err_data[f_var] == f]
        err_data.reset_index(drop=True, inplace=True)
        
        plt.plot(err_data[x_var], err_data['Error'], colors[i])
        plt.title(title, fontsize=14)
        plt.xlabel(x_var, fontsize=14)
        plt.ylabel("Error [%]", fontsize=14)
        plt.grid(True)
    
        labels.append(data[x])
        i += 1
        
    plt.legend(labels, fontsize=12)

    
def plotTimeErrGraph(case1, opt1, case2, opt2, t0, tf):
    colors = ['r', 'g', 'b', 'm', 'y', 'c', 'k', 'w']
    labels = []
    i = 0
    
    plt.figure(figsize=(6, 6))
    plt.axis([t0, tf, 0.1, 10**19])
    
    for x in sorted(os.listdir(case1)):
        data1 = pd.read_csv(case1 + "/" + x + "/log_err.csv")
        
        plt.plot(data1['t'], data1['e'], colors[i])
        plt.yscale("log")
        plt.title("t x e", fontsize=14)
        plt.xlabel("Time [s]", fontsize=14)
        plt.ylabel("Error [%]", fontsize=14)
        plt.grid(True) 
        
        labels.append(x + opt1)
        i += 1
    
    labels.append(opt2)
    i = 0
    
    for x in sorted(os.listdir(case2)):
        data2 = pd.read_csv(case2 + "/" + x + "/log_err.csv")
        
        plt.plot(data2['t'], data2['e'], "s", 
                 color=colors[i], markeredgecolor='k')
        plt.yscale("log")
        plt.title("t x e", fontsize=14)
        plt.xlabel("Time [s]", fontsize=14)
        plt.ylabel("Error [%]", fontsize=14)
        plt.grid(True) 
        i += 1  
    
    plt.legend(labels, fontsize=12)
    
def plotRadialProfs(case, data, coeffs, data_y, dns_y, axes, y_label, save_path):
    plt.figure(figsize=(6, 6))      
    
    if data[0].empty == False: plt.plot(data[0]['x'], data[0][data_y], "m")             
    if data[6].empty == False: plt.plot(data[6]['x'], data[6][data_y], "g")    
    if data[12].empty == False: plt.plot(data[12]['x'], data[12][data_y], "b")   
    if data[18].empty == False: plt.plot(data[18]['x'], data[18][data_y], "r") 
    if data[3].empty == False: plt.plot(data[3]['r/R'], data[3][dns_y], "m--") 
    if data[9].empty == False: plt.plot(data[9]['r/R'], data[9][dns_y], "g--")  
    if data[15].empty == False: plt.plot(data[15]['r/R'], data[15][dns_y], "b--")    
    if data[21].empty == False: plt.plot(data[21]['r/R'], data[21][dns_y], "r--")

    plt.axis(axes)
    
    title = "r/R x " + data_y + " [Case " + case + "]\n"
    
    for x in range(0, len(coeffs), 2):  
        title += (coeffs[x] + "=" + str(coeffs[x+1]) + "  ")
    
    plt.title(title, fontsize=14)
    plt.xlabel("r/R", fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.legend(("z = 0D", "z = 5D", "z = 15D", "z = 27.5D", "DNS"), fontsize=12)
    plt.grid(True)
    plt.savefig(save_path + "/best_rR_" + data_y + ".png")

def plotAxialProfs(case, data, coeffs, data_y, dns_y, axes, y_label, save_path):
    plt.figure(figsize=(6, 6))      
    
    if data[24].empty == False: plt.plot(data[24]['z'], data[24][data_y], 
           "c", linewidth=3)
    if data[27].empty == False: plt.plot(data[27]['z/D'], data[27][dns_y], 
           "k", linewidth=3)

    plt.axis(axes)
    
    title = "z/D x " + data_y + " [Case " + case + "]\n"
    
    for x in range(0, len(coeffs), 2):  
        title += (coeffs[x] + "=" + str(coeffs[x+1]) + "  ")
    
    plt.title(title, fontsize=14)
    plt.xlabel("z/D", fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.legend(("MK", "DNS"), fontsize=12)
    plt.grid(True)
    plt.savefig(save_path + "/best_zD_" + data_y + ".png")
    
def testCoeffErrGraph():
    case = "42F"
    other_vars = ["Astar", 70, "Cd", 5]
    err_data = pd.read_csv("_42F_2_4_Prt=0.85_Astar/log_err.csv")
    plotCoeffErrGraph(case, "Astar", other_vars, err_data, 
                      "_42F_2_4_Prt=0.85_Astar")

def testCoeffErrGraphEx():
    data = ["37U", "../Sampling/_37U_2_4_Prt=0.85_Astar_Cd/log_err.csv",
            "38D", "../Sampling/_38D_2_4_Prt=0.85_Astar_Cd/log_err.csv",
            "39F", "../Sampling/_39F_2_4_Prt=0.85_Astar_Cd/log_err.csv",
            "40U", "../Sampling/_40U_2_4_Prt=0.85_Astar_Cd/log_err.csv",
            "41D", "../Sampling/_41D_2_4_Prt=0.85_Astar_Cd/log_err.csv",
            "42F", "../Sampling/_42F_2_4_Prt=0.85_Astar_Cd/log_err.csv"]

    plotCoeffErrGraphEx(data, "Astar", "Cd", 6.5, 50, 600, 0, 2)  
    
def testTimeErrGraph():
    plotTimeErrGraph("../Optimization/Co=0.01", " (Co=0.01)",
                     "../Optimization/Co=0.1", "Co=0.1",
                     1, 5)
    
def testProfGraph():
    case = "42F"
    dns_path = "../DNS/42F"
    case_path = "../Sampling/_42F_2_4_Prt=0.85_Astar/_Astar=70"
    coeffs = ["Astar", 70, "Cd", 5]
    tau_w0 = 0.076
    qw = -61.74
    
    #
    # Read the DNS data.
    #
    dns_d0 = pd.read_csv(dns_path + "/DNS_" + case + "_z=0D.csv")
    dns_d5 = pd.read_csv(dns_path + "/DNS_" + case + "_z=5D.csv")
    dns_d15 = pd.read_csv(dns_path + "/DNS_" + case + "_z=15D.csv")
    dns_d27 = pd.read_csv(dns_path + "/DNS_" + case + "_z=27D.csv")
    dns_wall = pd.read_csv(dns_path + "/DNS_" + case + "_x.csv")
    
    #
    # Convert the DNS TKE values to relative values.
    #    
    dns_d0['TKE'] = dns_d0['TKE'] / tau_w0
    dns_d5['TKE'] = dns_d5['TKE'] / tau_w0
    dns_d15['TKE'] = dns_d15['TKE'] / tau_w0
    dns_d27['TKE'] = dns_d27['TKE'] / tau_w0
    
    #
    # Read the sample's data.
    #
    d0 = pd.DataFrame()
    d5 = pd.DataFrame()
    d15 = pd.DataFrame()
    d27 = pd.DataFrame()
    wall = pd.DataFrame()
    
    for y in sorted(os.listdir(case_path)):
        frame = pd.read_csv(case_path + "/" + y)

        if y.startswith("sampleLineD0_") == True: 
            d0 = frame if d0.empty == True else pd.merge(d0, frame, on='x')
        elif y.startswith("sampleLineD5_") == True: 
            d5 = frame if d5.empty == True else pd.merge(d5, frame, on='x')
        elif y.startswith("sampleLineD15_") == True: 
            d15 = frame if d15.empty == True else pd.merge(d15, frame, on='x')
        elif y.startswith("sampleLineD27_") == True: 
            d27 = frame if d27.empty == True else pd.merge(d27, frame, on='x')
        elif y.startswith("sampleLineWall_") == True: 
            wall = frame if wall.empty == True else pd.merge(wall, frame, on='z')

    #
    # Invert the data rows, so the values match the DNS data.   
    #
    d0 = d0.iloc[::-1]
    d0.reset_index(drop=True, inplace=True)
    d5 = d5.iloc[::-1]
    d5.reset_index(drop=True, inplace=True)
    d15 = d15.iloc[::-1]
    d15.reset_index(drop=True, inplace=True)
    d27 = d27.iloc[::-1]
    d27.reset_index(drop=True, inplace=True)
    
    #
    # Convert the sample's radii values to relative values and delete
    # negative z values.
    #    
    d0['x'] = d0['x'] / 0.001 # R=0.001m
    d5['x'] = d5['x'] / 0.001 # R=0.001m
    d15['x'] = d15['x'] / 0.001 # R=0.001m
    d27['x'] = d27['x'] / 0.001 # R=0.001m
    wall['z'] = wall['z'] / 0.002 # D=0.002m
    
    wall = wall[wall['z'] >= 0]
    wall.reset_index(drop=True, inplace=True)
 
    #
    # Convert the sample's TKE values to relative values.
    #    
    d0['k'] = d0['k'] * d0['rho'] / tau_w0
    d5['k'] = d5['k'] * d5['rho'] / tau_w0
    d15['k'] = d15['k'] * d15['rho'] / tau_w0
    d27['k'] = d27['k'] * d27['rho'] / tau_w0
    
    data = [d0, 'x', 'T', dns_d0, 'r/R', 'Tmean', 
            d5, 'x', 'T', dns_d5, 'r/R', 'Tmean', 
            d15, 'x', 'T', dns_d15, 'r/R', 'Tmean',
            d27, 'x', 'T', dns_d27, 'r/R', 'Tmean',
            wall, 'z', 'T', dns_wall, 'z/D', 'Tw']
    
    plotRadialProfs(case, data, coeffs, 'U_2', 'UzMean', [0, 1, -0.1, 0.5], 
                    "Velocity [m/s]", ".")
    plotRadialProfs(case, data, coeffs, 'T', 'Tmean', [0, 1, 260, 350], 
                    "Temperature [K]", ".")
    plotAxialProfs(case, data, coeffs, 'T', 'Tw', [0, 35, 250, 350], 
                   "Wall temperature [K]", ".")
    plotRadialProfs(case, data, coeffs, 'k', 'TKE', [0, 1, 0, 20], 
                    "k*ρ/τw,0 [-]", ".")

def testProfGraph2():
    case = "50F"
    dns_path = "../DNS/50F"
    case_path = "../Sampling/_50F_2_4_Prt=0.85_Astar/_Astar=775"
    coeffs = ["Astar", 775, "Cd", 5]
    
    dns_wall = pd.read_csv(dns_path + "/DNS_" + case + "_x.csv")
    wall = pd.read_csv(case_path + 
                       "/sampleLineWall_T_h_k_product_ass_k_buoyancy_ass_k_alphat_myPrt_rho.csv")
    
    wall['z'] = wall['z'] / 0.002 # D=0.002m
    
    data = [0, 'x', 'T', 0, 'r/R', 'Tmean', 
            0, 'x', 'T', 0, 'r/R', 'Tmean', 
            0, 'x', 'T', 0, 'r/R', 'Tmean',
            0, 'x', 'T', 0, 'r/R', 'Tmean',
            wall, 'z', 'T', dns_wall, 'z/D', 'Tw']   

    plotAxialProfs(case, data, coeffs, 'T', 'Tw', [0, 35, 250, 350], 
                   "Wall temperature [K]", ".")   
#
# For testing.
#
#testProfGraph()
#testProfGraph2()
#testTimeErrGraph()
#testCoeffErrGraphEx()

