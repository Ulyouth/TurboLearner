#!/usr/bin/env python
"""
NAME:   optzSamples.py
AUTHOR: André Maia
DATE:   31.05.2020
DESC:   A script to compare the evolution of a sample's results over 
        time and determine when they stop varying (when the flow becomes
        fully-developed.

USAGE:  optzSamples.py [CASE] [START] [END]
        Where:
        CASE   is a folder containing the OpenFOAM sample.
        START  is the initial time to be considered.
        END    is the final time to be considered (used as control).
"""

import sys
import os
import matplotlib.pyplot as plt
import pandas as pd

case = sys.argv[1]
t0 = float(sys.argv[2])
tf = float(sys.argv[3])

ctrl_path = case + "/postProcessing/sampleDict/" + sys.argv[3]
print("Input case: " + case + "\n")  

res = []
err = []
t = []
    
# List all folders in the case folder.
for x in sorted(os.listdir(case)):   
    # Check if the current folder is a time folder.
    if x.replace('.', '', 1).isdigit() == False:
        continue
    
    # Check if the time of the current folder is not smaller than 
    # the start time
    if float(x) < t0:
        continue
    
    # Check if the current time folder is not the control sample.
    if float(x) == tf:
        continue
    
    smpl_path = case + "/postProcessing/sampleDict/" + x
    e = 0
    count = 0
    
    # List all archives in the post-process folder.
    for y in sorted(os.listdir(smpl_path)):
        smpl = pd.read_csv(smpl_path + "/" + y)
        ctrl = pd.read_csv(ctrl_path + "/" + y)
        
        for i in smpl: 
            # Check if the current column is not a coordinate.
            if i == 'x' or i == 'z':
                continue
            
            for j in range(len(smpl)):
                # Load the CSV data of the current & control samples.
                smpl_data = float(smpl[i][j])
                ctrl_data = float(ctrl[i][j])

                # Calculate the error of the current sample.
                delta = abs((smpl_data - ctrl_data))
                base = (1 if ctrl_data == 0 else abs(ctrl_data))
                e += ((delta / base) * 100)
                count += 1
                
    e /= count
    
    err.append(e)
    t.append(float(x))
    res.append([float(x), e])
    
    print("t=" + str(x) + "  e=" + str(e))   

# Plot the resulting graph
plt.figure(figsize=(6, 6))
plt.plot(t[:], err[:])
plt.yscale("log")
plt.axis([1, 5, 0.1, max(err)])
plt.title("t x e [Case " + case + "]", fontsize=14)
plt.xlabel("Time [s]", fontsize=14)
plt.ylabel("Error [%]", fontsize=14)
plt.grid(True) 

# Saves the error data to a CSV file.
log_csv = pd.DataFrame(res, columns=["t", "e"])
log_csv.sort_values(by=["t", "e"], inplace=True)
log_csv.to_csv(case + "/log_err.csv", encoding="utf-8", index=False)
    