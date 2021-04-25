#!/usr/bin/env python
"""
Created on Tue May  5 16:22:46 2020

@author: andre
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp

m_T = np.zeros([201, 5])
m_cp = np.zeros([201, 5])
m_d = np.zeros([201, 5])

for p in range(80, 121, 10):
    for T in range(250, 451):
        x = T-250
        y = ((p-70)//10)-1
        
        m_T[x,y] = T
        m_cp[x,y] = cp.PropsSI("C", "T", T, "P", p*10**5, "CO2") / 1000
        m_d[x,y] = cp.PropsSI("D", "T", T, "P", p*10**5, "CO2")
        
plt.figure(figsize=(6, 6))
plt.plot(m_T[:,0], m_cp[:,0], "r",
         m_T[:,1], m_cp[:,1], "g",
         m_T[:,2], m_cp[:,2], "b",
         m_T[:,3], m_cp[:,3], "m",
         m_T[:,4], m_cp[:,4], "c")

plt.title("cp x T graph for CO2", fontsize=14)
plt.xlabel("Temperature [K]", fontsize=14)
plt.ylabel("Specific Heat [kJ/kg]", fontsize=14)
plt.legend(("p = 80 bar", "p = 90 bar", "p = 100 bar", "p = 110 bar", "p = 120 bar"), fontsize=14)
plt.grid(True)
plt.show()

plt.figure(figsize=(6, 6))
plt.plot(m_T[:,0], m_d[:,0], "r",
         m_T[:,1], m_d[:,1], "g",
         m_T[:,2], m_d[:,2], "b",
         m_T[:,3], m_d[:,3], "m",
         m_T[:,4], m_d[:,4], "c")

plt.title("ρ x T graph for CO2", fontsize=14)
plt.xlabel("Temperature [K]", fontsize=14)
plt.ylabel("Density [kg/m³]", fontsize=14)
plt.legend(("p = 80 bar", "p = 90 bar", "p = 100 bar", "p = 110 bar", "p = 120 bar"), fontsize=14)
plt.grid(True)
plt.show()
    
