# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 13:50:55 2020

@author: Lydia Meyer
"""


import numpy as np
from matplotlib import pyplot as plt
from math import exp, log

F = 96485e3    # Faraday's constant, C/kmol equivalent charge
R = 8314.5     # Gas constant, J/kmol-K
T=298          #
k_fwd_star = 4.16307062e+7 # Chemical forward rate constant, m^4/kmol^2/s

C_elyte = 46.05    # Total (reference) elyte concentration, kmol/m3
C_Ni_s = 2.6e-08   # Total (reference) concentration of Ni surface sites, kmol/m2

X_H_Ni = 0.6
X_H2O_Ni = 0.2
X_Vac_Ni = 0.2
X_Vac_elyte = 0.08
X_Ox_elyte = 0.92

"Species standard-state thermo"
g_H_Ni_o = -7.109209e+07      # standard-state gibbs energy for H adsorbed on Ni surface (J/kmol)
g_H2O_Ni_o = -3.97403035e+08  # standard-state gibbs energy for H2O adsorbed on Ni surface (J/kmol)
g_Vac_Ni_o = 0.0              # standard-state gibbs energy for Ni surface vacancy (J/kmol)
g_Vac_elyte_o = 0.0           # standard-state gibbs energy for electrolyte oxide vacancy (J/kmol)
g_Ox_elyte_o = -2.1392135e+08 # standard-state gibbs energy for electrolyte oxide O2- (J/kmol)



# Loop over these delta phi = phi_anode - phi_elyte values:
phi_anode=0;
phi_elyte=phi_anode+delta_phi_dl_anode
phi_cathode=phi_elyte
delta_phi = np.linspace(-0.9,0.05,100)
i_elementary = np.zeros_like(delta_phi)

delta_g_o=((g_Vac_elyte_o+g_Vac_Ni_o+g_H2O_Ni_o)-(2*g_H_Ni_o+g_Ox_elyte_o))
delta_g_rxn=delta_g_o+(R*T*log((X_Vac_elyte*X_H2O_Ni*X_Vac_Ni)/(X_Ox_elyte*X_H_Ni**2)))
c_term=((X_Vac_elyte*X_H2O_Ni*X_Vac_Ni)/(X_Ox_elyte*X_H_Ni**2))
k_rev_star=k_fwd_star/(exp(-delta_g_rxn/(R*T))*c_term)
beta=0.5
n=-2

for ind, E in enumerate(delta_phi):
    k_fwd=k_fwd_star*exp((-beta*n*F*delta_phi[ind])/(R*T))
    k_rev=k_rev_star*exp(((1-beta)*n*F*delta_phi[ind])/(R*T))
    i_elementary[ind] = n*F*(k_fwd*((C_elyte*X_Ox_elyte)*(C_Ni_s*X_H_Ni)**2)-k_rev*((C_elyte*X_Vac_elyte)*(C_Ni_s*X_H2O_Ni*X_Vac_Ni*C_Ni_s)))
    
fig1, ax1 = plt.subplots()
ax1.plot(delta_phi,i_elementary,linewidth = 1.5,color = 'k')
ax1.plot(E_validate,i_validate,'ro',linewidth = 1.5)
ax1.set_xlabel('Electric Potential [V]',family='Times New Roman',fontsize=14)
ax1.set_ylabel('Current [A/m2]',family='Times New Roman',fontsize=14)
fig2, ax2 = plt.subplots()
ax2.semilogy(delta_phi,abs(i_elementary),linewidth = 1.5,color = 'k')
ax2.semilogy(E_validate,abs(i_validate),'ro',linewidth = 1.5)
ax2.set_xlabel('Electric Potential [V]',family='Times New Roman',fontsize=14)
ax2.set_ylabel('Current [A/m2]',family='Times New Roman',fontsize=14)

i_BV = np.zeros_like(delta_phi)
i_o=n*F*(k_fwd_star**(1-beta))*(k_rev_star**beta)*(((C_elyte*X_Ox_elyte)*(C_Ni_s*X_H_Ni)**2)**(1-beta))*((C_elyte*X_Vac_elyte)*(C_Ni_s*X_H2O_Ni*X_Vac_Ni*C_Ni_s))**beta
for ind, E in enumerate(delta_phi):
    eta=delta_phi[ind]-(-delta_g_rxn/(n*F))
    i_BV[ind] = i_o*(exp((-beta*n*F*eta)/(R*T))-exp(((1-beta)*n*F*eta)/(R*T)))
plt.close('all')
plt.plot(delta_phi,i_elementary,linewidth = 1.5,color = 'k');
plt.plot(delta_phi,i_BV,'ro',markerfacecolor='none');
plt.xlabel('Overpotential [V]',family='Times New Roman',fontsize=14)
plt.ylabel('Current [A/m2]',family='Times New Roman',fontsize=14)
plt.legend(['Mass Action','Butler Volmer'],frameon=False)
plt.show()

i_Tafel = np.zeros_like(delta_phi)

for ind, E in enumerate(delta_phi):
    eta=delta_phi[ind]-(-delta_g_rxn/(n*F))
    i_Tafel[ind] =-i_o*exp((beta*n*F*eta)/(R*T))
plt.plot(delta_phi,i_elementary,linewidth = 1.5,color = 'k');
plt.plot(delta_phi,i_BV,'o',markeredgecolor='r',markerfacecolor='none');
plt.plot(delta_phi,i_Tafel,'x',markeredgecolor='b',markerfacecolor='none');
plt.xlabel('Overpotential [V]',family='Times New Roman',fontsize=14)
plt.ylabel('Current [A/m2]',family='Times New Roman',fontsize=14)
plt.legend(['Mass Action','Butler Volmer','Tafel'],frameon=False)
plt.show()


y=1
x=2
