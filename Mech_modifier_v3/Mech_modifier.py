from __future__ import division
from __future__ import print_function
import sys
sys.path.append(r'/home/yang/anaconda3/lib/python3.6/site-packages')

import cantera as ct
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import utils
import copy
import write_reactions
import soln2ck
import soln2cti
##############################################################################
# inputs
##############################################################################

INPUT_IDT = r'./input_IDTs'
INPUT_UNCERTAINTIES = r'./input_uncertainties'
INPUT_MECH = './nC12-PAH_mech/mech.cti'
INDICATOR = 'oh'

OUTPUT_REACTONS_FILE_NAME = './new_reactions'
OUTPUT_FACTORS_FILE_NAME = './factors'
OUTPUT_CK = True
OUTPUT_CTI = True

MIN_LOSS = 0.001
MAX_GAUSS = 5

# output info
print('===============================================================')
print('Inputs summary: ')
print("INPUT_IDT:                   ", INPUT_IDT)
print("INPUT_UNCERTAINTIES:         ", INPUT_UNCERTAINTIES)
print("INPUT_MECH:                  ", INPUT_MECH)
print("INDICATOR:                   ", INDICATOR)
print("OUTPUT_REACTONS_FILE_NAME:   ", OUTPUT_REACTONS_FILE_NAME)
print("OUTPUT_FACTORS_FILE_NAME:    ", OUTPUT_FACTORS_FILE_NAME)
print("OUTPUT_CK:                   ", OUTPUT_CK)
print("OUTPUT_CTI:                  ", OUTPUT_CTI)
print("MIN_LOSS:                    ", MIN_LOSS)
print("MAX_GAUSS:                   ", MAX_GAUSS)

# read information
info, idts = utils.sparse_exprimental_conditions(INPUT_IDT)
uncertainties = utils.sparse_uncertainties(INPUT_UNCERTAINTIES)

reactions = ct.Reaction.listFromFile(INPUT_MECH)
species = ct.Species.listFromFile(INPUT_MECH)

##############################################################################
# main loop
##############################################################################

# calc original predictions, loop through every gt exprimental point
original_taus = []
for ind, condition in info.items():
    gas = ct.Solution(INPUT_MECH)
    gas.TPX = condition
    reactor = ct.Reactor(contents=gas)
    reactorNetwork = ct.ReactorNet([reactor])
    
    # efficiency
    indicator_index = utils.get_indicator_index(reactor, INDICATOR)

    tau, delta_t = utils.calculate_idt(reactorNetwork, INDICATOR, indicator_index)
    original_taus.append(tau)
        
        
# modify mech, calc new predictions, compare with gt and save the best        
counter = 1

# place holders
best_reactions = None
best_factors = None
best_taus = None

# initialize loss with original loss
original_loss = utils.l2_loss(idts, original_taus)
loss = original_loss
print('===============================================================')
print('Original loss: %.5f' % loss)

# random loop
t0 = time.time()
while True:
    # generate a mech
    t_iteration_0 = time.time()
    reactions_var = ct.Reaction.listFromFile(INPUT_MECH)  
    
    # use Gaussian ditribution after 4 loops
    if counter <= 4:
        means = None
        print('Uniform search, ', end='')
    else:
        means = best_factors
        print('Gaussian search, ', end='')
        
    reactions_new, factors = utils.generate_new_reactions(reactions_var, uncertainties, means)
    
    # make predictions with new mech
    taus = []
    for ind, condition in info.items():
        gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                        species=species, reactions=reactions_new)
        gas.TPX = condition
        reactor = ct.Reactor(contents=gas)
        reactorNetwork = ct.ReactorNet([reactor])
        
        # efficiency
        indicator_index = utils.get_indicator_index(reactor, INDICATOR)

        tau, delta_t = utils.calculate_idt(reactorNetwork, INDICATOR, indicator_index)
        taus.append(tau)
    
    # calc idt loss
    l = utils.l2_loss(idts, taus)
 
    # keep the best
    if l < loss:
        loss = l
        best_reactions = reactions_new
        best_factors = factors
        best_taus = taus
        write_reactions.write_reactions(best_reactions, OUTPUT_REACTONS_FILE_NAME)
        write_reactions.write_factors(best_factors, OUTPUT_FACTORS_FILE_NAME)
        
        # use solution to save ck or cti file
        if OUTPUT_CK:
            gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                            species=species, reactions=best_reactions)            
            soln2ck.write(gas)
        if OUTPUT_CTI:
            gas = ct.Solution(thermo='IdealGas', kinetics='GasKinetics',
                            species=species, reactions=best_reactions)            
            soln2cti.write(gas)
     
    t_iteration_1 = time.time()    
    #print('Counter: %d, Loss: %.5f, Time cost: %.2f s' % (counter, loss, (t_iteration_1 - t_iteration_0)))
    num_space = 5 - len(str(counter))
    print('Counter:%s%d, Loss: %.5f, Time cost: %.2f s' % (num_space*' ', counter, loss, (t_iteration_1 - t_iteration_0)))

    # end metrics
    if loss <= MIN_LOSS or counter >= MAX_GAUSS:
        if loss == original_loss:
            print('Did not find a better mech, please increase MAX_GAUSS and retry.')
            exit()
        else:
            break

    counter += 1
t1 = time.time()
print('Total time cost: %.2f min' % (float(t1 - t0)/60.))
print('===============================================================')
##############################################################################
# plot results
##############################################################################
# plot tau-T figure, only valid when the input data share all the parameters but temperature.

# get Temperatures
Temperatures = [ condition[0] for ind, condition in info.items()]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.semilogy(1000./np.array(Temperatures), np.array(idts), 'o', color='b')
ax.semilogy(1000./np.array(Temperatures), np.array(original_taus), color='r')
ax.semilogy(1000./np.array(Temperatures), np.array(best_taus), color='g')
ax.legend(['Experimental data','Original predctions','New predictions'], loc='lower right')

ax.set_ylabel('Ignition Delay (s)', fontname='Times New Roman', fontsize=16)
ax.set_xlabel(r'$\mathdefault{1000/T\, (K^{-1})}$', fontname='Times New Roman', fontsize=16)

# Add a second axis on top to plot the temperature for better readability
ax2 = ax.twiny()
ticks = ax.get_xticks()
ax2.set_xticks(ticks)
ax2.set_xticklabels((1000/ticks).round(1))
ax2.set_xlim(ax.get_xlim())
ax2.set_xlabel('Temperature (K)',fontname='Times New Roman',fontsize=16)


# plot pred-true figure, valid for abitary inputs
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.loglog(np.array(idts), np.array(original_taus), 's', color='r')
ax.loglog(np.array(idts), np.array(best_taus), 'o', color='g')
ax.legend(['Exp-New','Exp-Old'], loc='lower right')

# set identity limits
low_lim = np.minimum(ax.get_xlim()[0], ax.get_ylim()[0])
high_lim = np.maximum(ax.get_xlim()[1], ax.get_ylim()[1])
lim = (low_lim, high_lim)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.plot(lim, lim, color='b')

# add labels
ax.set_ylabel('Experimental data', fontname='Times New Roman', fontsize=16)
ax.set_xlabel('Mech predictions',  fontname='Times New Roman', fontsize=16)

plt.show()

















