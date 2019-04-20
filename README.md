# CCC-Model-Optimizer
Optimize Your Kinetic Model with One Click

## Requirments:  
  1. Cantera  
  2. Python3  
  3. Pandas  

## Files:  
  Mech_modifier.py:        main funtion  
  utils.py:                utilities  
  write_reactions.py:      wirte optimized model  
  soln2ck.py:              convert solution to chemkin mechanism  
  soln2cti.py:             convert solution to cantera mechanism  
  input_IDTs:              input ignition experimental data  
  input_uncertainties:     input uncertainties for each reaction  
  nC12-PAH_mech：           test mechanism  
  
## Results:

### 1. Using only uniform sampling  
 <img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/10_iterations.png" height="300"><img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/20_iterations.png" height="300">  
 <img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/30_iterations.png" height="300"><img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/40_iterations.png" height="300">  

### 2. Using only uniform sampling   
Uniform search is slow for fine tuning    
<img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/uniform_search.png" height="300">   

Combine uniform sampling with Gaussian sampling  
<img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/uniform_and%20gaussian_search.png" height="300">   
### 3. Using only uniform and Gaussian sampling 
<img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/10_iterations_with_gaussian.png" height="300"><img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/compare_uniform_and_gaussian_search.png" height="300">   
### 4. Data visualization
For data from different sources  
<img src="https://github.com/fyang235/CCC-Model-Optimizer/blob/master/Mech_modifier_v3/Images/data_visulization.png" height="300">   

## Enjoy!

