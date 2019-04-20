from __future__ import division
from __future__ import print_function
import sys
sys.path.append(r'/home/yang/anaconda3/lib/python3.6/site-packages')

import cantera as ct
import numpy as np

# note the mech is in kcal/mol but cantera reactions are in J
J_to_KCAL = 4184
THOUSAND = 1000

def dict2str(dic):
    '''
    convert dictionary to string
    '''
    diclist = ['{}:{}'.format(k, v) for k, v in dic.items()]
    string = ''
    for i in diclist:
        string += ' '+i
    
    string = string[1:]
    return string

def write_reactions(reactions, filename):
    with open(filename, 'w') as f:
        for i, r in enumerate(reactions):
            # elementary reaction
            if r.reaction_type == 1:
            # eg:
            ## Reaction 378
            #reaction('ch2 + o2 <=> co + h2o', [7.280000e+19, -2.54, 1809.03])
                f.write('# Reaction {:s} \n'.format(str(i + 1)))
                f.write('''raction('{:s}', [{:.3e}, {:3.2f}, {:3.2f}]'''.format(r.equation,
                                                                                r.rate.pre_exponential_factor*THOUSAND,
                                                                                r.rate.temperature_exponent,
                                                                                r.rate.activation_energy/J_to_KCAL))
                f.write(')\n\n')
            # three body reaction
            elif r.reaction_type == 2:
            # eg:
            # Reaction 384
            #three_body_reaction('h + co + M => hco + M', [6.467000e+13, 0.0, -441.92],
                                #efficiencies='co:1.9 co2:3.8 h2:2.5 h2o:12.0')
                f.write('# Reaction {:s} \n'.format(str(i + 1)))
                f.write('''three_body_reaction('{:s}', [{:.3e}, {:3.2f}, {:.3e}]'''.format(r.equation,
                                                                                r.rate.pre_exponential_factor*THOUSAND,
                                                                                r.rate.temperature_exponent,
                                                                                r.rate.activation_energy/J_to_KCAL))
                if len(r.efficiencies) != 0:
                    f.write(', \n' + 23*' ' + '''efficiencies=\'{:s}\''''.format(dict2str(r.efficiencies)))
                f.write(')\n\n')
                
            # fall-off and chemical reaction
            elif r.reaction_type == 4:
            # eg:
            ## Reaction 345         
            #falloff_reaction('h + c2h4 (+ M) <=> c2h5 (+ M)',
                            #kf=[1.081000e+12, 0.45, 1821.94],
                            #kf0=[1.112000e+34, -5.0, 4447.9],
                            #efficiencies='co:2.0 co2:3.0 h2:2.0 h2o:5.0',
                            #falloff=Troe(A=1.0, T3=1e-15, T1=95.0, T2=200.0))
                f.write('# Reaction {:s} \n'.format(str(i + 1)))
                if len(r.efficiencies) != 0:
                    f.write('''falloff_reaction('{:s}',
                            kf=[{:.3e}, {:3.2f}, {:.3e}],
                            kf0=[{:.3e}, {:3.2f}, {:.3e}]'''.format(r.equation,
                                                                        r.high_rate.pre_exponential_factor*THOUSAND,
                                                                        r.high_rate.temperature_exponent,
                                                                        r.high_rate.activation_energy/J_to_KCAL,
                                                                        r.low_rate.pre_exponential_factor*THOUSAND*THOUSAND,
                                                                        r.low_rate.temperature_exponent,
                                                                        r.low_rate.activation_energy/J_to_KCAL))
                else:
                    f.write('''falloff_reaction('{:s}',
                            kf=[{:.3e}, {:3.2f}, {:.3e}],
                            kf0=[{:.3e}, {:3.2f}, {:.3e}]'''.format(r.equation,
                                                                        r.high_rate.pre_exponential_factor,
                                                                        r.high_rate.temperature_exponent,
                                                                        r.high_rate.activation_energy/J_to_KCAL,
                                                                        r.low_rate.pre_exponential_factor*THOUSAND,
                                                                        r.low_rate.temperature_exponent,
                                                                        r.low_rate.activation_energy/J_to_KCAL))
                        
                if len(r.efficiencies) != 0:
                    f.write(', \n' + 23*' ' + '''efficiencies=\'{:s}\''''.format(dict2str(r.efficiencies)))
                    
                if len(r.falloff.parameters ) != 0:
                    f.write(', \n' + 23*' ' + \
                        '''falloff={:s}(A={:.3e}, T3={:.3e}, T1={:.3e}, T2={:.3e})'''.format(r.falloff.type,
                                                                                        r.falloff.parameters[0],
                                                                                        r.falloff.parameters[1],
                                                                                        r.falloff.parameters[2],
                                                                                        r.falloff.parameters[3],))        
                f.write(')\n\n')
            else:
                print('The raction type is {} !'.format(r.reaction_type))
            
def write_factors(factors, filename):
    with open(filename, 'w') as f:
        for k, v in factors.items():
            f.write('{:d}, {:3.2f} \n'.format(k, v))
            
            
            
            
            
            