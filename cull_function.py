import numpy as np
sigma = 4
fnType = 'Power'
A = 0.095
B = 0
C = 0
decision = 'LEN'
thresh = 50

def cull_function(template, decision_string = None, p = (A,B,C)):
    if fnType == 'Random':
        return np.random.normal()
        
    if decision_string is None: decision_string = decision
    if '+' in decision_string:
        decision_list = decision.split('+')
        return sum([cull_function(template, dec_str, p) for dec_str in decision_list])
    
    # calc the impact of decision
    substr = ''.join([b for b in decision if b in 'ACGU'])
    if decision == 'LEN':
        impact = len(template)
    elif decision[0] == 'f': # find first
        impact = template.find(substr)
    elif decision[0] == 'l': # find last
        impact = len(template) - template.rfind(substr)
    else:
        impact = template.count(substr)
    if '%' in decision:
        impact = impact/len(template)
    if '$' in decision:
        impact = abs(impact)
            
    #calc output value    
    if fnType == 'Linear':
        return p[0]+p[1]*impact
    elif fnType == 'Quadratic':
        return p[2]*impact**2+p[1]*impact+p[0]
    elif fnType == 'Power':
        return impact**p[0]
    else:
        return 0
        