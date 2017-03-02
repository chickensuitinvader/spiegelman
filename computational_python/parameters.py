"""
sptParameters
Default Parameters for Simulation
"""

parameters = {
    # for SpSim (base)
    'InitialTemplates' : 1000,
    'TransferPercent' : 10,
    'SeedLength' : (1000, 50), # Mean Length, Std. Deviation
    'MinLength' : 50,
    'Cycles' : 400,
    'MaxReplications' : 100000,
    'Replicators' : 10,
    'Pairings' : str.maketrans({'A':'U', 'U':'A', 'C':'G', 'G':'C'}),
    'InitialPool' : {'A':200000, 'C':200000, 'G':200000, 'U':200000},
    'EmptyPool' : 0.001,
    'PointMutations' : {'substitution' : 0.01, 'addition' : 0.01, 'deletion' : 0.01},
    'BlockMutations' : {'addition' : 5e-6, 'deletion' : 5e-6},
    # for SptSim (Spatial Functional Behaviour)
    'Epochs' : 400,
    'ShufflePercent' : 100,
    'ShuffleType': 'ComplexPref',
    'Size' : (2,2),
    'Changes' : ('+A','+G','+C','+U'),
    'ChangeValue' : 1.5
    }
