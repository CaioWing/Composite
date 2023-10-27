import numpy as np

composite_stress = ['sigma1', 'sigma2', 'sigma12']
composite_deform = ['epsilon1', 'epsilon2', 'gamma12']
names_stress = ['Xt', 'Xc', 'Yt', 'Yc', 'S12']
names_deformation = ["Xt'", "Xc'", "Yt'", "Yc'", "S12'"]

def print_message(status : bool, criteria : str):
    """
    The function "print_message" prints a message indicating the status of a given criteria.
    
    :param status: A boolean value indicating the status of the check. True means the check passed,
    while False means the check failed
    :type status: bool
    :param criteria: A string that represents the criteria being checked. It could be anything that
    needs to be evaluated or verified
    :type criteria: str
    """
    print("Check: {} -> {}".format(criteria, status))

# data : [Xt, Xc, Yt, Yc, S12] or/and [Xt', Xc', Yt', Yc', S12']  

def max_stress(blades : dict, data : dict = {}):
    """
    The function `max_stress` checks if stress and deformation values for blades meet certain criteria.
    
    :param blades: The `blades` parameter is a dictionary that contains information about different
    blades. Each blade is represented by a key in the dictionary, and the value associated with each key
    is another dictionary that contains specific data about that blade
    :type blades: dict
    :param data: The `data` parameter is a dictionary that contains the maximum stress and maximum
    deformation values for each blade. It is an optional parameter with a default value of an empty
    dictionary
    :type data: dict
    """
    for blade in blades:
        print('\n', blade)
        for i in [0, 2, 4]:
            if 'max_stress' in list(data.keys()):
                print_message(blades[blade]['Global sigma'][i//2] < data['max_stress'][i],
                        composite_stress[i//2] + " < " + names_stress[i])
                print("MS value: {}".format(abs(data['max_stress'][i]/blades[blade]['Global sigma'][i//2]) - 1))

            if 'max_deformation' in list(data.keys()):
                print_message(blades[blade]['Deformation'][i//2] < data['max_deformation'][i],
                        composite_deform[i//2] + " < " + names_deformation[i])
                print("MS value: {}".format(abs(data['max_deformation'][i]/blades[blade]['Deformation'][i//2]) - 1))

def tsai_wo(blades : dict, data : dict = {}):
    """
    The `tsai_wo` function calculates the MS (Marginal Safety) factor for each blade in a given
    dictionary of blades, based on their local stress values and maximum stress values provided in a
    data dictionary.
    
    :param blades: The `blades` parameter is a dictionary that contains information about different
    blades. Each blade is represented by a key in the dictionary, and the corresponding value is another
    dictionary that contains the following information:
    :type blades: dict
    :param data: The `data` parameter is a dictionary that contains the following keys and values:
    :type data: dict
    """
    f1 = 1/data['max_stress'][0] + 1/data['max_stress'][1]
    f11 = -1/(data['max_stress'][0]*data['max_stress'][1])
    f2 = 1/data['max_stress'][2] + 1/data['max_stress'][3] 
    f22 = -1/(data['max_stress'][2]*data['max_stress'][3])
    f66 = 1/data['max_stress'][4]**2
    f12 = -1/2*np.sqrt(f11*f22)
    
    print('\n Tsai Wo criteria: \n')
    for blade in blades:
        A = (f11*blades[blade]['Global sigma'][0]**2 + f22*blades[blade]['Global sigma'][1]**2 +
             f66*blades[blade]['Global sigma'][2]**2 +
             2*f12*blades[blade]['Global sigma'][0]*blades[blade]['Global sigma'][1])
        B = f1*blades[blade]['Global sigma'][0] + f2*blades[blade]['Global sigma'][1]
        
        Sf_plus = (-B + np.sqrt(B**2 + 4*A))/(2*A)
        Sf_minus = abs((-B - np.sqrt(B**2 + 4*A))/(2*A))
        MS = 1 - max(Sf_minus, Sf_plus)
        
        print("{} -> MS factor: {}".format(blade, MS))