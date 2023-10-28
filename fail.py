import numpy as np

composite_stress = ['sigma1', 'sigma2', 'sigma12']
composite_deform = ['epsilon1', 'epsilon2', 'gamma12']
names_stress = ['Xt', 'Xc', 'Yt', 'Yc', 'S12']
names_deformation = ["Xt'", "Xc'", "Yt'", "Yc'", "S12'"]

def print_message(status : bool, criteria : str, MS : float = None):
    """
    The function "print_message" prints a message indicating the status of a given criteria.
    
    :param status: A boolean value indicating the status of the check. True means the check passed,
    while False means the check failed
    :type status: bool
    :param criteria: A string that represents the criteria being checked. It could be anything that
    needs to be evaluated or verified
    :type criteria: str
    """
    print("Check: {} -> {}. MS: {}".format(criteria, status, MS))

# data : [Xt, Xc, Yt, Yc, S12] or/and [Xt', Xc', Yt', Yc', S12']  

def calculate_MS(blades : dict, 
                 data : dict = {},
                 sigma_name : str = "Local sigma"):
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
        if 'max_stress' in list(data.keys()):
            if blades[blade][sigma_name][0] > 0:
                condition = blades[blade][sigma_name][0] < data['max_stress']["Xt"]
                print_message(condition,
                        "Sigma1" + " < " + "Xt",
                        MS = data['max_stress']["Xt"]/blades[blade][sigma_name][0])
            else:
                condition = blades[blade][sigma_name][0] < data['max_stress']["Xc"]
                print_message(condition,
                        "Sigma1" + " > " + "Xc",
                        MS = data['max_stress']["Xc"]/blades[blade][sigma_name][0])
                
            if blades[blade][sigma_name][1] > 0:
                condition = blades[blade][sigma_name][1] < data['max_stress']["Yt"]
                print_message(condition,
                        "Sigma2" + " < " + "Yt",
                        MS = data['max_stress']["Yt"]/blades[blade][sigma_name][1])
            else:
                condition = blades[blade][sigma_name][1] > data['max_stress']["Yc"]
                print_message(condition,
                        "Sigma2" + " > " + "Yc",
                        MS = data['max_stress']["Yc"]/blades[blade][sigma_name][1])
                
            condition = abs(blades[blade][sigma_name][2]) <= data['max_stress']["S12"]
            print_message(condition,
                    "Sigma12" + " >= " + "S12",
                    MS = data['max_stress']["S12"]/abs(blades[blade][sigma_name][2]))

        if 'max_deformation' in list(data.keys()):
            if blades[blade]['Deformation'][0] > 0:
                condition = abs(blades[blade]['Deformation'][0]) < data['max_deformation']["Xt'"]
                print_message(condition,
                        "Epsilon1" + " < " + "Xt'")
            else:
                condition = blades[blade]['Deformation'][0] < data['max_deformation']["Xc'"]
                print_message(condition,
                        "Epsilon1" + " > " + "Xc'")
                
            if blades[blade]['Global sigma'][1] > 0:
                condition = blades[blade]['Deformation'][1] < data['max_deformation']["Yt'"]
                print_message(condition,
                        "Epsilon2" + " < " + "Yt'")
            else:
                condition = blades[blade]['Deformation'][1] > data['max_deformation']["Yc'"]
                print_message(condition,
                        "Epsilon2" + " > " + "Yc'")
                
            condition = abs(blades[blade]['Global sigma'][2]) <= data['max_deformation']["S12'"]
            print_message(condition,
                    "Epsilon12" + " >= " + "S12'")

def tsai_wo(blades : dict, data : dict = {}, sigma_name : str = "Local sigma"):
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
    f1 = 1/data['max_stress']['Xt'] + 1/data['max_stress']['Xc']
    f11 = -1/(data['max_stress']['Xt']*data['max_stress']['Xc'])
    f2 = 1/data['max_stress']['Yt'] + 1/data['max_stress']['Yc'] 
    f22 = -1/(data['max_stress']['Yt']*data['max_stress']['Yc'])
    f66 = 1/data['max_stress']['S12']**2
    f12 = -1/2*np.sqrt(f11*f22)
    
    print('\n Tsai Wo criteria: \n')
    for blade in blades:
        A = (f11*blades[blade][sigma_name][0]**2 + f22*blades[blade][sigma_name][1]**2 +
             f66*blades[blade][sigma_name][2]**2 + 2*f12*blades[blade][sigma_name][0]*blades[blade][sigma_name][1])
        B = f1*blades[blade][sigma_name][0] + f2*blades[blade][sigma_name][1]
        
        Sf_plus = (-B + np.sqrt(B**2 + 4*A))/(2*A)
        Sf_minus = abs((-B - np.sqrt(B**2 + 4*A))/(2*A))
        MS = max(Sf_minus, Sf_plus) - 1
        
        print("{} -> MS factor: {}".format(blade, MS))