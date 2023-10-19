composite_stress = ['sigma1', 'sigma2', 'sigma12']
composite_deform = ['epsilon1', 'epsilon2', 'gamma12']
names_stress_T = ['Xt', 'Yt', 'S12']
names_deformation_T = ["Xt'", "Yt'", "S12'"]
names_stress_C = ['Xc', 'Yc', 'S12']
names_deformation_C = ["Xc'", "Yc'", "S12'"]

def print_message(status : bool, criteria : str):
    print("Check: {} -> {}".format(criteria, status))

# data : [Xt, Xc, Yt, Yc, S12] or/and [Xt', Xc', Yt', Yc', S12']  

def max_stress(blades : dict, data : dict = {"max_deformation" : [],
                                            "max_stress" : []}):
    for blade in blades:
        print('\n', blade)
        for i in range(3):
            if 'max_stress' in list(data.keys()):
                print_message(blades[blade]['Local sigma'][i] <= data['max_stress'][i],
                            composite_stress[i] + " <= " + names_stress_T[i])
                print(blades[blade]['Local sigma'][i], data['max_stress'][i])
            if 'max_deformation' in list(data.keys()):
                print_message(blades[blade]['Deformation'][i] <= data['max_deformation'][i],
                            composite_deform[i] + " <= " + names_deformation_T[i])
        