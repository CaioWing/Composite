import composites as composite
import fail as criteria

if __name__ == '__main__':
    composite_data = {"properties" : {'E11' : 155*10**9, 
                                    'E22' : 12.1*10**9, 
                                    'v12' : 0.35, 
                                    'G12' : 4.4*10**9},
                    "max_stress" : {"Xt" : 250*10**6, 
                                    "Xc" : -180*10**6,
                                    "Yt" : 40*10**6,
                                    "Yc" : -45*10**6,
                                    "S12" : 47*10**6}
                    # "max_deformation": {"Xt'" : 0.020,
                    #                 "Xc'" : -0.018,
                    #                 "Yt'" : 0.007,
                    #                 "Yc'" : -0.012,
                    #                 "S12'" : 0.01,}
                    }
    custom_layer = [
        {
            "name" : "nucleo",
            "thickness" : 30,
            "properties" : {'E11' : 75*10**6, 
                            'E22' : 75*10**6, 
                            'v12' : 0.875, 
                            'G12' : 20*10**6},
            "max_stress" : {"Xt" : 1.8*10**6, 
                            "Xc" : 0.9*10**6,
                            # "Yt" : 40*10**6,
                            # "Yc" : -45*10**6,
                            "S12" : 0.76*10**6}
        }
    ]
    test = composite.Composite(t = 0.1, 
                     angles = [0, 45, 45, 0],
                     normal= [1717.6/10**-3, 3435.2/10**-3, 0],
                     momentum= [100, 0, 0],
                     custom_layers= custom_layer,
                     **composite_data)
    
    # material_data = {'E11' : 19.76*10**9, 
    #                 'E22' : 1.97*10**9, 
    #                 'v12' : 0.35, 
    #                 'G12' : 0.7*10**9}
    
    # test = composite.Composite(t = 3, 
    #                  angles = [45, 0, 0, 45],
    #                  normal= [1000/10**-3, 200/10**-3, 0],
    #                  **material_data)

    test.run()
    test.save_data()
    print(test.h)
    # criteria.calculate_MS(blades= test.blades_data, 
    #                      data = composite_data)
    
    # criteria.tsai_wo(blades= test.blades_data, 
    #                     data= composite_data)
    
    # criteria.tsai_hill(blades= test.blades_data, 
    #                     data= composite_data)