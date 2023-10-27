import composites as composite
import fail as criteria

if __name__ == '__main__':
    # material_data = {'E11' : 77*10**9, 
    #                 'E22' : 75*10**9, 
    #                 'v12' : 0.06, 
    #                 'G12' : 6.5*10**9}
    
    # test = composite.Composite(t = 0.29, 
    #                  angles = [0, 45, 0, 45, 0, 45, 45, 0, 45, 0, 45, 0],
    #                  normal= [1069.6701/10**-3, 0, 0],
    #                  **material_data)
    material_data = {'E11' : 19.76*10**9, 
                    'E22' : 1.97*10**9, 
                    'v12' : 0.35, 
                    'G12' : 0.7*10**9}
    
    test = composite.Composite(t = 3, 
                     angles = [45, 0, 0, 45],
                     normal= [1000/10**-3, 200/10**-3, 0],
                     **material_data)
            
    test.run()
    test.save_data()
    criteria.max_stress(blades= test.blades_data, 
                         data= {"max_stress" : [963*10**6, -856*10**6, 900*10**6, -900*10**6, 71*10**6]})
    
    criteria.tsai_wo(blades= test.blades_data, 
                        data= {"max_stress" : [963*10**6, -856*10**6, 900*10**6, -900*10**6, 71*10**6]})