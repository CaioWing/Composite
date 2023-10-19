import composites as composite
import fail as criteria

if __name__ == '__main__':
    material_data = {'E11' : 19.76*10**9, 
                    'E22' : 1.97*10**9, 
                    'v12' : 0.35, 
                    'G12' : 0.7*10**9}
    
    test = composite.Composite(t = 3, 
                     angles = [0, 45, 0, 45, 0, 45, 45, 0, 45, 0, 45, 0],
                     normal= [1000/10**-3, 200, 0],
                     **material_data)
            
    test.run()
    criteria.max_stress(blades= test.blades_data, 
                        data= {"max_stress" : [963*10**6, -856*10**6, 900*10**6, -900*10**6, 71*10**6]})