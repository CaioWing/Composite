
# Composite calculator

These scripts can describe the behavior of composite materials that have multiple layers, it uses the following steps:

    The `Composite` class represents a composite material with multiple layers and calculates various
    properties and stress distributions within the material.
    
    Step 1: Calculate the Reduced Stiffness Matrix of a laminate
    with respect to the Local Coordinate System

    Step 2: Calculate the Transformed Reduced Stiffness Matrix of each laminate
    with respect to the Global Coordinate System

    Step 3: Calculate the sub-matrices A, B, and D of the stiffness matrix of the laminate
    with respect to the Global Coordinate System

    Step 4: Calculate mid-plane strains and curvatures
    with respect to the Global Coordinate System

    Step 5: Calculate the stresses acting on a given layer k
    with respect to the Global Coordinate System

    Step 6: Calculate the stresses acting on a given layer k
    with respect to the Local Coordinate System



## Example

The code just need the specification of the layer's angle, the material data and the external influence.

    material_data = {'E11' : 19.76*10**9, 
                    'E22' : 1.97*10**9, 
                    'v12' : 0.35, 
                    'G12' : 0.7*10**9}
    
    test = Composite(t = 3, 
                     angles = [45, 0, 0, 45],
                     normal= [1000/10**-3, 200/10**-3, 0],
                     **material_data)
    
    test.run()
    test.save_data()
## Install

Just need to clone the project! 

```bash
  git clone https://github.com/CaioWing/Composite
```
    
