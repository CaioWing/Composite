import numpy as np
import json

class Composite():
    """
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
        
    Returns:
        _dict_ `blades_data`: have the global and local stresses of each layers 
    """
    A = np.zeros((3,3))
    B = np.zeros((3,3))
    D = np.zeros((3,3))
    sigma_global = []
    sigma_local = []
    blades_data = {}
    system_properties = {}
    system_properties = {}
    blades_data = {}
    Q_ = []
    epsilon_local = []
    layers = []
    angles = []
    
    def __init__(self,
                 custom_layers : dict = None,
                 t : float = 3,
                 angles = [45, 0, 0, 45], 
                 normal : list = [0, 0, 0], 
                 momentum : list = [0, 0, 0],
                 **composite_data) -> None:
        """
        The function initializes an object with various attributes and assigns values to them.
        
        :param t: The parameter `t` is a float that represents the thickness of each layer in [mm]. It has
        a default value of 3, defaults to 3
        :type t: float (optional)
        :param angles: The `angles` parameter is a list that represents the angles of incidence for each
        layer in a multilayer system. The angles are specified in degrees
        :param normal: The `normal` parameter is a list that represents the normal stress aplied to the structure [N/m].
        It has three elements, where each element represents the x, y, and xy components of the normal vector,
        respectively
        :type normal: list
        :param momentum: The `momentum` parameter is a list that represents the momentum aplied to the structure [N/mm]. It has 
        three components, corresponding to the momentum in the x, y, and xy directions respectively
        :type momentum: list
        """
        self.custom_layers = custom_layers
        self.t = t/1000
        self.N = normal
        self.M = momentum
        self.composite_angles = angles
        self.n_layers = len(angles) if custom_layers is None else len(angles) + len(custom_layers) 
        self.composite_data = composite_data
        self.calculate_layer_thicknesses()

    def calculate_layer_thicknesses(self):
        half_layers = self.n_layers // 2
        if self.custom_layers is None:
            self.h = [-(half_layers) * self.t + i * self.t for i in range(half_layers + 1)]
        else:
            self.h = [-(half_layers) * self.t + i * self.t - self.custom_layers[0]['thickness']/2 for i in range(half_layers + 1)]
        self.layers = [{"name": "blade " + str(i), "properties": self.composite_data["properties"]} for i in range(len(self.h))]
        self.angles = self.composite_angles[:half_layers]

        if self.custom_layers is not None:
            for layer in self.custom_layers:
                position = layer.get("position", "center")  # Assume "center" if not specified
                thickness = layer.get("thickness", self.t)  # Assume self.t if not specified
                angle = layer.get("angle", 0)
                if position == 'center':
                    self.h.append(self.h[-1] + thickness)
                    self.layers.append({"name": layer["name"], "properties": layer["properties"]})
                    self.angles.append(angle)

        # Corrigindo a linha abaixo
        self.h.extend([self.h[-1] + i * self.t for i in range(1, half_layers + 1)])
        self.layers.extend([{"name": "blade " + str(i), "properties": self.composite_data["properties"]} for i in range(half_layers + 1, self.n_layers + 1)])
        self.angles.extend(self.composite_angles[half_layers:])

# Certifique-se de que os ângulos na segunda metade sejam pegos a partir da metade da lista original


                
    def calc_RRT(self, E11 : float,
                 E22 : float,
                 v12 : float,
                 G12 : float) -> list:
        """
        Calculate the Rigid Reduce Transform
        """
        v21 = v12
        Q66 = G12
        
        try:
            Q11 = E11**2/(E11 - v12**2*E22)
            Q22 = E11*E22/(E11 - v12**2*E22)
            Q1221 = v12*E11*E22/(E11 - v12**2*E22)

        except:
            Q11 = E11/(1 - v12*v21)
            Q22 = E22/(1 - v12*v21)
            Q1221 = v12*E22/(1 - v12*v21)
        
        return [Q11, Q22, Q1221, Q66]
    
    def calc_MRRT(self, theta : float, index : int = 0):
        """
        Calculate the Transformed Reduced Stiffness Matrix of each laminate
        with respect to the Global Coordinate System
        
        Args:
            _float_ theta : layer's angle 
        """
        m = np.cos(theta * np.pi / 180)
        n = np.sin(theta * np.pi / 180)
        E11 = self.layers[index]['properties']['E11']
        E22 = self.layers[index]['properties']['E22']
        v12 = self.layers[index]['properties']['v12']
        G12 = self.layers[index]['properties']['G12']
        
        [Q11, Q22, Q1221, Q66] = self.calc_RRT(E11, E22, v12, G12)
        
        Q11_ = Q11*m**4 + 2*m**2*n**2*(Q1221 + 2*Q66) + Q22*n**4
        Q12_ = (Q11 + Q22 - 4*Q66)*n**2*m**2 + Q1221*(n**4 + m**4)
        Q22_ = Q11*n**4 + 2*(Q1221 + 2*Q66)*n**2*m**2 + Q22*m**4
        Q16_ = (Q11 - Q1221)*n*m**3 + (Q1221 - Q22)*n**3*m - 2*m*n*(m**2 - n**2)*Q66
        Q26_ = (Q11 - Q1221)*n**3*m + (Q1221 - Q22)*n*m**3 + 2*m*n*(m**2-n**2)*Q66
        Q66_ = (Q11 + Q22 - 2*Q1221 - 2*Q66)*n**2*m**2 + Q66*(n**4+m**4)

        Q_ =    [[Q11_, Q12_, Q16_],
                [Q12_, Q22_, Q26_],
                [Q16_, Q26_, Q66_]]
        
        self.Q_.append(Q_)
     
    def calc_matrix(self):
        """
        The function iterates over a list of angles, calculates the MRRT for each angle, and then calculates
        a matrix.
        """
        for i in range(1, len(self.h)):
            self.A += np.array(self.Q_[i-1])*(self.h[i] - self.h[i-1])
            self.B += 1/2*np.array(self.Q_[i-1])*(self.h[i]**2 - self.h[i-1]**2)
            self.D += 1/3*np.array(self.Q_[i-1])*(self.h[i]**3 - self.h[i-1]**3)
        
        #Equações Constitutivas Completamente Invertidas
        
        self.A_star = np.linalg.inv(self.A)
        self.B_star = -np.dot(self.A_star, self.B)
        self.C_star = np.dot(self.B, self.A_star)
        self.D_star = self.D - np.dot(np.dot(self.B, self.A_star), self.B)
        
        self.D_line = np.linalg.inv(self.D_star)
        self.A_line = self.A_star + np.dot(np.dot(self.B_star, self.D_line), np.transpose(self.B_star))
        self.B_line = np.dot(self.B_star, self.D_line)
        self.C_line = self.B_line
        
        #deformações em relação aos eixos globais. Mesma para todas as lâminas
        self.epsilon_global = np.dot(self.A_line, np.array(self.N) + np.dot(self.B_line, self.M))
        self.K_global = np.dot(self.C_line, np.array(self.N) + np.dot(self.D_line, self.M))
        
        self.system_properties['properties'] = {'A' : self.A,
                                                'B' : self.B,
                                                'D' : self.D,
                                                'Global deformation' : self.epsilon_global}
        
        #calculando as tensões globais na lâmina
        for i in range(1, len(self.h)):
            m = np.cos(self.angles[i-1] * np.pi / 180)
            n = np.sin(self.angles[i-1] * np.pi / 180)
            
            T = np.array([[m**2, n**2, 2*m*n],
                          [n**2, m**2, -2*m*n],
                          [-m*n, m*n, m**2 - n**2]])
            
            self.sigma_global.append(np.dot(self.Q_[i-1], self.epsilon_global + (self.h[i] - self.h[i-1])*self.K_global))
            self.sigma_local.append(np.dot(T, self.sigma_global[i-1]))
            self.epsilon_local.append(np.dot(T, self.epsilon_global))
            self.blades_data[self.layers[i]["name"]] = {'Global sigma' : self.sigma_global[i-1].tolist(),
                                                  'Local sigma' : self.sigma_local[i-1].tolist(),
                                                  'Local deformation' : self.epsilon_local[i-1]}
        
    def run(self):
        """
        Run simulation
        """
        for n in range(self.n_layers):
            self.calc_MRRT(self.angles[n], index = n)
        self.calc_matrix()
        
    def save_data(self, filename: str = 'results.json'):
        combined_data = {
            **self.system_properties,
            **self.blades_data,
        }

        # Convert NumPy arrays to Python lists
        def convert_to_json_serializable(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            raise TypeError(f"Object of type {obj.__class__.__name__} is not JSON serializable")

        # Serialize the data to JSON
        with open(filename, "w") as f:
            json.dump(combined_data, f, indent=4, default=convert_to_json_serializable)
                    
            
if __name__ == "__main__":
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