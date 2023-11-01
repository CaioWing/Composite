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
    
    def __init__(self, t : float = 3,
                 angles = [45, 0, 0, 45], 
                 normal : list = [0, 0, 0], 
                 momentum : list = [0, 0, 0],
                 **kwargs) -> None:
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
        self.n_layers = len(angles)
        self.N = normal
        self.M = momentum
        self.epsilon_local = []
        self.h = [(-(self.n_layers//2)+i)*t/1000 for i in range(self.n_layers + 1)]
        self.angles = angles
        self.Q_ = []
        
        for key, value in kwargs.items():
            try:
                if type(value) == dict:
                    value = value['value']
                setattr(self, key, value)
            except:
                continue
                
    def calc_RRT(self):
        """
        Calculate the Rigid Reduce Transform
        """
        Q66 = self.G12
        
        try:
            Q11 = self.E11**2/(self.E11 - self.v12**2*self.E22)
            Q22 = self.E11*self.E22/(self.E11 - self.v12**2*self.E22)
            Q1221 = self.v12*self.E11*self.E22/(self.E11 - self.v12**2*self.E22)

        except:
            Q11 = self.E11/(1 - self.v12*self.v21)
            Q22 = self.E22/(1 - self.v12*self.v21)
            Q1221 = self.v12*self.E22/(1 - self.v12*self.v21)
        
        return [Q11, Q22, Q1221, Q66]
    
    def calc_MRRT(self, theta : float):
        """
        Calculate the Transformed Reduced Stiffness Matrix of each laminate
        with respect to the Global Coordinate System
        
        Args:
            _float_ theta : layer's angle 
        """
        m = np.cos(theta * np.pi / 180)
        n = np.sin(theta * np.pi / 180)
        
        [Q11, Q22, Q1221, Q66] = self.calc_RRT()
        
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
                                                'A_line' : self.A_line,
                                                'B_line' : self.B_line,
                                                'C_line' : self.C_line,
                                                'D_line' : self.D_line,
                                                'A_star' : self.A_star,
                                                'B_star' : self.B_star,
                                                'C_star' : self.C_star,
                                                'D_star' : self.D_star,
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
            
            self.blades_data['blade-' + str(i)] = {'Global sigma' : self.sigma_global[i-1].tolist(),
                                                  'Local sigma' : self.sigma_local[i-1].tolist(),
                                                  'Local deformation' : self.epsilon_local[i-1]}
        
    def run(self):
        """
        Run simulation
        """
        n_blades = len(self.h) - 1
        for n in range(n_blades):
            self.calc_MRRT(self.angles[n])
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