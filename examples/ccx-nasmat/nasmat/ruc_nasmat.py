"""
Solve RVE using NASMAT: https://software.nasa.gov/software/LEW-20244-1
"""
import numpy as np
from ctypes import cdll, byref, c_double, POINTER
import time

from pynasmat import Model, Constituent, Ruc
from pynasmat.constants import MAT_TRANSVERSE_ISOTROPIC, RUC_ARCHID_SQPACK_1FIBER, RUC_MODID_GMC_2D

class MicroSimulation:

    def __init__(self, sim_id):
        """
        Constructor of MicroSimulation class.
        """
        self._dims = 3
        self._sim_id = sim_id

        print("MicroSimulation object created with ID: ", sim_id)

    def solve(self, macro_data, dt):
        assert dt != 0

        print(macro_data)

        # Processing input data
        rve_id = int(macro_data["input_id"])
        mod_id = 102
        ruc_size = 16
        input_strain = np.zeros((6))
        for i in range(3):
            input_strain[i] = float(macro_data["strain1to3"][i])
            input_strain[i + 3] = float(macro_data["strain4to6"][i])

        # Building NASMAT Model
        model = self.build_nasmat_model(rve_id, mod_id, ruc_size)

        # Solving RVE using NASMAT
        cmat = model.homogenize(print_output=0)
        print("Homogenization completed for RVE ID: ", self._sim_id)
    
        print(cmat)


        # Stress calculation
        stresses = np.zeros((6))
        stresses = np.dot(cmat, input_strain)


        return {"stress1to3": stresses[0:3], "stress4to6": stresses[3:6],
                "cmat1": cmat[0][0:3], "cmat2": cmat[0][3:6], 
                "cmat3": cmat[1][1:4], "cmat4": np.array(list(cmat[1][4:6])+list(cmat[2][2:3])),
                "cmat5": cmat[2][3:6], "cmat6": cmat[3][3:6],
                "cmat7": np.array(list(cmat[4][4:6])+list(cmat[5][5:6])), 
                "conv_flag": 1}

    def get_state(self):
        # ID is returned as it is trivial. In real case, this method should return the state of the simulation.
        return self._sim_id

    def set_state(self, state):
        # ID is set as it is trivial. In real case, this method should set the state of the simulation.
        self._sim_id = state

    def build_nasmat_model(self, rve_id, mod_id, ruc_size):

        # Building NASMAT Model
        model = Model(name='example_01', id=rve_id) 

        # Adding constituents
        const1 = Constituent(name='Fiber',
                     model=MAT_TRANSVERSE_ISOTROPIC, 
                     id=1, 
                     data=[388.2E3,  7.6E3,0.41,0.45,  14.9E3,-0.68E-6,9.74E-6], # Units in MPa 
                     comments='Fiber')
        model.add_constituent(const1)
        const2 = Constituent(name='Matrix',
                     model=MAT_TRANSVERSE_ISOTROPIC, 
                     id=2, 
                     data=[388.2E3,  7.6E3,0.41,0.45,  14.9E3,-0.68E-6,9.74E-6], # Units in MPa 
                     comments='Matrix')
        model.add_constituent(const2)

        # Adding RUC
        if (mod_id%10 == 2): # 2D Model
            ruc1 = Ruc( arch_id=RUC_ARCHID_SQPACK_1FIBER, 
                        mod_id=mod_id,
                        nb=ruc_size, ng=ruc_size, 
                        fiber_id=1, matrix_id=2, vf=0.6, 
                        lh=1.0, ll=1.0)
        elif (mod_id%10 == 3): # #D Model
            ruc1 = Ruc(arch_id=RUC_ARCHID_SQPACK_1FIBER, 
                        mod_id=mod_id,
                        na = ruc_size, nb=ruc_size, ng=ruc_size,  
                        fiber_id=1, matrix_id=2, vf=0.6, 
                        ld = 1.0, lh=1.0, ll=1.0)
        model.add_ruc(ruc1)
        return model

def main():
    for i in range(1, 3):
        sim = MicroSimulation(i)
        strains = { "rve_id": i,
                    "mod_id": RUC_MODID_GMC_2D,
                    "ruc_size": 32,
                    "strains1to3": [0.1, 0.2, 0.3], 
                    "strains4to6": [0.4, 0.5, 0.6]}
        dt = 0.1
        print(sim.solve(strains, dt))

if __name__ == "__main__":
    main()
