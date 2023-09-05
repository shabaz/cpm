from distutils.core import setup, Extension
import numpy as np

module1 = Extension('cpm',
                    sources = ['src/python_wrapper.cpp', 'src/cpm.cpp', 'src/lattice_2d.cpp',
                        'src/lattice_3d.cpp', 'src/hamiltonian.cpp', 'src/simulation.cpp',
                        'src/dice_set.cpp', 'src/ranxoshi256.cpp', 'src/cell_states.cpp',
                        'src/linalg.cpp', 'src/centroids.cpp'],
                    include_dirs = [np.get_include(),'src'],
                    extra_compile_args=['-std=c++17', '-O3'], )

if __name__ == "__main__": setup( 
	name="speedy-cpm",
	py_modules=['gpucpm'],
	ext_modules = [module1] )

