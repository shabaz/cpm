import cpm 
import random
import numpy as np

def add_cell(sim, lattice, celltype, dimension):
    while True:
        x = random.randint(0, dimension-1)
        y = random.randint(0, dimension-1)
        z = random.randint(0, dimension-1)
        if lattice[x,y,z] == 0:
            break
    sim.add_cell(celltype, x,y,z)

def make_sim(dimension):
    number_of_types = 2
    temperature = 7

    sim = cpm.Cpm3d(dimension, number_of_types, temperature)

    sim.set_constraints(cell_type = 1, lambda_area = 25, target_area = 1800)
    sim.set_constraints(cell_type = 1, lambda_perimeter = 0.1, target_perimeter = 8600)
    sim.set_constraints(cell_type = 1, other_cell_type = 1, adhesion = 0)
    sim.set_constraints(cell_type = 0, other_cell_type = 1, adhesion = 5)

    # act constraint should still work
    sim.set_constraints(cell_type = 1, lambda_act = 55, max_act = 30)


    #persistence constraint
    #sim.set_constraints(cell_type = 1, lambda_persistence = 200.0, persistence_diffusion = 0.99, persistence_time = 20)

    lattice = sim.get_state()

    for i in range(100):
        add_cell(sim, lattice, 1, dimension)

    #make a cell type 'fixed', i.e. it won't move, so can be used for e.g. FRC
    #sim.set_constraints(cell_type = 1, fixed=1)

    
    #sim.add_cell(dimension//2,dimension//2,dimension//2,1)
    return sim


dimension = 512

sim = make_sim(dimension)
state = sim.get_state()
act = sim.get_act_state()

total_mcs = 5
mcs_per_celltrack_snapshot = 205

#for i in range(total_mcs//mcs_per_celltrack_snapshot):
sim.run(total_mcs)

centroids = sim.get_centroids()
print(centroids)

#density = np.sum(state>0)/(dimension**3)
#print("density: ", density)
