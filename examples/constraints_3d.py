maxact = 10

sim.set_constraints(cell_type = 1, lambda_area = 25, target_area = 150)
sim.set_constraints(cell_type = 1, lambda_perimeter = 0.2, target_perimeter = 1400)
sim.set_constraints(cell_type = 1, other_cell_type = 1, adhesion = 10)
sim.set_constraints(cell_type = 0, other_cell_type = 1, adhesion = 5)
#sim.set_constraints(cell_type = 1, lambda_act = 0, max_act = maxact)
sim.set_constraints(cell_type = 1, lambda_persistence = 100.0, persistence_diffusion = 0.5, persistence_time = 20)
sim.set_constraints(cell_type = 1, fixed=0)


#visualisation offset
xoff = 0
yoff = 0
zoff = -37

#mcs/frame
runrate = 2

#how many cells do we render fully (instead of a sphere)
cells_visualised = 10

#number of active cells of type 1
simulated_cells = 100


