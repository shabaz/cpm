max_act = 40
sim.set_constraints(cell_type = 1, lambda_area = 30, target_area = 500)
sim.set_constraints(cell_type = 1, lambda_perimeter = 2, target_perimeter = 300)
sim.set_constraints(cell_type = 0, other_cell_type = 1, adhesion = 5)
sim.set_constraints(cell_type = 1, other_cell_type = 1, adhesion = 20)
sim.set_constraints(cell_type = 1, max_act = max_act, lambda_act=150)
sim.set_constraints(cell_type = 1, other_cell_type = 2, adhesion = 40)
tick_count = 10
total_cells = 4


