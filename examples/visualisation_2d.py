import cpm
from skimage.segmentation import mark_boundaries
import numpy as np
import pygame
import numpy as np
from pygame import surfarray
import os
from timeit import default_timer as timer
import random

def get_color(p):
    palette = [[78,78,78], [190, 30, 45], [241, 90, 41], [247, 148, 29], [254, 206, 7]]

    nr = len(palette)
    intervals = nr-1
    interval_length = 1.0/intervals
    left = (p / interval_length).astype(np.uint32)
    lerp = (p - (left* interval_length))/interval_length
    right = left + 1
    right[right > nr-1] = left[right>nr-1]

    c = np.zeros(p.shape+(3,), dtype=np.uint32)
    for i in range(intervals):
        c[left==i, 0] = palette[i][0] * (1 - lerp[left==i]) + palette[i+1][0] * lerp[left==i]
        c[left==i, 1] = palette[i][1] * (1 - lerp[left==i]) + palette[i+1][1] * lerp[left==i]
        c[left==i, 2] = palette[i][2] * (1 - lerp[left==i]) + palette[i+1][2] * lerp[left==i] 
    c[right==nr-1] = palette[-1]
    return c


def create_chemo_simulation(dimension):
    number_of_types = 2
    temperature = 10
    simulation = cpm.Cpm2d(dimension, number_of_types, temperature)

    
    simulation.set_constraints(cell_type = 1, lambda_area = 30, target_area = 500)
    simulation.set_constraints(cell_type = 1, lambda_perimeter = 2, target_perimeter = 300)
    simulation.set_constraints(cell_type = 0, other_cell_type = 1, adhesion = 5)
    simulation.set_constraints(cell_type = 1, other_cell_type = 1, adhesion = 20)
    simulation.set_constraints(cell_type = 1, max_act = 40, lambda_act=500)
    simulation.set_constraints(cell_type = 1, other_cell_type = 2, adhesion = 40)
    simulation.set_constraints(cell_type = 2, fixed=1)

    return simulation


pygame.init()

dimension = 256

display = pygame.display.set_mode((dimension, dimension))

img = np.zeros((dimension,dimension,3),dtype=np.uint8)


init_time = timer()
frames_displayed = 0



sim = create_chemo_simulation(dimension)


state = np.load("initial_state.npy")
types = state // 2**24
ids = state % 2**24
nr_of_cells = np.max(ids) + 1
print(nr_of_cells)

state = np.ones((dimension,dimension)).astype(np.uint32) * (1 + 2**24 * 2)
state[:, 128-8:128+8] = 0
sim.initialize_from_array(state, 1)

def add_cell(sim, dimension):
    lattice = sim.get_state()
    while True:
        x = random.randint(0, dimension-1)
        y = random.randint(0, 15) + 128-8
        if lattice[x,y] == 0:
            break
    print(x,y)
    sim.add_cell(1, y, x)
    
add_cell(sim, dimension)

active_cells = 1
total_cells = 1


max_act = 40



act = sim.get_act_state()
tick_count = 100

sim.run(tick_count)

exec(open("constraints_2d.py").read())
last_change = os.stat("constraints_2d.py").st_mtime

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    try:
        timestamp = os.stat("constraints_2d.py").st_mtime
    except:
        pass
    if timestamp > last_change:
        last_change = timestamp
        try:
            exec(open("constraints_2d.py").read())
        except:
            print("constraints file has error")
    if active_cells < total_cells:
        add_cell(sim, dimension)
        active_cells += 1

    sim.run(tick_count)
    state = sim.get_state()
    types = state // 2**24
    ids = state % 2**24
    img = np.zeros((dimension,dimension,3),dtype=np.uint8)
    img[ids!=0] = [78, 78, 78]

    act_vals = np.maximum(0, (act-np.max(act)+(max_act*2)))/(max_act*2)
    col = get_color(act_vals)
    img[act_vals>0] = col[act_vals>0]
    img[ids==0] = [255, 255, 255]

    img[types==2] = [150, 150, 150]

    #img = 255*mark_boundaries(img/255., ids!=0, color=(0,0,0), mode="outer")

    surfarray.blit_array(display, img)
    pygame.display.flip()

    frames_displayed+=1


print("average frame rate:", frames_displayed/(timer()-init_time), "fps")

pygame.quit()
