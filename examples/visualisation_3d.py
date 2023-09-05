import cpm
import pygame
import math
from OpenGL.GL import *
from OpenGL.GLU import *
import os
from skimage.measure import marching_cubes
import numpy as np
import random


import scipy.ndimage as nd

from OpenGL.GL import shaders

active_cells = 0
simulated_cells = 0

def init_opengl(width, height):
    glViewport(0, 0, width, height)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(90.0, float(width)/height, 0.1, 1000.0)
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    #glEnable(GL_DEPTH_TEST)
    glShadeModel(GL_SMOOTH)
    glClearColor(0.05, 0.0, 0.1, 0.0)
    glPointSize(3)

    vertexBuffer = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer)
    normalBuffer = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer)
    colorBuffer = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer)
    indexBuffer = glGenBuffers(1)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer)

    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)
    glEnableClientState( GL_COLOR_ARRAY )

    glEnable(GL_NORMALIZE)
    glLightfv(GL_LIGHT0,GL_POSITION,[ .0, 10.0, 10., 0. ] )
    glLightfv(GL_LIGHT0,GL_AMBIENT,[ .5, .5, .5, 1.0 ]);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,[ 0.85, 0.85, 0.85, 1.0 ]);
    glLightfv(GL_LIGHT0,GL_SPECULAR,[ 0.5, 0.5, 0.5, 1.0 ]);


    VERTEX_SHADER = shaders.compileShader("""#version 120
            varying vec3 N;
            varying vec3 v;
    
            void main() {
                //gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;

                
               v = vec3(gl_ModelViewMatrix * gl_Vertex);       
               N = normalize(gl_NormalMatrix * gl_Normal);

               gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;


            }""", GL_VERTEX_SHADER)
    
    FRAGMENT_SHADER = shaders.compileShader("""#version 120
            varying vec3 N;
            varying vec3 v;
    
            void main() {
               vec3 L = normalize(gl_LightSource[0].position.xyz - v);   
               vec3 E = normalize(-v); // we are in Eye Coordinates, so EyePos is (0,0,0)  
               vec3 R = normalize(-reflect(L,N));  
             
               //calculate Ambient Term:  
               vec4 Iamb = gl_FrontLightProduct[0].ambient;    

               //calculate Diffuse Term:  
               vec4 Idiff = gl_FrontLightProduct[0].diffuse * max(dot(N,L), 0.0);
               Idiff = clamp(Idiff, 0.0, 1.0);     
               
               // calculate Specular Term:
               vec4 Ispec = vec4(0.3,0.3,0.3,1.0) 
                            * pow(max(dot(R,E),0.0),40);
               Ispec = clamp(Ispec, 0.0, 1.0); 

               // write Total Color:  
               gl_FragColor = gl_FrontLightModelProduct.sceneColor + Iamb + Idiff + Ispec;     

        }""", GL_FRAGMENT_SHADER)

    shader = shaders.compileProgram(VERTEX_SHADER,FRAGMENT_SHADER)
    shaders.glUseProgram(shader)
    
    #glEnable(GL_LIGHT0)
    #glEnable(GL_DEPTH_TEST)
    

    return vertexBuffer, normalBuffer, indexBuffer, colorBuffer
    


def make_sim(dimension, maxact):
    number_of_types = 3
    temperature = 20

    sim = cpm.Cpm3d(dimension, number_of_types, temperature)

    return sim

def add_cell(sim, lattice, celltype, dimension):
    while True:
        x = random.randint(0, dimension-1)
        y = random.randint(0, dimension-1)
        z = random.randint(0, dimension-1)
        if lattice[x,y,z] == 0:
            break
    print(x,y,z)
    print(dimension, celltype)
    sim.add_cell(celltype, x,y,z)



def display_centroids(centroids, cellCount, path_history, max_path, dimension, sphere):
    glColor3f(0.3,0.3,0.4)

    for i, centroid in enumerate(centroids):
        if len(path_history) < i + 1:
            path_history.append([])
        path_history[i].append(centroid)

        if len(path_history[i]) > max_path:
            path_history[i] = path_history[i][-max_path:]

    #glBegin(GL_POINTS)
    a = 0
    for centroid in centroids:

        glPushMatrix()
        glTranslatef(centroid[2], centroid[1], centroid[0]) #Move to the place
        gluSphere(sphere, 0.5, 32, 16) #Draw sphere
        glPopMatrix()
        

        #glVertex3f(centroid[2], centroid[1], centroid[0])
    #glEnd()

    for path in path_history:
        glBegin(GL_LINE_STRIP)
        prev = None
        for v in path:
            if not prev is None:
                dx = abs(prev[0] - v[0])
                dy = abs(prev[1] - v[1])
                dz = abs(prev[2] - v[2])
                if dx > dimension/2 or dy > dimension/2 or dz > dimension/2:
                    glEnd()
                    glBegin(GL_LINE_STRIP)
            prev = v
            glVertex3f(v[2], v[1], v[0])
        glEnd()

def display_lattice(lattice, dimension, vertexBuffer, normalBuffer, indexBuffer, colorBuffer, act, maxact, cellIndex):
    colorindex = cellIndex % 3
    if colorindex == 0:
        glColor3f(0.3,0.3,0.4)
    elif colorindex == 1:
        glColor3f(0.7,0.3,0.4)
    elif colorindex == 2:
        glColor3f(0.3,0.6,0.4)
    lattice = np.copy(lattice)
    lattice[0] = 0
    lattice[-1] = 0
    lattice[:,0] = 0
    lattice[:,-1] = 0
    lattice[:,:,0] = 0
    lattice[:,:,-1] = 0

    #lattice = np.copy(lattice)
    lattice = lattice % 2**24

    if not np.any(lattice == (cellIndex+1)):
        return

    verts, faces, normals, values = marching_cubes((lattice==(cellIndex+1)).astype(np.float32), 0.5)

    xs = np.array([i[0] for i in verts]).astype(np.uint32)
    ys = np.array([i[1] for i in verts]).astype(np.uint32)
    zs = np.array([i[2] for i in verts]).astype(np.uint32)
    act_vals = act[xs,ys,zs].astype(np.float32)

    colors = np.ones((len(verts), 3), dtype=np.float32) *0.5

    colors[:,2] =  0.6
    #colors[:,0] =  np.maximum(maxact + act_vals -np.max(act_vals) , 0.0) / (maxact*2) + 0.5
    colors[:,0] =  0.5
    #colors[cellids % 2 == 0] = [1,1,1]

    #print(len(colors.ravel()))
    #print(np.array(verts).shape)
    #print(colors.shape)

    glBindBuffer(GL_ARRAY_BUFFER, colorBuffer)
    glBufferData(GL_ARRAY_BUFFER, colors, GL_STATIC_DRAW)
    glColorPointer(3, GL_FLOAT, 0, None)


    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer)
    glBufferData(GL_ARRAY_BUFFER, np.array(verts).astype(np.float32), GL_STATIC_DRAW)
    glVertexPointer(3, GL_FLOAT, 0, None)
    glBindBuffer(GL_ARRAY_BUFFER, normalBuffer)
    glBufferData(GL_ARRAY_BUFFER, (np.array(normals)*1).astype(np.float32), GL_STATIC_DRAW)
    glNormalPointer(GL_FLOAT, 0, None)
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer)
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, np.array(faces).astype(np.uint32), GL_STATIC_DRAW)

    glEnable(GL_LIGHTING)
    glEnable(GL_COLOR_MATERIAL)
    #glColor3f(1.0, 0.0,0.0)
    glDrawElements(GL_TRIANGLES, len(faces)*3, GL_UNSIGNED_INT, None)
    


if __name__ == "__main__":
    pygame.init()
    size = width, height = 800, 600

    maxact = 15

    screen = pygame.display.set_mode(size, pygame.DOUBLEBUF | pygame.OPENGL)
    vertexBuffer, normalBuffer, indexBuffer, colorBuffer = init_opengl(width,height)

    sphere = gluNewQuadric()

    font = pygame.font.Font(None, 32) 
    text = font.render('GeeksForGeeks', True, (1,0,0), (0,1,0)) 
    textRect = text.get_rect()  
    textRect.center = (width // 2, height // 2) 

    dimension = 32
    cells_visualised = 1000

    sim = make_sim(dimension, maxact)
    state = sim.get_state()
    act = sim.get_act_state()



    xoff = 0
    yoff = 0
    zoff = -40
    runrate = 1
    d = 0
    e = 0
    max_path = 30

    path_history = []

    exec(open("constraints_3d.py").read())
    last_change = os.stat("constraints_3d.py").st_mtime

    done = False

    pressed = False
    while not done:
        try:
            timestamp = os.stat("constraints_3d.py").st_mtime
        except:
            pass
        if timestamp > last_change:
            last_change = timestamp
            try:
                exec(open("constraints_3d.py").read())
            except:
                print("constraints file has error")

        if active_cells < simulated_cells:
            add_cell(sim, state, 1, dimension)
            active_cells += 1
        if runrate > 0:
            sim.run(runrate)
        centroids = sim.get_centroids()

        a = nd.center_of_mass(state>0)

        density = np.sum(state>0)/(dimension**3)
        print("density: ", density)


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        #d+=0.2
        glTranslatef(xoff,yoff,zoff)
        glRotate(e,1,0,0)
        glRotate(d,0,1,0)
        glTranslatef(-dimension//2, -dimension//2, -dimension//2)

        glEnable(GL_LIGHT0)
        glEnable(GL_DEPTH_TEST)
        for i in range(cells_visualised):
            display_lattice(state, dimension, vertexBuffer, normalBuffer, indexBuffer, colorBuffer, act, maxact, i)
        display_centroids(centroids, cells_visualised, path_history, max_path, dimension, sphere)


        pygame.display.flip()
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT: sys.exit()
            elif event.type == pygame.MOUSEWHEEL:
                zoff += event.y
            elif event.type == pygame.MOUSEBUTTONDOWN:
                pressed = True
            elif event.type == pygame.MOUSEBUTTONUP:
                pressed = False
            elif event.type == pygame.MOUSEMOTION and pressed:
                d += event.rel[0]
                e += event.rel[1]
                if e > 90:
                    e = 90
                if e < -90:
                    e = -90
        
    

