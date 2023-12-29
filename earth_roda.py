from vpython import *

scene = canvas(width=1280, height=720)

earth = sphere(pos = vector(0,0,0), radius=1, texture = textures.earth)

rotation_speed = 0.02

while True:
    rate(120) 

    earth.rotate(angle=rotation_speed,axis=vector(0,1,0)) 