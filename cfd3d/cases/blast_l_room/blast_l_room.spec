# L-shaped room with one obstacle and one vent.
# Bounding box: 1 m^3, 30^3 cells (dx = dy = dz ≈ 0.0333 m).
# Polygonal floor: 1 m x 1 m square with the [0.5,1] x [0.5,1] corner
# removed, leaving an L (bottom arm runs along +x at y<0.5, top arm
# runs along +y at x<0.5).
# Vent: 0.2 m x 0.2 m on the x_max face, centred in the bottom arm.
# Obstacle: a 0.15 m cube on the floor of the bottom arm near the bend.

GRID    30 30 30
BBOX    1.0 1.0 1.0
ROOM_Z  0.0 1.0

POLYGON 6
VERT  0.0  0.0
VERT  1.0  0.0
VERT  1.0  0.5
VERT  0.5  0.5
VERT  0.5  1.0
VERT  0.0  1.0

VENT    1   0.15 0.40   0.35 0.60

CUBE    0.60 0.15 0.40   0.75 0.30 0.55

END
