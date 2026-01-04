# Mathematics Behind the Program

[AI Conversation](https://chatgpt.com/share/695a4393-9898-800d-9277-c0dcd8628821)

## Mathematical ideas used
This program is based on basic coordinate geometry and vector mathematics.
The main mathematical ideas are:
### Distance between two points
Distance is computed using the Pythagorean theorem.
It measures how far two points are from each other in 2D or 3D space.
### Vectors
Subtracting two points produces a vector.
Vectors are used to represent directions, lengths, and perpendicular relationships.
### Dot product
The dot product is used to check whether two directions are perpendicular and to compute projections.
It is also used in calculating distances from a point to a line or a plane.
### Geometric object definitions
    - A line is represented by one point and one direction vector.
    - A circle/sphere is defined as all points that have the same distance from a center point.
    - A plane is defined by a point and a normal vector.

## Why the algorithms work
Each algorithm directly follows the mathematical definition of the geometric object.
### Distance algorithms follow the Euclidean distance formula derived from the Pythagorean theorem.
### Point–line and point–plane distances are calculated using vector projections, which always give the shortest (perpendicular) distance.
### Intersection algorithms work by solving two geometric definitions at the same time (for example, solving both the line equation and the circle or plane equation together).
Because the algorithms are direct implementations of geometric definitions, the computed results correctly represent the real geometric meanings.

## Connection between geometry theory and the code
Geometry theory defines what points, lines, circles, and planes are.
The code turns these definitions into computational steps.
### Points are stored as coordinates.
### Vectors are obtained from coordinate differences.
### Geometric objects (lines, circles, planes) are implemented using their mathematical definitions.
This makes the program a computational version of classical geometry.
