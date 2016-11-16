#Raycaster with Lighting


This program builds upon the previous ray casting projects to incorporate ray tracing to
allow for the modeling of refraction and reflection. By modifying the behavior of he ray caster
to allow for recursive calls to the trace method, we can represent reflection and refraction by
tracing additional rays that continue past the first collision with modifications to their direction
depending upon the amount and type of reflection or refraction modeled by each object


##Usage


To build the program, simply type ‘make’ and the project will build


Use the format ./raytrace width height object-file.json output-file.ppm


width and height cannot be zero


Please make sure the object-file is a valid .json file. See the example json file for a format example
