#Raycaster with Lighting


This program incorporates the simple lighting model into a basic ray caster. The ray caster reads in a json scene file and outputs a .ppm image of the objects in the scene. It includes support for point and spotlights, and uses specular and diffuse color components for shading of primative objects


##Usage


To build the program, simply type ‘make’ and the project will build


Use the format ./raycast width height object-file.json output-file.ppm


width and height cannot be zero


Please make sure the object-file is a valid .json file. See the example json file for a format example
