### recursive-raytracer

#### Compile: 
	$ make

#### Clean: 
	$ make clean

#### Execute: 
	./raytracer ./path/to/driverXX.txt ./path/to/diverXX.ppm
	
#### Format of driverXX.txt:

The first three lines are the location of the focal point (i.e., the eye point), the look at point, and the up vector. Then comes the focal length, i.e., the distance from the focal point to the image plane (near clipping plane). The 'bounds' values indicate the minimum and maximum extend of the bounded image rectangle on the infinite image plane in the camera horizontal and vertical directions respectively. Then the resolution values separately indicate the pixel sampling resolution across the horizontal and vertical dimensions of the bounded rectangle. One feature of this specification format is that you can generate intermediate cameras with low resolution, say 8 by 8 or even 4 by 4, when developing and debugging code. This speeds development considerably. The recursion depth is specified by the field "recursionLevel 3".

After resolution is specified, the next line specifies the ambient illumination in the scene. This is a low level 'white' light with values of 0.1 for red, green, and blue bands on a scale of zero to one. After specifying the amount of ambient light in the scene, zero or more light sources may be specified. The first four values given are the x, y, z, and w coordinates of the light source in world coordinates. The fourth value w is generally one, but a zero indicates a light source at infinity in the direction specified by x, y, and z. The last three values indicate the red, green, and blue levels of the light source on a zero to one scale.

Finally, zero or more polygonal models can be specified for inclusion in the scene. The format of polygonal model inclusion is "model wx wy wz theta scale tx ty tz model.obj". The keyword model identifies that the line is a description model and its transform. The triplet wx wy wz represents the axis about which to rotate followed by theta in degrees representing the angle by which to rotate (thus Axis-Angle format). Next, the scale is a uniform scaling factor to shrink or grow a model. Next, the triplet tx ty tz defines a model-to-world translation. Finally, model.obj is the model file on which you will be applying the defined 3D transformation.

Here is a very important detail not to overlook. The *.obj file for a model specifies a material file that your program must load to understand the material properties associated with the faces in the model. A material file has a standard format with .mtl extension.

There could be zero or more spheres. The first three values are the x, y, and z coordinates of the center of the sphere in world coordinates. The fourth value is the radius of the sphere. All the next values are the simplified material properties indicating the 'color' of the sphere in terms of red, green, and blue. First triplet after the radius represents (Ka_red, Ka_green, Ka_blue) ambient coefficients. Next triplet represents diffuse (Kd_red, Kd_green, Kd_blue) coefficients. Next we have (Ks_red, Ks_green, Ks_blue) specular coefficients. Finally, we have (Kr_red, Kr_green, Kr_blue) attenuation coefficients. For spheres the phong constant is taken as 16.


#### This program only supports outputing .ppm images. 

#### [Read more here](http://www.cs.colostate.edu/~cs410/yr2018fa/home_assignments.php)
