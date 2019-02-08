CC=g++
INCLUDE=-O2 -c -I -wall

all: raytracer

raytracer: raytracer.o vector3d.o matrix2d.o ray.o camera.o light.o plyobject.o model.o sphere.o point.o cameraFileParser.o sceneFileParser.o
	$(CC) raytracer.o vector3d.h matrix2d.h ray.h camera.h light.h plyobject.h model.h sphere.h point.h cameraFileParser.h sceneFileParser.h -o raytracer

raytracer.o: raytracer.cpp vector3d.h matrix2d.h ray.h camera.h light.h plyobject.h model.h sphere.h point.h cameraFileParser.cpp sceneFileParser.cpp
	$(CC) $(INCLUDE) raytracer.cpp

vector3d.o: vector3d.h
	$(CC) $(INCLUDE) vector3d.h

matrix2d.o: matrix2d.h vector3d.h
	$(CC) $(INCLUDE) matrix2d.h

ray.o: ray.h vector3d.h
	$(CC) $(INCLUDE) ray.h

camera.o: camera.h vector3d.h matrix2d.h
	$(CC) $(INCLUDE) camera.h

light.o: light.h vector3d.h
	$(CC) $(INCLUDE) light.h

plyobject.o: plyobject.h vector3d.h matrix2d.h
	$(CC) $(INCLUDE) plyobject.h

model.o: model.h vector3d.h matrix2d.h
	$(CC) $(INCLUDE) model.h

sphere.o: sphere.h vector3d.h
	$(CC) $(INCLUDE) sphere.h

point.o: point.h vector3d.h
	$(CC) $(INCLUDE) point.h

cameraFileParser.o: cameraFileParser.h cameraFileParser.cpp camera.h
	$(CC) $(INCLUDE) cameraFileParser.cpp

sceneFileParser.o: sceneFileParser.h sceneFileParser.cpp light.h model.h
	$(CC) $(INCLUDE) sceneFileParser.cpp

clean:
	rm -f *.o raytracer *.gch

tar:
	tar -cvf raytracer.tar *.cpp *.h *.txt *.mtl *.obj *.ppm Makefile
