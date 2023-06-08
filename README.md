# Lightweight Curvature Estimation on Point Clouds with Randomized Corrected Curvature Measures

This repository provides an example application that implements the curvature estimators on oriented point clouds presented in the paper presented at [Symposium on Geometry Processing 2023, Genova, Italy, July 3-7](https://sgp2023.github.io):

Jacques-Olivier Lachaud, David Coeurjolly, Céline Labart, Pascal Romon, Boris Thibert, **Lightweight Curvature Estimation on Point Clouds with Randomized Corrected Curvature Measures**, to appear in *Comput. Graph. Forum*, 2023.

## Building

Once cloned, proceed as follows on Linux/macos:

```
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j 8
```

It will automatically fetch and install the dependencies [polyscope](https://polyscope.run) and [eigen](https://eigen.tuxfamily.org/).

## Usage

* Computing curvatures of point clouds approximating simple shapes

```
./curvatures
```

Run the program with a GUI that allows you to generate point clouds approximating simple shapes (sphere, torus, cube, dodecahedron).

* Computing curvatures of an oriented point cloud (given as text file)

```
./curvatures ../data/bearded-man-xyz-nxyz.pts
```

Text file should be composed of lines of the form `x y z nx ny nz` for each point, determining the coordinates (x,y,z) of each point and the components (nx,ny,nz) of its oriented normal vector.

* Computing curvatures of an oriented point cloud (given as a wavefront OBJ file)

```
./curvatures ../data/bunnyhead.obj
```

The OBJ file should contained the vertices as `v x y z` and normal vectors as `vn nx ny nz`.


## Examples



