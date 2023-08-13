![teaser](images/teaser-s.png)

# Lightweight Curvature Estimation on Point Clouds with Randomized Corrected Curvature Measures

[![CI (linux/macOS)](https://github.com/JacquesOlivierLachaud/PointCloudCurvCNC/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/JacquesOlivierLachaud/PointCloudCurvCNC/actions/workflows/cmake-multi-platform.yml)

This repository provides an example application that implements the curvature estimators on oriented point clouds presented at [Symposium on Geometry Processing 2023, Genova, Italy, July 3-7](https://sgp2023.github.io):

Jacques-Olivier Lachaud, David Coeurjolly, Céline Labart, Pascal Romon, Boris Thibert, **Lightweight Curvature Estimation on Point Clouds with Randomized Corrected Curvature Measures**, *Comput. Graph. Forum*, 42(5), 2023.

**Authors of the code**:

* Jacques-Olivier Lachaud
* David Coeurjolly

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

### Computing curvatures of point clouds approximating simple shapes

```
./curvatures
```

Run the program with a GUI that allows you to generate point clouds approximating simple shapes (sphere, torus, cube, dodecahedron).

### Computing curvatures of an oriented point cloud (given as text file)

```
./curvatures ../data/bearded-man-xyz-nxyz.pts
```

Text file should be composed of lines of the form `x y z nx ny nz` for each point, determining the coordinates (x,y,z) of each point and the components (nx,ny,nz) of its oriented normal vector.

### Computing curvatures of an oriented point cloud (given as a wavefront OBJ file)

```
./curvatures ../data/bunnyhead.obj
```

The OBJ file should contained the vertices as `v x y z` and normal vectors as `vn nx ny nz`.

### Using the interface

* Shape generation
  - `Sphere` generates a sphere of size `R` with `N` points
  - `Torus` generates a torus of great radius `R` and small radius `r` with `N` points
  - `Cube` generates a cube of size `R` with `N` points
  - `Dodecahedron` generates a dodecahedron of size `R` with `N` points
  - `InputFile` generates a set of points from the given input file
  - all shapes can be perturbated in position (parameter `x`) and in normals (parameter `xi`)
  - you may toggle between fast display (better when N > 5 millions) or nice display

* Curvature computation
  - `Curvatures` computes all curvature information (mean, Gaussian, principal curvatures and directions.
  - parameter `K` is the chosen number of nearest neighbors
  - parameter `L` is the chosen number of triangles (only for `Uniform` generation)
  - you may choose for triangle random generation methods among `Uniform`, `Independent`, `Hexagram`, `Avg-Hexagram` (see paper for details)
  - method `Avg-Hexagram` is the fastest, very accurate while staying robust to noise
  - parameter `W` balances between the normal to the points and the average normals of its neighbors to define the local sampling plane, 0.5 works well.

* Information
  - the total computation time of curvatures is displayed
  - l2 and loo errors are displayed for `Sphere` and `Torus` shape.
  
# Replicate paper results

We provide additional datasets to be able to reproduce the paper results. To exactly get the paper visual after clicking on `Curvatures` the button, for each interested curvature quantity, change polyscope's colormap from `viridis` to  `coolwarm` and use the clamping values given in the following table (to specify the colormap min/max values, just CTRL-click to the min/max text area for the selected curvature quantities).


| Pointcloud | Paper figure | Params | curvature colormap clamping ranges |
| -----------| -------------| -------| ---------------- |
| ` ./curvatures ../data/goursat/goursat_025000.pts` | Fig6 |K=50, L=100 | Mean: [-0.107, 0.345], Gaussian: [-0.034, 0.119] |
| ` ./curvatures ../data/torus/torus_025000.pts` | supplem Fig2 |K=50, L=100 | Mean: [ 0.125, 0.32], Gaussian: [-0.125, 0.0625] |
| ` ./curvatures ../data/goursat/goursat_025000.pts` | supplem Fig3 |K=50, L=200 | Mean: [-0.107, 0.345], Gaussian: [-0.034, 0.119] |
| ` ./curvatures ../data/torus/torus_025000.pts` | supplem Fig4 |K=50, L=200 | Mean: [ 0.125, 0.32], Gaussian: [-0.125, 0.0625] |
| ` ./curvatures ../data/anisotropic/goursat-filter.pts` | supplem Fig5 |K=50, L=200 | Mean: [-0.107, 0.345], Gaussian: [-0.034, 0.119] |
| ` ./curvatures ../data/anisotropic/torus-filter.pts` | supplem Fig6 |K=50, L=200 | Mean: [ 0.125, 0.32], Gaussian: [-0.125, 0.0625] |
| ` ./curvatures ../data/lidarSim/goursat-lidar.pts` | supplem Fig5 |K=50, L=200 | Mean: [-0.107, 0.345], Gaussian: [-0.034, 0.119] |
| ` ./curvatures ../data/lidarSim/torus-lidar.pts` | supplem Fig6 |K=50, L=200 | Mean: [ 0.125, 0.32], Gaussian: [-0.125, 0.0625] |


## Examples

| ![bearded man scan](images/bearded-H-AvgHexagram-20-s.png) |
| --------- |
| Mean curvature (Avg-Hexagram, N=20) |

|![torus H](images/torus-H-AvgHexagram-20-s.png) | ![torus V1](images/torus-V1-AvgHexagram-20-s.png) | ![torus V2](images/torus-V2-AvgHexagram-20-s.png) |
| ---- | ---- | ---- |
| Mean curvature on torus (N=10000, N=20) | First principal direction | Second principal direction |



