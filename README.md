# Computer Graphics
This repository contains my coursework for CSE306 "Computer Graphics" @ Ecole Polytechnique.

## Ray Tracer

![Sample Render](https://raw.githubusercontent.com/joshuapjacob/computer-graphics/main/Ray%20Tracer/renders/final.png)
### Features

- Spheres
- Reflection
- Refraction (Solid & Hollow)
- Indirect Lighting
- Spherical Lighting (Soft Shadows)
- Basic Parallelization
- Antialiasing
- Triangle Meshes (.obj)
- Bounding Volume Hierarchies Acceleration

### Requirements
- [stb](https://github.com/nothings/stb) (to write image)

### Compilation &
```g++ -O3 main.cpp vector.cpp -fopenmp -lpthread -std=c++11 -o main.out```

## Sliced Optimal Transport Color Matching

Input | Model | Output
:----:|:-----:|:------:
![](https://...Dark.png)  |  ![](https://...Ocean.png) |

### Requirements
- [stb](https://github.com/nothings/stb) (to read/write image)
### Compilation & Execution
Last argument is the number of iterations.
```
g++ -O3 main.cpp vector.cpp -std=c++11 -o main.out
./main.out input.png model.png 100
```
