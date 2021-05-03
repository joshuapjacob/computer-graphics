# Computer Graphics
This repository contains my coursework for CSE306 "Computer Graphics" @ Ecole Polytechnique.

## Simple Ray-Tracer (from scratch)

![Sample Render](https://raw.githubusercontent.com/joshuapjacob/Computer-Graphics/main/renders/final.png)
### Features

- Spheres
- Reflection
- Refraction (Solid & Hollow)
- Indirect Lighting
- Spherical Lighting (Soft Shadows)
- Parallelization
- Antialiasing
- Triangle Meshes (.obj)
- Bounding Volume Hierarchies Acceleration

### Compilation

```g++ -O3  main.cpp vector.cpp -fopenmp -lpthread -std=c++11 -o main.out```
