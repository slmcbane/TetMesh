# TetMesh
This repository exists to provide mesh facilities for finite element simulations
in two and three dimensions on meshes composed of triangle and tetrahedra. The
mesh data structure is not explicitly tied to the FEM, but provides information
that should be useful in that context. It was primarily created to support my
research in topology optimization of lattice-type structures based on
substructuring.

Currently, the library interface is contained in `TetMesh.hpp`, which implements
a two dimensional mesh of triangles with additional boundary information.
Constructing and inspecting the mesh are implemented, with some simple sanity
checks. The primary missing functionality in this version is serialization;
I would like to be able to save the mesh in an HDF5 group but I haven't had the
time to figure out the HDF5 interface yet.

## Version
This is `v1.0.0`.

## Author
Sean McBane (<sean.mcbane@protonmail.com>)

## Copyright
Copyright 2020 The University of Texas at Austin.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal in the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.