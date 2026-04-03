// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

// Define box dimensions
L = 5; // Length (x-direction)
W = 4; // Width (y-direction)
H = 3; // Height (z-direction)

// Define corner points of the box
Point(1) = {0, 0, 0}; // Origin
Point(2) = {L, 0, 0};
Point(3) = {L, W, 0};
Point(4) = {0, W, 0};
Point(5) = {0, 0, H};
Point(6) = {L, 0, H};
Point(7) = {L, W, H};
Point(8) = {0, W, H};

// Connect points to create edges
Line(1) = {1, 2}; // Bottom face
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {1, 5}; // Vertical edges
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6}; // Top face
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 5};

// Create surfaces
Line Loop(1) = {1, 2, 3, 4}; // Bottom face
Plane Surface(1) = {1};

Line Loop(2) = {5, 9, -6, -1}; // Front face
Plane Surface(2) = {2};

Line Loop(3) = {6, 10, -7, -2}; // Right face
Plane Surface(3) = {3};

Line Loop(4) = {7, 11, -8, -3}; // Back face
Plane Surface(4) = {4};

Line Loop(5) = {8, 12, -5, -4}; // Left face
Plane Surface(5) = {5};

Line Loop(6) = {9, 10, 11, 12}; // Top face
Plane Surface(6) = {6};

// Define volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Set meshing parameters
Mesh.CharacteristicLengthMax = 0.4;
Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 4; // 4 = Netgen algorithm (good for tets)

// Generate 3D tetrahedral mesh
Physical Volume("BoxVolume", 1) = {1};
