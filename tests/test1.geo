// Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

SetFactory("OpenCASCADE");

// outer box
Lx0 = 5; // x
Ly0 = 4; // y
Lz0 = 3; // z
Box(1) = {0, 0, 0, Lx0, Ly0, Lz0};

// inner box
Lx1 = 2; // x
Ly1 = 2; // y
Lz1 = 1; // z
Box(2) = {1, 1, 1, Lx1, Ly1, Lz1};

// Subtract inner box from outer box
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; };

// physical volumes
Physical Volume("volume1", 1) = {3};
Physical Volume("volume2", 2) = {2};

// tolerance for surface selection
eps = 0.01;

// physical surfaces (using bounding boxes)
Physical Surface("face1", 1) = { Surface In 
    BoundingBox{0-eps, 0-eps, 0-eps, Lx0+eps, Ly0+eps, 0+eps}
};
Physical Surface("face2", 2) = { Surface In 
    BoundingBox{0-eps, 0-eps, 0-eps, Lx0+eps, 0+eps, Lz0+eps}
};
Physical Surface("face3", 3) = { Surface In 
    BoundingBox{Lx0-eps, 0-eps, 0-eps, Lx0+eps, Ly0+eps, Lz0+eps}
};
Physical Surface("face4", 4) = { Surface In 
    BoundingBox{0-eps, Ly0-eps, 0-eps, Lx0+eps, Ly0+eps, Lz0+eps}
};
Physical Surface("face5", 5) = { Surface In 
    BoundingBox{0-eps, 0-eps, 0-eps, 0+eps, Ly0+eps, Lz0+eps}
};
Physical Surface("face6", 6) = { Surface In 
    BoundingBox{0-eps, 0-eps, Lz0-eps, Lx0+eps, Ly0+eps, Lz0+eps}
};

// meshing config
Mesh.CharacteristicLengthMax = 0.5;
Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 4; // 4 = Netgen algorithm
