// Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

SetFactory("OpenCASCADE");

// outer box
Lx1 = 5;
Ly1 = 4;
Lz1 = 3;
Box(1) = {0, 0, 0, Lx1, Ly1, Lz1};

// inner box 2
Lx2 = 1;
Ly2 = 1;
Lz2 = 1;
ox2 = 1; oy2 = 1; oz2 = 1; // origin
Box(2) = {ox2, oy2, oz2, Lx2, Ly2, Lz2};

// inner box 3
Lx3 = 1;
Ly3 = 1;
Lz3 = 1;
ox3 = 3; oy3 = 1; oz3 = 1; // origin
Box(3) = {ox3, oy3, oz3, Lx3, Ly3, Lz3};

// Subtract inner boxes from outer box
BooleanDifference(4) = { Volume{1}; Delete; }{ Volume{2, 3}; };

// physical volumes
Physical Volume("volume1", 1) = {4};
Physical Volume("volume2", 2) = {2};
Physical Volume("volume3", 3) = {3};

// tolerance for surface selection
eps = 0.01;

// outer box physical surfaces 
Physical Surface("face1", 1) = { Surface In BoundingBox{
    0-eps, 0-eps, 0-eps, Lx1+eps, Ly1+eps, 0+eps
}};
Physical Surface("face2", 2) = { Surface In BoundingBox{
    0-eps, 0-eps, 0-eps, Lx1+eps, 0+eps, Lz1+eps
}};
Physical Surface("face3", 3) = { Surface In BoundingBox{
    Lx1-eps, 0-eps, 0-eps, Lx1+eps, Ly1+eps, Lz1+eps
}};
Physical Surface("face4", 4) = { Surface In BoundingBox{
    0-eps, Ly1-eps, 0-eps, Lx1+eps, Ly1+eps, Lz1+eps
}};
Physical Surface("face5", 5) = { Surface In BoundingBox{
    0-eps, 0-eps, 0-eps, 0+eps, Ly1+eps, Lz1+eps
}};
Physical Surface("face6", 6) = { Surface In BoundingBox{
    0-eps, 0-eps, Lz1-eps, Lx1+eps, Ly1+eps, Lz1+eps
}};

// inner box 2 physical surfaces
Physical Surface("box2_xmin", 7) = { Surface In BoundingBox{
    ox2-eps, oy2-eps, oz2-eps, ox2+eps, oy2+Ly2+eps, oz2+Lz2+eps
}};
Physical Surface("box2_xmax", 8) = { Surface In BoundingBox{
    ox2+Lx2-eps, oy2-eps, oz2-eps, ox2+Lx2+eps, oy2+Ly2+eps, oz2+Lz2+eps
}};
Physical Surface("box2_ymin", 9) = { Surface In BoundingBox{
    ox2-eps, oy2-eps, oz2-eps, ox2+Lx2+eps, oy2+eps, oz2+Lz2+eps
}};
Physical Surface("box2_ymax", 10) = { Surface In BoundingBox{
    ox2-eps, oy2+Ly2-eps, oz2-eps, ox2+Lx2+eps, oy2+Ly2+eps, oz2+Lz2+eps
}};
Physical Surface("box2_zmin", 11) = { Surface In BoundingBox{
    ox2-eps, oy2-eps, oz2-eps, ox2+Lx2+eps, oy2+Ly2+eps, oz2+eps
}};
Physical Surface("box2_zmax", 12) = { Surface In BoundingBox{
    ox2-eps, oy2-eps, oz2+Lz2-eps, ox2+Lx2+eps, oy2+Ly2+eps, oz2+Lz2+eps
}};

// inner box 3 physical surfaces
Physical Surface("box3_xmin", 13) = { Surface In BoundingBox{
    ox3-eps, oy3-eps, oz3-eps, ox3+eps, oy3+Ly3+eps, oz3+Lz3+eps
}};
Physical Surface("box3_xmax", 14) = { Surface In BoundingBox{
    ox3+Lx3-eps, oy3-eps, oz3-eps, ox3+Lx3+eps, oy3+Ly3+eps, oz3+Lz3+eps
}};
Physical Surface("box3_ymin", 15) = { Surface In BoundingBox{
    ox3-eps, oy3-eps, oz3-eps, ox3+Lx3+eps, oy3+eps, oz3+Lz3+eps
}};
Physical Surface("box3_ymax", 16) = { Surface In BoundingBox{
    ox3-eps, oy3+Ly3-eps, oz3-eps, ox3+Lx3+eps, oy3+Ly3+eps, oz3+Lz3+eps
}};
Physical Surface("box3_zmin", 17) = { Surface In BoundingBox{
    ox3-eps, oy3-eps, oz3-eps, ox3+Lx3+eps, oy3+Ly3+eps, oz3+eps
}};
Physical Surface("box3_zmax", 18) = { Surface In BoundingBox{
    ox3-eps, oy3-eps, oz3+Lz3-eps, ox3+Lx3+eps, oy3+Ly3+eps, oz3+Lz3+eps
}};

// meshing config
Mesh.CharacteristicLengthMin = 100;
Mesh.CharacteristicLengthMax = 100;
Mesh.ElementOrder = 1;
Mesh.Algorithm3D = 4; // 4 = Netgen algorithm

// Generate 3D mesh
Mesh 3;

// BDF export — small field format (0 = free, 1 = small, 2 = large)
Mesh.BdfFieldFormat = 1;
// Mesh.SaveAll = 1; // include all elements, not only those in Physical groups
Save "test1.bdf";
