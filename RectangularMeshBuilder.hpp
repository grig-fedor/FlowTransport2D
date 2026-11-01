#pragma once
#include "inmost.h"

class RectangularMeshBuilder {
public:
	static INMOST::Mesh* Generate(double Lx, double Ly, int Nx, int Ny, double x0, double y0);
};
