
#include <random>
#include <iostream>
#include <vector>

#ifndef SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATIONS


#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATIONS
#ifndef SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATIONS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATIONS_DECLARE

#ifndef HEATEQUATION3D_H
#define HEATEQUATION3D_H

// Grid & Simulation Configuration
const int N = 80;               // Grid size (N x N)
const float L = 2.5f;           // Domain boundaries (-L to L)
const float D = 1.0f;           // Thermal diffusivity
const float dx = (2.0f * L) / (N - 1);
const float dt = 0.125f * (dx * dx) / D; // CFL explicit limit

enum Shape { CRESCENT, CYLINDER, BOX };

// Grid representation
typedef std::vector<std::vector<float>> Grid;

// Initialization
void resetGrid(Grid& U, Shape currentShape, float maxInitialTemp);

// Physics Solvers & Boundaries
void stepExplicit(Grid& U, Grid& nextU);
void stepImplicit(Grid& U, Grid& nextU);
void applyBoundaryConditions(Grid& U, bool isNeumann);

#endif // HEATEQUATION3D_H

#endif
#endif


#endif
