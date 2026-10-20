/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATION_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATION_DEFINE

void resetGrid(Grid& U, Shape currentShape, float maxInitialTemp) 
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            float x = -L + i * dx;
            float y = -L + j * dx;
            if (currentShape == CRESCENT) {
                float c = 2.0f, d = -1.0f, e = 1.0f, f = 0.5f, k = 1.2f, shift = 10.0f;
                float pattern = std::pow(c*c - std::pow(x/e - d, 2) - std::pow(y/f, 2), 2) + k * std::pow(c + d - x/e, 3) - shift;
                U[i][j] = (pattern < 0.0f) ? maxInitialTemp : 0.0f;
            } else if (currentShape == CYLINDER) {
                U[i][j] = (x*x + y*y < 1.0f) ? maxInitialTemp : 0.0f;
            } else if (currentShape == BOX) {
                U[i][j] = (std::abs(x) < 0.8f && std::abs(y) < 0.8f) ? maxInitialTemp : 0.0f;
            }
        }
    }
}

void stepExplicit(Grid& U, Grid& nextU) 
{
    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            float laplacian = (U[i+1][j] + U[i-1][j] + U[i][j+1] + U[i][j-1] - 4.0f * U[i][j]) / (dx * dx);
            nextU[i][j] = U[i][j] + dt * D * laplacian;
        }
    }
    U = nextU; // Swapped outward safely after all cell states calculate
}

void stepImplicit(Grid& U, Grid& nextU) 
{
    float current_dt = dt * 3.0f;
    float r = (D * current_dt) / (2.0f * dx * dx);
    Grid B(N, std::vector<float>(N, 0.0f));

    for (int i = 1; i < N - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            B[i][j] = U[i][j] * (1.0f - 4.0f * r) + r * (U[i+1][j] + U[i-1][j] + U[i][j+1] + U[i][j-1]);
        }
    }
    
    nextU = U;
    for (int it = 0; it < 8; ++it) {
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                nextU[i][j] = (B[i][j] + r * (nextU[i+1][j] + nextU[i-1][j] + nextU[i][j+1] + nextU[i][j-1])) / (1.0f + 4.0f * r);
            }
        }
    }
    U = nextU;
}

void applyBoundaryConditions(Grid& U, bool isNeumann) 
{
    for (int k = 0; k < N; ++k) {
        if (isNeumann) {
            U[0][k]   = U[1][k];     // Left Boundary 
            U[N-1][k] = U[N-2][k];   // Right Boundary
            U[k][0]   = U[k][1];     // Bottom Boundary
            U[k][N-1] = U[k][N-2];   // Top Boundary
        } else {
            U[0][k]   = 0.0f;        // Left Boundary
            U[N-1][k] = 0.0f;        // Right Boundary
            U[k][0]   = 0.0f;        // Bottom Boundary
            U[k][N-1] = 0.0f;        // Top Boundary
        }
    }
}

#endif
#endif