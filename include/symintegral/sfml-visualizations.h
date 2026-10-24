
#include <random>
#include <iostream>
#include <vector>
#include <string>

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

#ifndef GENETICALGORITHMS2VARIABLES_H
#define GENETICALGORITHMS2VARIABLES_H

const float PI_VAL = 3.14159265f;
const int POP_SIZE = 60;

enum FunctionType { ROSENBROCK, RASTRIGIN, BOOTH, ACKLEY };

struct Domain {
    float xMin, xMax, yMin, yMax;
};

struct Point3D {
    float x, y, z;
};

struct Individual {
    float x, y;
    float fitness;
    std::string origin;
    void evaluate(FunctionType type);
};

// Global Utility Functions
Domain getDomain(FunctionType type);
std::string getFunctionName(FunctionType type);
float evaluateObjective(float x, float y, FunctionType type);
float getZHeight(float x, float y, FunctionType type);
float randFloat(float min, float max);

class GeneticAlgorithm {
public:
    std::vector<Individual> population;
    std::vector<Individual> historicalMinima;
    int generation;
    int countElite, countCrossover, countMutation;

    GeneticAlgorithm();
    void reset(FunctionType type);
    void evolve(FunctionType type, float mutationRate, float crossAlpha);

private:
    Individual tournamentSelect();
    void updateHistoricalMinima(FunctionType type, const Individual& bestCandidate);
};

#endif // GENETICALGORITHMS2VARIABLES_H

#pragma warning(disable : 4996)
#ifndef PSO2VARIABLES_H
#define PSO2VARIABLES_H

#include <vector>

// Structure representing a single particle in the 2D workspace
struct Particle {
    float x, y;         // Current space coordinates
    float vx, vy;       // Current component velocities
    float bestX, bestY; // Personal historical best coordinates
    float bestFitness;  // Lowest fitness value reached by this particle
};

// Main controller state container for the swarm
struct SwarmState {
    std::vector<Particle> particles;
    float globalBestX;
    float globalBestY;
    float globalBestFitness;
};

// PSO Hyperparameters
struct PSOHyperparams {
    float w;  // Inertia weight coefficient
    float c1; // Cognitive/Personal weight coefficient
    float c2; // Social/Global weight coefficient
};

// Typedef for generic 2D objective functions mapping (x,y) -> evaluation score
typedef float (*ObjectiveFunction)(float, float);

// Core mathematical optimization methods
void initSwarm(SwarmState& swarm, int numParticles, 
               float xMin, float xMax, float yMin, float yMax, 
               ObjectiveFunction f);

void updateSwarmMath(SwarmState& swarm, const PSOHyperparams& params,
                     float xMin, float xMax, float yMin, float yMax, 
                     ObjectiveFunction f);

#endif // PSO2VARIABLES_H

/*
#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <SFML/System/Vector2.hpp>
#include <SFML/Graphics/Color.hpp>
#include <vector>
#include <string>

enum MoleculeType {
    CAFFEINE = 0,
    ALCOHOL, // Ethanol
    WATER,
    GRAPHENE,
    DIAMOND,
    METHANE,
    COUNT
};

struct ElementProperties {
    std::string fullName;
    int atomicNumber;
    float atomicMass;
    std::string periodicGroup;
};

struct CustomColor {
    unsigned char r, g, b, a;
};

struct Atom {
    sf::Vector2f basePosition; // Original local coordinates
    sf::Vector2f position;     // Current transformed world coordinates
    sf::Color color;
    float radius;
    std::string label;
    ElementProperties props;
};

struct Bond {
    size_t atomA;
    size_t atomB;
    bool isDouble;
};

struct MoleculeData {
    std::string name;
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
};

// Core Chemistry Engine API
std::vector<MoleculeData> initializeMolecules();
void updateAtomRotations(std::vector<Atom>& atoms, float angle, const sf::Vector2f& center);
int checkAtomHover(const std::vector<Atom>& atoms, const sf::Vector2f& mousePos);

#endif // CHEMISTRY_H
*/
#endif
#endif


#endif
