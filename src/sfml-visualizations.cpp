/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATION_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_SFMLVISUALIZATION_DEFINE
/*

	Genetic Algorithms

*/


void Individual::evaluate(FunctionType type) {
    fitness = 1.0f / (1.0f + evaluateObjective(x, y, type));
}

Domain getDomain(FunctionType type) {
    switch (type) {
        case RASTRIGIN: return { -5.12f, 5.12f, -5.12f, 5.12f };
        case BOOTH:     return { -10.0f, 10.0f, -10.0f, 10.0f };
        case ACKLEY:    return { -5.0f,  5.0f,  -5.0f,  5.0f  };
        case ROSENBROCK:
        default:        return { -2.0f,  2.0f,  -1.0f,  3.0f  };
    }
}

std::string getFunctionName(FunctionType type) {
    switch (type) {
        case ROSENBROCK: return "Rosenbrock (Valley)";
        case RASTRIGIN:  return "Rastrigin (Highly Multimodal)";
        case BOOTH:      return "Booth (Plateau)";
        case ACKLEY:     return "Ackley (Deep Central Hole)";
        default:         return "Unknown";
    }
}

float evaluateObjective(float x, float y, FunctionType type) {
    switch (type) {
        case ROSENBROCK: 
            return std::pow(1.0f - x, 2) + 100.0f * std::pow(y - std::pow(x, 2), 2);
            
        case RASTRIGIN: 
            return 20.0f + (std::pow(x, 2) - 10.0f * std::cos(2.0f * PI_VAL * x)) 
                         + (std::pow(y, 2) - 10.0f * std::cos(2.0f * PI_VAL * y));
                         
        case BOOTH: 
            return std::pow(x + 2.0f * y - 7.0f, 2) + std::pow(2.0f * x + y - 5.0f, 2);
            
        case ACKLEY: {
            float term1 = -20.0f * std::exp(-0.2f * std::sqrt(0.5f * (std::pow(x, 2) + std::pow(y, 2))));
            float term2 = -std::exp(0.5f * (std::cos(2.0f * PI_VAL * x) + std::cos(2.0f * PI_VAL * y)));
            return term1 + term2 + 20.0f + std::exp(1.0f);
        }
        
        default: 
            return 0.0f;
    }
}

float getZHeight(float x, float y, FunctionType type) {
    float val = evaluateObjective(x, y, type);
    if (type == ROSENBROCK) return std::log1p(val) * 0.15f;
    if (type == RASTRIGIN)  return val * 0.025f;
    if (type == BOOTH)      return std::log1p(val) * 0.15f;
    if (type == ACKLEY)     return val * 0.08f;
    return val;
}

float randFloat(float min, float max) {
    return min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min)));
}

GeneticAlgorithm::GeneticAlgorithm() : generation(0), countElite(0), countCrossover(0), countMutation(0) {}

void GeneticAlgorithm::reset(FunctionType type) {
    generation = 0;
    population.clear();
    historicalMinima.clear();
    
    Domain d = getDomain(type);
    for (int i = 0; i < POP_SIZE; ++i) {
        Individual ind;
        ind.x = randFloat(d.xMin, d.xMax);
        ind.y = randFloat(d.yMin, d.yMax);
        ind.origin = "Initial";
        ind.evaluate(type);
        population.push_back(ind);
    }
    updateHistoricalMinima(type, population[0]);
}

void GeneticAlgorithm::evolve(FunctionType type, float mutationRate, float crossAlpha) {
    std::vector<Individual> newPopulation;
    countElite = 0;
    countCrossover = 0;
    countMutation = 0;
    Domain d = getDomain(type);
    
    std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness > b.fitness;
    });
    
    // Elitism
    for (int i = 0; i < 2; ++i) {
        Individual elite = population[i];
        elite.origin = "Elite";
        newPopulation.push_back(elite);
        countElite++;
    }
    
    // Crossover and Mutation
    while (newPopulation.size() < POP_SIZE) {
        Individual parent1 = tournamentSelect();
        Individual parent2 = tournamentSelect();
        
        Individual child;
        child.x = crossAlpha * parent1.x + (1.0f - crossAlpha) * parent2.x;
        child.y = crossAlpha * parent1.y + (1.0f - crossAlpha) * parent2.y;
        child.origin = "Crossover";
        countCrossover++;
        
        if (randFloat(0.0f, 1.0f) < mutationRate) {
            float rangeX = (d.xMax - d.xMin) * 0.05f;
            float rangeY = (d.yMax - d.yMin) * 0.05f;
            child.x += randFloat(-rangeX, rangeX);
            child.y += randFloat(-rangeY, rangeY);
            
            // Clamp within boundary limits
            child.x = std::max(d.xMin, std::min(d.xMax, child.x));
            child.y = std::max(d.yMin, std::min(d.yMax, child.y));
            child.origin = "Mutated";
            countMutation++;
        }
        
        child.evaluate(type);
        newPopulation.push_back(child);
    }
    
    population = newPopulation;
    generation++;
    
    std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness > b.fitness;
    });
    
    updateHistoricalMinima(type, population[0]);
}

Individual GeneticAlgorithm::tournamentSelect() {
    int i1 = rand() % POP_SIZE;
    int i2 = rand() % POP_SIZE;
    return (population[i1].fitness > population[i2].fitness) ? population[i1] : population[i2];
}

void GeneticAlgorithm::updateHistoricalMinima(FunctionType type, const Individual& bestCandidate) {
    if (historicalMinima.empty()) {
        historicalMinima.push_back(bestCandidate);
        return;
    }
    
    bool isUnique = true;
    for(const auto& entry : historicalMinima) {
        float dist = std::sqrt(std::pow(entry.x - bestCandidate.x, 2) + std::pow(entry.y - bestCandidate.y, 2));
        if (dist < 0.05f) {
            isUnique = false;
            break;
        }
    }
    
    if (isUnique || evaluateObjective(bestCandidate.x, bestCandidate.y, type) < evaluateObjective(historicalMinima.back().x, historicalMinima.back().y, type)) {
        historicalMinima.push_back(bestCandidate);
        
        std::sort(historicalMinima.begin(), historicalMinima.end(), [type](const Individual& a, const Individual& b) {
            return evaluateObjective(a.x, a.y, type) < evaluateObjective(b.x, b.y, type);
        });
        
        if (historicalMinima.size() > 3) {
            historicalMinima.resize(3);
        }
    }
}


/*

	PSO

*/

void initSwarm(SwarmState& swarm, int numParticles, 
               float xMin, float xMax, float yMin, float yMax, 
               ObjectiveFunction f) {
    
    swarm.particles.resize(numParticles);
    swarm.globalBestFitness = 1e9f; // Seeking global minimums
    swarm.globalBestX = 0.0f;
    swarm.globalBestY = 0.0f;

    for (int i = 0; i < numParticles; ++i) {
        // Uniform random positioning across the mathematical workspace bounds
        swarm.particles[i].x = xMin + (static_cast<float>(std::rand()) / RAND_MAX) * (xMax - xMin);
        swarm.particles[i].y = yMin + (static_cast<float>(std::rand()) / RAND_MAX) * (yMax - yMin);
        
        // Minor randomized initial trajectories
        swarm.particles[i].vx = ((static_cast<float>(std::rand()) / RAND_MAX) - 0.5f) * 0.5f;
        swarm.particles[i].vy = ((static_cast<float>(std::rand()) / RAND_MAX) - 0.5f) * 0.5f;
        
        swarm.particles[i].bestX = swarm.particles[i].x;
        swarm.particles[i].bestY = swarm.particles[i].y;
        swarm.particles[i].bestFitness = f(swarm.particles[i].x, swarm.particles[i].y);

        if (swarm.particles[i].bestFitness < swarm.globalBestFitness) {
            swarm.globalBestFitness = swarm.particles[i].bestFitness;
            swarm.globalBestX = swarm.particles[i].bestX;
            swarm.globalBestY = swarm.particles[i].bestY;
        }
    }
}

void updateSwarmMath(SwarmState& swarm, const PSOHyperparams& params,
                     float xMin, float xMax, float yMin, float yMax, 
                     ObjectiveFunction f) {
    
    for (size_t i = 0; i < swarm.particles.size(); ++i) {
        float r1 = static_cast<float>(std::rand()) / RAND_MAX;
        float r2 = static_cast<float>(std::rand()) / RAND_MAX;

        // Velocity kinematic adjustments incorporating personal and swarm best metrics
        swarm.particles[i].vx = params.w * swarm.particles[i].vx 
                                + params.c1 * r1 * (swarm.particles[i].bestX - swarm.particles[i].x) 
                                + params.c2 * r2 * (swarm.globalBestX - swarm.particles[i].x);
        
        swarm.particles[i].vy = params.w * swarm.particles[i].vy 
                                + params.c1 * r1 * (swarm.particles[i].bestY - swarm.particles[i].y) 
                                + params.c2 * r2 * (swarm.globalBestY - swarm.particles[i].y);

        // Position translation step
        swarm.particles[i].x += swarm.particles[i].vx;
        swarm.particles[i].y += swarm.particles[i].vy;

        // Strictly enforce physical system boundary constraints
        if (swarm.particles[i].x < xMin) { swarm.particles[i].x = xMin; swarm.particles[i].vx = 0; }
        if (swarm.particles[i].x > xMax) { swarm.particles[i].x = xMax; swarm.particles[i].vx = 0; }
        if (swarm.particles[i].y < yMin) { swarm.particles[i].y = yMin; swarm.particles[i].vy = 0; }
        if (swarm.particles[i].y > yMax) { swarm.particles[i].y = yMax; swarm.particles[i].vy = 0; }

        float currentFitness = f(swarm.particles[i].x, swarm.particles[i].y);

        // Check and evaluate for personal baseline updates
        if (currentFitness < swarm.particles[i].bestFitness) {
            swarm.particles[i].bestFitness = currentFitness;
            swarm.particles[i].bestX = swarm.particles[i].x;
            swarm.particles[i].bestY = swarm.particles[i].y;
        }

        // Check and evaluate for collective global baseline updates
        if (currentFitness < swarm.globalBestFitness) {
            swarm.globalBestFitness = currentFitness;
            swarm.globalBestX = swarm.particles[i].x;
            swarm.globalBestY = swarm.particles[i].y;
        }
    }
}


/*

	PDE

*/
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

/*

	CHEMISTRY



const sf::Color COLOR_CARBON(50, 50, 50);
const sf::Color COLOR_NITROGEN(30, 144, 255);
const sf::Color COLOR_OXYGEN(255, 50, 50);
const sf::Color COLOR_HYDROGEN(235, 235, 235);

const ElementProperties PROP_HYDROGEN = {"Hydrogen", 1, 1.008f, "Group 1 (Reactive Non-metal)"};
const ElementProperties PROP_CARBON   = {"Carbon", 6, 12.011f, "Group 14 (Carbon Group)"};
const ElementProperties PROP_NITROGEN = {"Nitrogen", 7, 14.007f, "Group 15 (Pnictogens)"};
const ElementProperties PROP_OXYGEN   = {"Oxygen", 8, 15.999f, "Group 16 (Chalcogens)"};

std::vector<MoleculeData> initializeMolecules() {
    std::vector<MoleculeData> dataList(MoleculeType::COUNT);

    // 0. CAFFEINE
    dataList[CAFFEINE].name = "Caffeine (C8H10N4O2)";
    auto& cafAtoms = dataList[CAFFEINE].atoms;
    auto& cafBonds = dataList[CAFFEINE].bonds;
    cafAtoms.push_back({sf::Vector2f(-50.f, -86.6f),  {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(50.f, -86.6f),   {0,0}, COLOR_NITROGEN, 18.f, "N", PROP_NITROGEN});
    cafAtoms.push_back({sf::Vector2f(100.f, 0.f),     {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(50.f, 86.6f),    {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(-50.f, 86.6f),   {0,0}, COLOR_NITROGEN, 18.f, "N", PROP_NITROGEN});
    cafAtoms.push_back({sf::Vector2f(-100.f, 0.f),    {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(180.f, -40.f),   {0,0}, COLOR_NITROGEN, 18.f, "N", PROP_NITROGEN});
    cafAtoms.push_back({sf::Vector2f(230.f, 43.3f),   {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(150.f, 110.f),   {0,0}, COLOR_NITROGEN, 18.f, "N", PROP_NITROGEN});
    cafAtoms.push_back({sf::Vector2f(-200.f, 0.f),    {0,0}, COLOR_OXYGEN,   22.f, "O", PROP_OXYGEN});
    cafAtoms.push_back({sf::Vector2f(-90.f, -160.f),  {0,0}, COLOR_OXYGEN,   22.f, "O", PROP_OXYGEN});
    cafAtoms.push_back({sf::Vector2f(70.f, -200.f),   {0,0}, COLOR_CARBON,   20.f, "CH3", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(-90.f, 180.f),   {0,0}, COLOR_CARBON,   20.f, "CH3", PROP_CARBON});
    cafAtoms.push_back({sf::Vector2f(230.f, -140.f),  {0,0}, COLOR_CARBON,   20.f, "CH3", PROP_CARBON});
    cafBonds.push_back({0, 1, false}); cafBonds.push_back({1, 2, false}); cafBonds.push_back({2, 3, true});
    cafBonds.push_back({3, 4, false}); cafBonds.push_back({4, 5, false}); cafBonds.push_back({5, 0, false});
    cafBonds.push_back({2, 6, false}); cafBonds.push_back({6, 7, false}); cafBonds.push_back({7, 8, true});
    cafBonds.push_back({8, 3, false}); cafBonds.push_back({5, 9, true});  cafBonds.push_back({0, 10, true});
    cafBonds.push_back({1, 11, false});cafBonds.push_back({4, 12, false});cafBonds.push_back({6, 13, false});

    // 1. ALCOHOL (Ethanol - C2H5OH)
    dataList[ALCOHOL].name = "Alcohol / Ethanol (C2H5OH)";
    auto& alcAtoms = dataList[ALCOHOL].atoms;
    auto& alcBonds = dataList[ALCOHOL].bonds;
    alcAtoms.push_back({sf::Vector2f(-80.f, 0.f),     {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    alcAtoms.push_back({sf::Vector2f(40.f, 0.f),      {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    alcAtoms.push_back({sf::Vector2f(140.f, 70.f),    {0,0}, COLOR_OXYGEN,   22.f, "O", PROP_OXYGEN});
    alcAtoms.push_back({sf::Vector2f(230.f, 40.f),    {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcAtoms.push_back({sf::Vector2f(-160.f, 0.f),    {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcAtoms.push_back({sf::Vector2f(-80.f, -80.f),   {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcAtoms.push_back({sf::Vector2f(-80.f, 80.f),    {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcAtoms.push_back({sf::Vector2f(40.f, -80.f),    {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcAtoms.push_back({sf::Vector2f(40.f, 80.f),     {0,0}, COLOR_HYDROGEN, 14.f, "H", PROP_HYDROGEN});
    alcBonds.push_back({0, 1, false}); alcBonds.push_back({1, 2, false}); alcBonds.push_back({2, 3, false});
    alcBonds.push_back({0, 4, false}); alcBonds.push_back({0, 5, false}); alcBonds.push_back({0, 6, false});
    alcBonds.push_back({1, 7, false}); alcBonds.push_back({1, 8, false});

    // 2. WATER (H2O)
    dataList[WATER].name = "Water (H2O)";
    auto& watAtoms = dataList[WATER].atoms;
    auto& watBonds = dataList[WATER].bonds;
    watAtoms.push_back({sf::Vector2f(0.f, -30.f),     {0,0}, COLOR_OXYGEN,   24.f, "O", PROP_OXYGEN});
    watAtoms.push_back({sf::Vector2f(-90.f, 50.f),    {0,0}, COLOR_HYDROGEN, 16.f, "H", PROP_HYDROGEN});
    watAtoms.push_back({sf::Vector2f(90.f, 50.f),     {0,0}, COLOR_HYDROGEN, 16.f, "H", PROP_HYDROGEN});
    watBonds.push_back({0, 1, false}); watBonds.push_back({0, 2, false});

    // 3. GRAPHENE
    dataList[GRAPHENE].name = "Graphene Sheet Segment (Carbon Lattice)";
    auto& graAtoms = dataList[GRAPHENE].atoms;
    auto& graBonds = dataList[GRAPHENE].bonds;
    float r = 50.f; 
    for (int row = -1; row <= 1; ++row) {
        for (int col = -2; col <= 2; ++col) {
            float cx = col * 1.5f * r;
            float cy = row * std::sqrt(3.f) * r + (std::abs(col) % 2 == 1 ? std::sqrt(3.f)*0.5f*r : 0.f);
            if (cx*cx + cy*cy < 250.f*250.f) { 
                graAtoms.push_back({sf::Vector2f(cx, cy), {0,0}, COLOR_CARBON, 16.f, "C", PROP_CARBON});
            }
        }
    }
    for (size_t i = 0; i < graAtoms.size(); ++i) {
        for (size_t j = i + 1; j < graAtoms.size(); ++j) {
            sf::Vector2f d = graAtoms[i].basePosition - graAtoms[j].basePosition;
            float dist = std::sqrt(d.x*d.x + d.y*d.y);
            if (dist > r * 0.95f && dist < r * 1.05f) {
                graBonds.push_back({i, j, (i % 3 == 0 && j % 2 == 0)}); 
            }
        }
    }

    // 4. DIAMOND
    dataList[DIAMOND].name = "Diamond Crystal Lattice Segment";
    auto& diaAtoms = dataList[DIAMOND].atoms;
    auto& diaBonds = dataList[DIAMOND].bonds;
    diaAtoms.push_back({sf::Vector2f(0.f, 0.f),       {0,0}, COLOR_CARBON,   20.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(0.f, -90.f),     {0,0}, COLOR_CARBON,   18.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(-80.f, 40.f),    {0,0}, COLOR_CARBON,   18.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(80.f, 40.f),     {0,0}, COLOR_CARBON,   18.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(20.f, 60.f),     {0,0}, COLOR_CARBON,   16.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(-160.f, 80.f),   {0,0}, COLOR_CARBON,   16.f, "C", PROP_CARBON});
    diaAtoms.push_back({sf::Vector2f(160.f, 80.f),    {0,0}, COLOR_CARBON,   16.f, "C", PROP_CARBON});
    diaBonds.push_back({0, 1, false}); diaBonds.push_back({0, 2, false}); diaBonds.push_back({0, 3, false}); 
    diaBonds.push_back({0, 4, false}); diaBonds.push_back({2, 5, false}); diaBonds.push_back({3, 6, false});

    // 5. METHANE
    dataList[METHANE].name = "Methane (CH4)";
    auto& metAtoms = dataList[METHANE].atoms;
    auto& metBonds = dataList[METHANE].bonds;
    metAtoms.push_back({sf::Vector2f(0.f, 0.f),       {0,0}, COLOR_CARBON,   22.f, "C", PROP_CARBON});
    metAtoms.push_back({sf::Vector2f(0.f, -100.f),    {0,0}, COLOR_HYDROGEN, 15.f, "H", PROP_HYDROGEN});
    metAtoms.push_back({sf::Vector2f(-95.f, 45.f),    {0,0}, COLOR_HYDROGEN, 15.f, "H", PROP_HYDROGEN});
    metAtoms.push_back({sf::Vector2f(95.f, 45.f),     {0,0}, COLOR_HYDROGEN, 15.f, "H", PROP_HYDROGEN});
    metAtoms.push_back({sf::Vector2f(30.f, 35.f),     {0,0}, COLOR_HYDROGEN, 13.f, "H", PROP_HYDROGEN});
    metBonds.push_back({0, 1, false}); metBonds.push_back({0, 2, false}); 
    metBonds.push_back({0, 3, false}); metBonds.push_back({0, 4, false});

    return dataList;
}

// Math/Geometry Transformation Logic
void updateAtomRotations(std::vector<Atom>& atoms, float angle, const sf::Vector2f& center) {
    float cosA = std::cos(angle);
    float sinA = std::sin(angle);
    for (auto& atom : atoms) {
        float rx = atom.basePosition.x * cosA - atom.basePosition.y * sinA;
        float ry = atom.basePosition.x * sinA + atom.basePosition.y * cosA;
        atom.position = center + sf::Vector2f(rx, ry);
    }
}

// Distance checking lookup validation
int checkAtomHover(const std::vector<Atom>& atoms, const sf::Vector2f& mousePos) {
    for (size_t i = 0; i < atoms.size(); ++i) {
        sf::Vector2f diff = mousePos - atoms[i].position;
        if (std::sqrt(diff.x * diff.x + diff.y * diff.y) <= atoms[i].radius) {
            return static_cast<int>(i);
        }
    }
    return -1;
}
*/

#endif
#endif