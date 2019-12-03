#pragma once

#include "RandGenerator.h"
#include<list>

class Population
{
public:

    Population();

    ~Population();
   

    Population(const Population &p); 

    unsigned short int s = std::numeric_limits<unsigned short int>::max();
    unsigned short int fg = std::numeric_limits<unsigned short int>::max();
    size_t idx = std::numeric_limits<size_t>::max();

    double EC50s = 0.;
    double Growth_adj_tolerance = 0.;
    std::pair<double, double>* normParams;
    double size = 0.;
    size_t initPop = 0;
    size_t ts_length = 0;
    double epsilon = 0.;
    double growth = 0.;
    RandGenerator* normal;
    RandGenerator* uniform;
    double probabilityOfDying;
    
    static const size_t potsize = 10000;
    
    // std::array<double,potsize> c;
    std::list<double> c;

    void fillUntil(size_t pos, double& val);

    void initRandPop();
    
    double getMass();    

    size_t getNumAlive();
    
    void doTheEvolution();
};