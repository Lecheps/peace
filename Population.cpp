#include "Population.h"

Population::Population() :
normal(RandGenerator::make_generator(DISTRIBUTION::NORMAL)),
uniform(RandGenerator::make_generator(DISTRIBUTION::UNIFORM))
{ 
        normParams = nullptr; 
        // normal->start_async_gen();
        // uniform->start_async_gen();
}

Population::~Population()
{
    delete normal;
    delete uniform;
};

Population::Population(const Population &p) : 
    normal(RandGenerator::make_generator(DISTRIBUTION::NORMAL)),
    uniform(RandGenerator::make_generator(DISTRIBUTION::UNIFORM))
{
    s = p.s;
    fg = p.fg;
    idx = p.idx;
    EC50s = p.EC50s;
    normParams = p.normParams;
    size = p.size;
    initPop = p.initPop;
    ts_length = p.ts_length;
    epsilon = p.epsilon;
    growth = p.growth;
}

 void Population::fillUntil(size_t pos, double& val)
{
    std::list<double> dummy(pos,val);
    c = std::move(dummy);
}


void Population::initRandPop()
{
    std::list<double> dummy(initPop,size);
    for (auto& i :  dummy) i = i * 2. / 3. * uniform->get() + i / 3.;
    c = std::move(dummy);    
}

double Population::getMass()
{
    double result = 0.;
    for (const auto& i : c) result +=i;
    return result;
}

size_t Population::getNumAlive()
{
    return c.size();
}

void Population::doTheEvolution()
{                
    for (auto it = c.begin(); it != c.end(); ++it)
    {
        *it += normParams->first / 2. * growth;
        if ( *it > (normParams->second  * normal->get() + normParams->first))
        {
            *it /= 2.;
            if (c.size() < potsize && uniform->get() >= probabilityOfDying)
            {
                c.insert(it,*it);
            } 
        }     
        
        if (uniform->get() < probabilityOfDying) 
        {
            it = c.erase(it);
            it--; 
        }
    }
}