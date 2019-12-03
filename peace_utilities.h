#pragma once

#ifdef DEBUG_ENABLE
#define D(x) x
#else 
#define D(x)
#endif

#include "RandGenerator.h"
#include "Parameters.h"
#include "Population.h"
#include <boost/range/combine.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/lexical_cast.hpp>

#include <math.h>
#include <array>
#include <chrono>
#include <limits>
#include <thread>
#include <vector>
#include <list>
#include <iomanip> 
#include <omp.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>


using Container = boost::multi_index_container<
    Population*, // the data type stored
    boost::multi_index::indexed_by
    <
        boost::multi_index::hashed_non_unique
        <  
            boost::multi_index::tag<struct s>, 
            boost::multi_index::member<Population, unsigned short int, & Population::s> 
        >,
        boost::multi_index::hashed_non_unique
        <  
            boost::multi_index::tag<struct fg>, 
            boost::multi_index::member<Population, unsigned short int, & Population::fg >
        >,
        boost::multi_index::hashed_unique
        <
            boost::multi_index::tag<struct idx>,
            boost::multi_index::member<Population, size_t, & Population::idx >
        >,
        boost::multi_index::ordered_unique
        <
            boost::multi_index::tag<struct s_and_fg>,
            boost::multi_index::composite_key
            <
                Population*,
                boost::multi_index::member<Population,unsigned short int,&Population::s>,
                boost::multi_index::member<Population,unsigned short int,&Population::fg>
            >,
            boost::multi_index::composite_key_compare
            <
                std::less<unsigned short int>,     
                std::less<unsigned short int> 
            >
        >
    >
>;


struct ParB
{
    unsigned short int s;
    unsigned short int fg;
    size_t idx;
    size_t t;
    size_t reit;
    double mass;
    size_t num_alive;

    ParB(const ParB& p)
    {
        s = p.s;
        fg = p.fg;
        idx = p.idx;
        t = p.t;
        reit = p.reit;
        mass = p.mass;
        num_alive = p.num_alive;
    }

    ParB(unsigned short int s_init,
         unsigned short int fg_init,
         size_t t_init,
         size_t reit_init) :
         s(s_init), fg(fg_init), t(t_init), reit(reit_init)  
    {
        mass = 0.;
        num_alive = 0.;
    }
};



using ParBContainer = boost::multi_index_container
<
    ParB*, // the data type stored
    boost::multi_index::indexed_by
        <
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct s_par>,
                boost::multi_index::member<ParB, unsigned short int, &ParB::s>
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct fg_par>,
                boost::multi_index::member<ParB, unsigned short int, &ParB::fg >
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct idx_par>,
                boost::multi_index::member<ParB, size_t, &ParB::idx >
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct t_par>,
                boost::multi_index::member<ParB, size_t, &ParB::t >
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct t_and_reit>,
                boost::multi_index::composite_key
                <
                    ParB,
                    boost::multi_index::member<ParB,size_t,&ParB::t>,
                    boost::multi_index::member<ParB,size_t,&ParB::reit>
                    
                >//,
                // boost::multi_index::composite_key_compare
                // <
                //     std::less<size_t>,
                //     std::less<size_t>
                // >
            >,
            boost::multi_index::hashed_unique
            <
                boost::multi_index::tag<struct all>,
                boost::multi_index::composite_key
                <
                    ParB,
                    boost::multi_index::member<ParB,unsigned short int,&ParB::s>,
                    boost::multi_index::member<ParB,unsigned short int,&ParB::fg>,
                    boost::multi_index::member<ParB,size_t,&ParB::t>,
                    boost::multi_index::member<ParB,size_t,&ParB::reit>
                    
                >//,
                // boost::multi_index::composite_key_compare
                // <
                //     std::less<size_t>,
                //     std::less<size_t>
                // >
            >
        >
>;


class ScopedGILRelease
{
public:
    inline ScopedGILRelease()
    {
        m_thread_state = PyEval_SaveThread();
    }

    inline ~ScopedGILRelease()
    {
        PyEval_RestoreThread(m_thread_state);
        m_thread_state = NULL;
    }

private:
    PyThreadState * m_thread_state;
};



std::vector<ParB> run(Parameters par,bool fromPython=true)
{
    if (fromPython) ScopedGILRelease scoped; 
    auto start = std::chrono::high_resolution_clock::now();
    double L1 = 1.0;
    par.setDeathProbabilities();
    std::vector<double> Unbiased_GR = par.sizeClasses;
    for (auto& i : Unbiased_GR ) i = par.Fi_1 * pow(i + par.b, par.g);
    par.setDeathProbabilities();

    D(for (auto& i : Unbiased_GR) std::cout << i << " ";)
   
    std::vector<std::pair<double, double>> distParams{par.sc};
    std::transform(par.sizeClasses.begin(),
            par.sizeClasses.end(),
            distParams.begin(),
            distParams.begin(),
            [](double& v, std::pair<double, double>& p) {
                    return std::make_pair(v * 0.92, v * 0.3 * 0.92);
            });

	D(for (auto& i : distParams) std::cout << i.first << " " << i.second << " / " ;)
	
      
    double kEC50 = 1. / (par.FG_EC50Orders + 1.);
    //----------------------------------------------------------------------------------------
    //Allocating space for results and results
    const size_t ts = par.length / par.step;
    const size_t array_size = par.FGnum * par.sc * ts;

    std::vector<ParB> parB_array;
    parB_array.reserve(array_size *  par.reiterations);
    
    //----------------------------------------------------------------------------------------
    //Work work work work work!
    #pragma omp parallel for
    for (int reit = 0; reit < par.reiterations; ++reit)
    {

    //----------------------------------------------------------------------------------------
    // Initializing population    
        Container C;
        Population* c_array = new Population[par.sc * par.FGnum];
        size_t cnt = 0;
        for (auto sc = 0; sc < par.sc; ++sc)
        {		
            for (auto fg = 0; fg < par.FGnum; ++fg,++cnt)
            {			
                c_array[cnt].s = sc;
                c_array[cnt].fg = fg;
                c_array[cnt].idx = cnt; //sc << 16 | fg;
                c_array[cnt].normParams = &distParams[sc];
                c_array[cnt].size = par.sizeClasses[sc];
                c_array[cnt].initPop = par.InitN[sc];
                c_array[cnt].probabilityOfDying = par.Pd_C[sc];

                c_array[cnt].fillUntil(par.InitN[sc], par.sizeClasses[sc]);
                c_array[cnt].EC50s = fg == 0 ? 4. + par.size_EC50Orders /
                        (double(par.sizeClasses.size()) - 1.) * double(sc) :
                        c_array[cnt-1].EC50s + (par.FG_EC50Orders / (par.FGnum - 1.));

                D(std::cout << c_array[cnt].EC50s << " ";)
                C.insert(&c_array[cnt]);
            }
                D(std::cout << std::endl;)
        }
        double val = (**(C.get<s_and_fg>().find(std::make_tuple(0,0)))).EC50s;
        for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
        {
            (**it).Growth_adj_tolerance = 1. - pow(kEC50 * ( (**it).EC50s - val ) , par.q); 
        }

       
        for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
        {
            (**it).initRandPop();
        }

    //----------------------------------------------------------------------------------------
    // Computing initial biomass and storing result
                
        double old_biomass = 0.;
        for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
        {
            ParB result((**it).s,
                        (**it).fg,
                        0,
                        reit
                        );      
            result.mass = (**it).getMass();
            result.num_alive = (**it).getNumAlive();
            old_biomass += result.mass;
            #pragma omp critical
            {
                parB_array.push_back(result);
            }            
            
            
        }

    //----------------------------------------------------------------------------------------
    // Do the evolution        
        for (auto i = 1; i < ts; ++i)
        {
            double L = i < par.t_stressor/par.step ? L1 : L1 + log10(pow(10,par.L2) * exp((-log(2) / par.tao ) * (i - par.t_stressor/par.step)));
            double new_biomass = 0.;
            for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
            {
                (**it).epsilon = (**it).Growth_adj_tolerance *
                        (1. - (1. / (1. + (pow(10, (**it).EC50s) / pow(pow(10, L), par.nlogit)))));

                (**it).growth = ((Unbiased_GR[(**   it).s] -
                        (old_biomass * par.BM_inhibit)) * par.step) * (**it).epsilon;

                if ((**it).growth < 0.) (**it).growth = 0.;

                (**it).doTheEvolution();

                //Storing results
                ParB result((**it).s,
                        (**it).fg,
                        i,
                        reit
                        );          
                result.mass = (**it).getMass();
                result.num_alive = (**it).getNumAlive();
                new_biomass += result.mass;
                #pragma omp critical
                {
                    parB_array.push_back(result);
                }
                
                
            }    
            old_biomass = new_biomass;
        }
       


    //----------------------------------------------------------------------------------------
    // Cleaning up
        delete [] c_array;        
    }

    auto stop = std::chrono::high_resolution_clock().now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds > (stop - start);

	
	std::cout << "Time taken by peace4U: " << duration.count() << " milliseconds " << std::endl;
	
    return parB_array;

}


boost::python::dict getAll(boost::python::dict pythonPar)
{   
    
    Parameters par = Parameters::setParameters(pythonPar);
    const size_t ts = par.length / par.step;

    auto parB_array = run(par);

    ParBContainer parB_container;
    for (auto &i : parB_array) parB_container.insert(&i);

    using dict = boost::python::dict;
    using list = boost::python::list;

    list reitList, tsList, massList, aliveList;


    for (auto i = 0; i < par.reiterations; ++i)
    {
        for (auto j = 0; j < ts; ++j)
        {
            double biomass = 0.;
            size_t alive =  0;
            for (auto it : boost::make_iterator_range(parB_container.get<t_and_reit>().equal_range(std::make_tuple(j,i))))
            {
                biomass += it->mass;
                alive += it ->num_alive;
            }
            reitList.append(i);
            tsList.append(j);
            massList.append(biomass);
            aliveList.append(alive);
            // std::cout << biomass << " " << alive << "\t";
        }
        // std::cout << std::endl << std::endl;
    }

    dict result;
    result["reit"] = reitList;
    result["ts"] = tsList;
    result["mass"] = massList;
    result["alive"] = aliveList;
    

	return result;

}


boost::python::dict getLimits(boost::python::dict pythonPar)
{   
    Parameters par = Parameters::setParameters(pythonPar);
    const size_t ts = par.length / par.step;
    auto parB_array = run(par);

    ParBContainer parB_container;
    for (auto &i : parB_array) parB_container.insert(&i);

    using dict = boost::python::dict;
    using list = boost::python::list;
    namespace ba = boost::accumulators;

    typedef ba::accumulator_set<double,ba::stats<ba::tag::extended_p_square> > accumulator_t;
    std::vector<double> p1 = {0.01, 0.05, 0.5, 0.95, 0.99};
    accumulator_t acc(ba::tag::extended_p_square::probabilities=p1);
    
    std::vector<std::string> headerP1 = {"0.01","0.05","0.5","0.95","0.99"};
    std::vector<accumulator_t> accBio, accAlive;
    for (size_t i = 0; i < ts; ++i)
    {
        accBio.push_back(acc);
        accAlive.push_back(acc);
    }
    std::vector<list> quantileMass{p1.size()}, quantileAlive{p1.size()};


    for (auto i = 0; i < par.reiterations; ++i)
    {
        for (auto j = 0; j < ts; ++j)
        {
            double biomass = 0.;
            size_t alive =  0;
            for (auto it : boost::make_iterator_range(parB_container.get<t_and_reit>().equal_range(std::make_tuple(j,i))))
            {
                biomass += it->mass;
                alive += it ->num_alive;
            }
            accBio[j](biomass);
            accAlive[j](boost::lexical_cast<double>(alive));
            // std::cout << biomass << " " << alive << "\t";
        }
        // std::cout << std::endl << std::endl;
    }

    dict bioDict,aliveDict;
    for (auto &i : accBio)
    {
        for (auto j = 0; j < p1.size(); ++j)
        {
            double dummy = ba::extended_p_square(i)[j];
            quantileMass[j].append(dummy);
        }        
    }
    for (auto &i : accAlive)
    {
        for (auto j = 0; j < p1.size(); ++j)
        {
            double dummy = ba::extended_p_square(i)[j];
            quantileAlive[j].append(dummy);
        }        
    }



    for (auto i = 0; i < p1.size(); ++i)
    {
        bioDict[headerP1[i]]=quantileMass[i];
        aliveDict[headerP1[i]]=quantileAlive[i];
    }   

    boost::python::dict result;
    result["biomass"] = bioDict;
    result["alive"] = aliveDict;




    return result;
}







