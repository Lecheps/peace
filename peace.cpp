// peace.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#ifdef DEBUG_ENABLE
#define D(x) x
#else 
#define D(x)
#endif

#include <iostream>
#include "RandGenerator.h"
// #include <boost/generator_iterator.hpp>
#include <boost/range/combine.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
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
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numpy.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>
#include <boost/lexical_cast.hpp>

struct Parameters {
    size_t FGnum;
    size_t sc;
    std::vector<double> sizeClasses;

    double Pd = 0.05;
    double g = -0.0035;
    double b = -4.5195;
    double Fi_1, Fi_1_Pd;
    double b_Pd = -4.0;
    double g_Pd = -0.00215;
    std::vector<double> Pd_C;

    void setDeathProbabilities()
    {
        Fi_1 = 0.3314 * exp((-1.001 * log(b + 5) - 0.0366) * g);
        Fi_1_Pd = 0.3314 * exp((-1.001 * log(b + 5) - 0.0366) * g_Pd);
        Pd_C = sizeClasses;
        for (auto& i : Pd_C) i = Fi_1_Pd * pow(i + b_Pd, g_Pd);            
    }

    double BM_inhibit = 2e-6;
    size_t length = 160;
    size_t step = 1;
    double t_stressor = 60;
    double L2 = 1.0;
    double tao = 7.0;
    double nlogit = 2.0;
    size_t size_EC50Orders = 2;
    size_t FG_EC50Orders = 2;
    double q = 3.;
    std::vector<size_t> InitN;
    size_t reiterations = 3;

    Parameters();
    Parameters(size_t,size_t);
    // Parameters(int FG, int SC) : FGnum(FG), sc(SC), sizeClasses(std::array<double, SC>{}, initN(std::array<size_t,FG>{}) {};

    boost::python::dict getParameters();
    static Parameters setParameters(boost::python::dict);
    

};


 
 Parameters::Parameters() : FGnum(6), sc(6), sizeClasses({ 5,10,25,50,100,200 }), InitN({ 9000,750,200,130,20,10 }) 
 {    
 };

 Parameters::Parameters(size_t FG, size_t SC) : FGnum(FG), sc(SC) 
 {    
    sizeClasses.resize(SC);
    InitN.resize(FG);
 };


boost::python::dict Parameters::getParameters()
    {
        
        boost::python::dict result;
        boost::python::list sizeList, initList;
        for (auto& i : sizeClasses) sizeList.append(i);
        for (auto& i : InitN) initList.append(i);
        result["FGnum"]             = FGnum;
        result["sc"]                = sc;
        result["Pd"]                = Pd;
        result["g"]                 = g;
        result["b"]                 = b;
        result["Fi_1"]              = Fi_1;
        result["Fi_1_pd"]           = Fi_1_Pd;
        result["b_Pd"]              = b_Pd;
        result["g_Pd"]              = g_Pd;
        result["BM_inhibit"]        = BM_inhibit;
        result["length"]            = length;
        result["step"]              = step;
        result["t_stressor"]        = t_stressor;
        result["L2"]                = L2;
        result["tao"]               = tao;
        result["nlogit"]            = nlogit;
        result["size_EC50Orders"]   = size_EC50Orders;
        result["FG_EC50Orders"]     = FG_EC50Orders;
        result["q"]                 = q;
        result["reiterations"]      = reiterations;
        result["Pd"]                = Pd;
        result["sizeClasses"]       = sizeList;
        result["InitN"]             = initList;

        return result;
    }

Parameters Parameters::setParameters(boost::python::dict par)
{
    using get_size_t = boost::python::extract<size_t>; 
    using get_double = boost::python::extract<double>;
    using list       = boost::python::list;
    using get_list   = boost::python::extract<list>;
    
    Parameters result;
    result = Parameters(get_size_t(par["FGnum"]),get_size_t(par["sc"]));
    result.Pd                   = get_double(par["Pd"]);
    result.g                    = get_double(par["g"]);
    result.b                    = get_double(par["b"]);
    result.Fi_1                 = get_double(par["Fi_1"]);
    result.Fi_1_Pd              = get_double(par["Fi_1_pd"]);
    result.b_Pd                 = get_double(par["b_Pd"]);
    result.g_Pd                 = get_double(par["g_Pd"]);
    result.BM_inhibit           = get_double(par["BM_inhibit"]);
    result.length               = get_size_t(par["length"]);
    result.step                 = get_size_t(par["step"]);
    result.t_stressor           = get_double(par["t_stressor"]);
    result.L2                   = get_double(par["L2"]);
    result.tao                  = get_double(par["tao"]);
    result.nlogit               = get_double(par["nlogit"]);
    result.size_EC50Orders      = get_size_t(par["size_EC50Orders"]);
    result.FG_EC50Orders        = get_size_t(par["FG_EC50Orders"]);
    result.q                    = get_double(par["q"]);
    result.reiterations         = get_size_t(par["reiterations"]);
    result.Pd                   = get_double(par["Pd"]);
    list InitN                  = get_list(par["InitN"]);
    list sizeClasses            = get_list(par["sizeClasses"]);
    size_t cnt = 0;
    for(auto it = result.InitN.begin(); it < result.InitN.end(); ++it, ++cnt)
    {
        *it = get_size_t(InitN[cnt]);
    } 
    cnt = 0;
    for(auto it = result.sizeClasses.begin(); it < result.sizeClasses.end(); ++it, ++cnt)
    {
        *it = get_double(sizeClasses[cnt]);
    }
    return result;
}


struct Population
{

    Population() :
    normal(RandGenerator::make_generator(DISTRIBUTION::NORMAL)),
    uniform(RandGenerator::make_generator(DISTRIBUTION::UNIFORM))
    { 
            normParams = nullptr; 
            // normal->start_async_gen();
            // uniform->start_async_gen();
    }

    // Population(bool dummy) : normal(nullptr),uniform(nullptr)
    // { 
    //         normParams = nullptr; 
    //         // normal->start_async_gen();
    //         // uniform->start_async_gen();
    // }

    ~Population()
    {
        // normal->stop_async_gen();
        // uniform->stop_async_gen();
    }

    Population(const Population &p) : 
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
    
    static const size_t potsize = 10000;
    
    // std::array<double,potsize> c;
    std::list<double> c;

    void fillUntil(size_t pos, double& val)
    {
        // c.clear();
        // for (size_t i = 0; i < pos; ++i)
        // {
        //         // c[i] =val;
        //         c.push_back(val);
        // }
        std::list<double> dummy(pos,val);
        c = std::move(dummy);
    }

    void initRandPop()
    {
        std::list<double> dummy(initPop,size);
        for (auto& i :  dummy) i = i * 2. / 3. * uniform->get() + i / 3.;
        c = std::move(dummy);
        // c.clear();
        // for (auto i = 0; i < initPop; ++i) c.push_back(size * 2. / 3. * genie->get_uni() + size / 3.);
    }

    double getMass()
    {
        double result = 0.;
        for (const auto& i : c) result +=i;
        return result;
    }

    size_t getNumAlive()
    {
        // size_t result = 0;
        // for (const auto& i : c) if (i != 0) result++;
        // return result;
        return c.size();
    }

    void doTheEvolution()
    {                
        for (auto it = c.begin(); it != c.end(); ++it)
        {
            *it += normParams->first / 2. * growth;
            if ( *it > (normParams->second  * normal->get() + normParams->first))
            {
                *it /= 2.;
                if (c.size() < potsize && uniform->get() >= 0.05)
                {
                    c.insert(it,*it);
                } 
            }     
            
            if (uniform->get() < 0.05)
            {
                it = c.erase(it);
                it--;
            }
        }
    }
  
    // void doTheEvolution_o()
    // {                
    //     for (auto it = c.begin(); it != c.end(); ++it)
    //     {
    //         *it += normParams->first / 2. * growth;
    //         if (*it > (normParams->second  * genie->get_rand_norm() + normParams->first))
    //         {
    //             *it /= 2.;
    //             if (c.size() < potsize)
    //             {
    //                 c.insert(it,*it);
    //             } 
    //         }            
    //     }
        
    //     for (auto it = c.begin(); it != c.end();)
    //     {
    //         if (genie->get_rand_uni() < 0.05)
    //         {
    //             it = c.erase(it);
    //         } 
    //         else ++it;
    //     }
    // }
};

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




// boost::python::dict runFromPython(boost::python::dict pythonPar)
// {
//     auto start = std::chrono::high_resolution_clock::now();
//     Parameters par = Parameters::setParameters(pythonPar);
//     double L1 = 1.0;
//     par.setDeathProbabilities();
//     std::vector<double> Unbiased_GR = par.sizeClasses;
//     for (auto& i : Unbiased_GR ) i = par.Fi_1 * pow(i + par.b, par.g);
//     par.setDeathProbabilities();

//     D(for (auto& i : Unbiased_GR) std::cout << i << " ";)
   
//     std::vector<std::pair<double, double>> distParams{par.sc};
//     std::transform(par.sizeClasses.begin(),
//             par.sizeClasses.end(),
//             distParams.begin(),
//             distParams.begin(),
//             [](double& v, std::pair<double, double>& p) {
//                     return std::make_pair(v * 0.92, v * 0.3 * 0.92);
//             });

// 	D(for (auto& i : distParams) std::cout << i.first << " " << i.second << " / " ;)
	
      
//     double kEC50 = 1. / (par.FG_EC50Orders + 1.);
//     //----------------------------------------------------------------------------------------
//     //Allocating space for results and results
//     const size_t ts = par.length / par.step;
//     const size_t array_size = par.FGnum * par.sc * ts;

//     std::vector<ParB> parB_array;
//     parB_array.reserve(array_size *  par.reiterations);
    
//     //----------------------------------------------------------------------------------------
//     //Work work work work work!
//     #pragma omp parallel for
//     for (int reit = 0; reit < par.reiterations; ++reit)
//     {

//     //----------------------------------------------------------------------------------------
//     // Initializing population    
//         Container C;
//         Population* c_array = new Population[par.sc * par.FGnum];
//         size_t cnt = 0;
//         for (auto sc = 0; sc < par.sc; ++sc)
//         {		
//             for (auto fg = 0; fg < par.FGnum; ++fg,++cnt)
//             {			
//                 c_array[cnt].s = sc;
//                 c_array[cnt].fg = fg;
//                 c_array[cnt].idx = cnt; //sc << 16 | fg;
//                 c_array[cnt].normParams = &distParams[sc];
//                 c_array[cnt].size = par.sizeClasses[sc];
//                 c_array[cnt].initPop = par.InitN[sc];

//                 c_array[cnt].fillUntil(par.InitN[sc], par.sizeClasses[sc]);
//                 c_array[cnt].EC50s = fg == 0 ? 4. + par.size_EC50Orders /
//                         (double(par.sizeClasses.size()) - 1.) * double(sc) :
//                         c_array[cnt-1].EC50s + (par.FG_EC50Orders / (par.FGnum - 1.));

//                 D(std::cout << c_array[cnt].EC50s << " ";)
//                 C.insert(&c_array[cnt]);
//             }
//                 D(std::cout << std::endl;)
//         }
//         double val = (**(C.get<s_and_fg>().find(std::make_tuple(0,0)))).EC50s;
//         for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
//         {
//             (**it).Growth_adj_tolerance = 1. - pow(kEC50 * ( (**it).EC50s - val ) , par.q); 
//         }

       
//         for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
//         {
//             (**it).initRandPop();
//         }

//     //----------------------------------------------------------------------------------------
//     // Computing initial biomass and storing result
                
//         double old_biomass = 0.;
//         for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
//         {
//             ParB result((**it).s,
//                         (**it).fg,
//                         0,
//                         reit
//                         );      
//             result.mass = (**it).getMass();
//             result.num_alive = (**it).getNumAlive();
//             old_biomass += result.mass;
//             #pragma omp critical
//             {
//                 parB_array.push_back(result);
//             }            
            
            
//         }

//     //----------------------------------------------------------------------------------------
//     // Do the evolution        
//         for (auto i = 1; i < ts; ++i)
//         {
//             double L = i < par.t_stressor/par.step ? L1 : L1 + log10(pow(10,par.L2) * exp((-log(2) / par.tao ) * (i - par.t_stressor/par.step)));
//             double new_biomass = 0.;
//             for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
//             {
//                 (**it).epsilon = (**it).Growth_adj_tolerance *
//                         (1. - (1. / (1. + (pow(10, (**it).EC50s) / pow(pow(10, L), par.nlogit)))));

//                 (**it).growth = ((Unbiased_GR[(**   it).s] -
//                         (old_biomass * par.BM_inhibit)) * par.step) * (**it).epsilon;

//                 if ((**it).growth < 0.) (**it).growth = 0.;

//                 (**it).doTheEvolution();

//                 //Storing results
//                 ParB result((**it).s,
//                         (**it).fg,
//                         i,
//                         reit
//                         );          
//                 result.mass = (**it).getMass();
//                 result.num_alive = (**it).getNumAlive();
//                 new_biomass += result.mass;
//                 #pragma omp critical
//                 {
//                     parB_array.push_back(result);
//                 }
                
                
//             }    
//             old_biomass = new_biomass;
//         }
       


//     //----------------------------------------------------------------------------------------
//     // Cleaning up
//         delete [] c_array;        
//     }
	
//     //----------------------------------------------------------------------------------------
//     // Showing results
//     ParBContainer parB_container;
//     for (auto &i : parB_array) parB_container.insert(&i);

//     using dict = boost::python::dict;
//     using list = boost::python::list;

//     list reitList, tsList, massList, aliveList;


//     for (auto i = 0; i < par.reiterations; ++i)
//     {
//         for (auto j = 0; j < ts; ++j)
//         {
//             double biomass = 0.;
//             size_t alive =  0;
//             for (auto it : boost::make_iterator_range(parB_container.get<t_and_reit>().equal_range(std::make_tuple(j,i))))
//             {
//                 biomass += it->mass;
//                 alive += it ->num_alive;
//             }
//             reitList.append(i);
//             tsList.append(j);
//             massList.append(biomass);
//             aliveList.append(alive);
//             // std::cout << biomass << " " << alive << "\t";
//         }
//         // std::cout << std::endl << std::endl;
//     }

//     dict result;
//     result["reit"] = reitList;
//     result["ts"] = tsList;
//     result["mass"] = massList;
//     result["alive"] = aliveList;
    
// 	auto stop = std::chrono::high_resolution_clock().now();
// 	auto duration = std::chrono::duration_cast<std::chrono::milliseconds > (stop - start);

	
// 	std::cout << "Time taken by peace4U: " << duration.count() << " milliseconds " << std::endl;
// 	return result;
// }



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



std::vector<ParB> run(Parameters par)
{
    ScopedGILRelease scoped; 
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





BOOST_PYTHON_MODULE(peace)
{
    using namespace boost::python;
    class_< Parameters>("Parameters",init<>())
        .def(init<size_t,size_t>())
        .def("getParameters",&Parameters::getParameters)
        .def("setParameters",&Parameters::setParameters)
    ;

    def("run",getAll);
    def("runLimits",getLimits);

}

