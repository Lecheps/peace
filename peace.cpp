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

struct Parameters {
    static const size_t FGnum = 6;
    static const size_t sc =  6;
    std::array<double,sc> sizeClasses = { 5,10,25,50,100,200 };

    double Pd = 0.05;
    double g = -0.0035;
    double b = -4.5195;
    double Fi_1, Fi_1_Pd;
    double b_Pd = -4.0;
    double g_Pd = -0.00215;
    std::array<double,sc> Pd_C;

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
    std::array<size_t,FGnum> InitN = { 9000,750,200,130,20,10 };
    size_t reiterations = 3;

};

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

struct parB
{
    unsigned short int s = std::numeric_limits<unsigned short int>::max();
    unsigned short int fg = std::numeric_limits<unsigned short int>::max();
    size_t idx = std::numeric_limits<size_t>::max();
    size_t t = std::numeric_limits<size_t>::max();
    double mass = 0.;
    // size_t num_alive = 0;
    std::vector<size_t>parD;
};

using parBContainer = boost::multi_index_container
<
    parB*, // the data type stored
    boost::multi_index::indexed_by
        <
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct s_par>,
                boost::multi_index::member<parB, unsigned short int, &parB::s>
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct fg_par>,
                boost::multi_index::member<parB, unsigned short int, &parB::fg >
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct idx_par>,
                boost::multi_index::member<parB, size_t, &parB::idx >
            >,
            boost::multi_index::hashed_non_unique
            <
                boost::multi_index::tag<struct t_par>,
                boost::multi_index::member<parB, size_t, &parB::t >
            >,
            boost::multi_index::ordered_non_unique
            <
                boost::multi_index::tag<struct s_and_t>,
                boost::multi_index::composite_key
                <
                    parB*,
                    boost::multi_index::member<parB,unsigned short int,&parB::s>,
                    boost::multi_index::member<parB,size_t,&parB::t>
                >,
                boost::multi_index::composite_key_compare
                <
                    std::less<unsigned short int>,     
                    std::less<size_t> 
                >
            >
        >
>;


int main()
{
    
    //----------------------------------------------------------------------------------------
	//Starting timer
    auto start = std::chrono::high_resolution_clock::now();
   
    //----------------------------------------------------------------------------------------
    /*
    Initial calculations
    */
    Parameters par; //Structure containing the default parameter set
    
    double L1 = 1.0;
    double Fi_1 = 0.3314 * exp((-1.001 * log(par.b + 5) - 0.0366) * par.g);
    std::array<double,par.sc> Unbiased_GR;
    Unbiased_GR = par.sizeClasses;
    for (auto& i : Unbiased_GR ) i = Fi_1 * pow(i + par.b, par.g);
    par.setDeathProbabilities();

    D(for (auto& i : Unbiased_GR) std::cout << i << " ";)
   
    std::array<std::pair<double, double>,par.sc> distParams;
    std::transform(par.sizeClasses.begin(),
            par.sizeClasses.end(),
            distParams.begin(),
            distParams.begin(),
            [](double& v, std::pair<double, double>& p) {
                    return std::make_pair(v * 0.92, v * 0.3 * 0.92);
            });

	D(for (auto& i : distParams) std::cout << i.first << " " << i.second << " / " ;)
	
    Container C;
    Population* c_array = new Population[par.sc * par.FGnum];
    
    
    double kEC50 = 1. / (par.FG_EC50Orders + 1.);
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

    //----------------------------------------------------------------------------------------
    //Allocating space for temporary variables and results
    const size_t ts = par.length / par.step;
    const size_t array_size = par.FGnum * par.sc * ts;

    parB* parB_array = new parB[array_size];
    for (auto i = 0; i < array_size; ++i) parB_array[i].parD.resize(par.reiterations);
    parBContainer parB_container;
    
    // parB*** parB_array = new parB**[par.reiterations];
    // for (auto i = 0; i < par.reiterations; ++i)
    // {
    //     parB_array[i] = new parB*[ts];
    //     for (auto j =  0; j < ts; ++j) parB_array[i][j]=new parB[par.sc];
    // }
    // parBContainer parB_container;

      

    double*** FQ_sc = new double**[par.reiterations];
    for (auto i = 0; i < par.reiterations; ++i)
    {
        FQ_sc[i] = new double*[ts];
        for (auto j =  0; j < ts; ++j) FQ_sc[i][j]=new double[par.sc];
    }

    // for (auto i = 0; i < par.sc; ++i) FQ_sc[i] = new double[ts];
    double* B_fin = new double[ts];
    double** totD = new double* [par.reiterations];
    for (auto i = 0; i < par.reiterations; ++i) totD[i] = new double[ts];

    // //----------------------------------------------------------------------------------------
    //Work work work work work!
    #pragma omp parallel for
    for (int reit = 0; reit < par.reiterations; ++reit)
    {
        double* B = new double[ts];
        parB* parB_array = new parB[array_size];
        for (auto i = 0; i < array_size; ++i) parB_array[i].parD.resize(par.reiterations);
        parBContainer parB_container;

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
                
        size_t parB_cnt = 0;
        for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it, ++parB_cnt)
        {
            parB_array[parB_cnt].mass = (**it).getMass();
            parB_array[parB_cnt].parD[reit] = (**it).getNumAlive();
            parB_array[parB_cnt].t = 0;
            // if (reit == 100)
            // {
                parB_array[parB_cnt].s = (unsigned short int) (**it).s; 
                parB_array[parB_cnt].fg = (unsigned short int) (**it).fg;
                parB_array[parB_cnt].idx = (**it).idx;
                parB_container.insert(&(parB_array[parB_cnt]));
            // } 
        }

        // Getting total number of organisms that are alive
        size_t tot_alive = 0;
        for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it) tot_alive += (**it).getNumAlive(); 

        
        for (auto i = 0; i < par.sc; ++i)
        {
            double fq_sc = 0;
            for (auto it : boost::make_iterator_range(parB_container.get<s_and_t>().equal_range( std::make_tuple(i,0) )))
            {
                fq_sc += it->parD[reit];
            }
            FQ_sc[reit][0][i] =  fq_sc / tot_alive;
            D(std::cout << "Size class " << i << " FQ_sc " << FQ_sc[reit][0][i] << std::endl;)

        }


        B[0] = 0.;
        totD[reit][0] = 0.;
        for (auto it : boost::make_iterator_range(parB_container.get<t_par>().equal_range(0)))
        {
                B[0] += it->mass;
                totD[reit][0] += it->parD[reit];
        }


        //Do the evolution baby!

        for (auto i = 1; i < ts; ++i)
        {
            // std::cout << i << " ";
            double L = i < par.t_stressor/par.step ? L1 : L1 + log10(pow(10,par.L2) * exp((-log(2) / par.tao ) * (i - par.t_stressor/par.step)));
            
            // for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
            // {
            //     (**it).epsilon = (**it).Growth_adj_tolerance *
            //             (1. - (1. / (1. + (pow(10, (**it).EC50s) / pow(pow(10, L), par.nlogit)))));

            //     (**it).growth = ((Unbiased_GR[(**   it).s] -
            //             (B[i - 1] * par.BM_inhibit)) * par.step) * (**it).epsilon;

            //     if ((**it).growth < 0.) (**it).growth = 0.;

            //     // std::cout << (**it).s << " " << (**it).fg << std::endl;

            //     (**it).doTheEvolution();
            // }

            // #pragma omp parallel for
            for (int k = 0; k< 36 ; ++k)
            {
                c_array[k].epsilon = c_array[k].Growth_adj_tolerance *
                        (1. - (1. / (1. + (pow(10,c_array[k].EC50s) / pow(pow(10, L), par.nlogit)))));

                c_array[k].growth = ((Unbiased_GR[c_array[k].s] -
                        (B[i - 1] * par.BM_inhibit)) * par.step) * c_array[k].epsilon;

                if (c_array[k].growth < 0.) c_array[k].growth = 0.;

                // std::cout << c_array[k].s << " " << c_array[k].fg << std::endl;
                c_array[k].doTheEvolution();
            }




            for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it, ++parB_cnt)
            {
                
                parB_array[parB_cnt].mass = (**it).getMass();
                parB_array[parB_cnt].parD[reit] = (**it).getNumAlive();
                parB_array[parB_cnt].t = i;
                // if(reit == 100)
                // {
                    parB_array[parB_cnt].s = (unsigned short int) (**it).s; 
                    parB_array[parB_cnt].fg = (unsigned short int) (**it).fg;
                    parB_array[parB_cnt].idx = (**it).idx;
                    parB_container.insert(&(parB_array[parB_cnt]));
                // } 
            }

            size_t tot_alive = 0;
            for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it) tot_alive += (**it).getNumAlive(); 

            
            for (auto j = 0; j < par.sc; ++j)
            {
                double fq_sc = 0;
                for (auto it : boost::make_iterator_range(parB_container.get<s_and_t>().equal_range( std::make_tuple(j,i) )))
                {
                    fq_sc += it->parD[reit];
                }
                FQ_sc[reit][i][j] =  fq_sc / tot_alive;
                D(std::cout << "Size class " << j << " FQ_sc " << FQ_sc[reit][i][j] << std::endl;)
            }

            B[i] = 0.;
            totD[reit][i] = 0.;
            for (auto it : boost::make_iterator_range(parB_container.get<t_par>().equal_range(i)))
            {
                    B[i] += it->mass;
                    totD[reit][i] += it->parD[reit];
            }
            // std::cout << std::fixed << std::setprecision(2) << B[i] << " " << totD[reit][i] << "\t";







        }
        // std::cout << std::endl << std::endl << std::endl;
        
        delete [] parB_array;
        delete [] c_array;
        for(auto i = 0; i<ts;++i) B_fin[i] = B[i];
        delete [] B;

	}
	
    for(auto i = 0; i<ts;++i) std::cout << B_fin[i] << " ";
    std::cout << std::endl;

	auto stop = std::chrono::high_resolution_clock().now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds > (stop - start);

	
	//Freeing memory
    delete [] parB_array;
    delete [] c_array;
    
    for (auto i = 0; i < par.reiterations; ++i)
    {
        for (auto j = 0; j < ts; ++j) delete[] FQ_sc[i][j];
        delete[] FQ_sc[i];
    }
    delete FQ_sc;


    delete [] B_fin;	
    for (auto i = 0; i < par.reiterations; ++i) delete[] totD[i];
    delete [] totD;


	std::cout << "Time taken by peace4U: " << duration.count() << " milliseconds " << std::endl;
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file


