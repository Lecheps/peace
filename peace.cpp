// peace.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#ifdef DEBUG_ENABLE
#define D(x) x
#else 
#define D(x)
#endif

#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/range/combine.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/key_extractors.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/lockfree/queue.hpp>
#include <math.h>
#include <array>
#include <chrono>
#include <limits>
#include <thread>
#include <atomic>
#include <vector>
#include <list>
#include <iomanip> 

typedef boost::mt11213b base_generator_type;
//typedef boost::mt19937 base_generator_type;

typedef boost::normal_distribution<> norm_dist;
typedef boost::uniform_real<> uni_dist;
typedef boost::variate_generator<base_generator_type&, norm_dist> norm_gen;
typedef boost::variate_generator<base_generator_type&, uni_dist> uni_gen;
const size_t randBuffer = 10000000;

struct randGenerator{

	static base_generator_type gen_norm, gen_uni;
	static norm_dist norm;
	static uni_dist uni;
	static norm_gen get_rand_norm;
	static uni_gen get_rand_uni;
	static std::thread thread_uni, thread_norm;
    static std::atomic<bool> cancellation_token_uni, cancellation_token_norm;
	 	
	static boost::lockfree::spsc_queue<double, boost::lockfree::capacity<randBuffer>> queue_norm, queue_uni;

	static void async_norm()
	{
            while (cancellation_token_norm) queue_norm.push(get_rand_norm());		
	}

	static void async_uni()
	{
            while (cancellation_token_uni) queue_uni.push(get_rand_uni());		
	}
        
    static double get_uni(){
        double uni;
        while (!queue_uni.pop(uni)) {};
        return uni;
    }
    
    static double get_norm(){
        double norm;
        while (!queue_norm.pop(norm)) {};
        return norm;
    }
};

base_generator_type randGenerator::gen_norm(42), randGenerator::gen_uni(84);
norm_dist randGenerator::norm(0.0, 1.0);
uni_dist randGenerator::uni(0.0, 1.0);
uni_gen randGenerator::get_rand_uni(gen_uni, uni);
norm_gen  randGenerator::get_rand_norm(gen_norm, norm);// , randGenerator::get_rand_norm1(gen1, norm1);
boost::lockfree::spsc_queue<double, boost::lockfree::capacity<randBuffer>> randGenerator::queue_norm;
boost::lockfree::spsc_queue<double, boost::lockfree::capacity<randBuffer>> randGenerator::queue_uni;
std::thread randGenerator::thread_norm, randGenerator::thread_uni;
std::atomic<bool> randGenerator::cancellation_token_uni, randGenerator::cancellation_token_norm;

struct Parameters {
    static const size_t FGnum = 6;
    static const size_t sc =  6;
    std::array<double,sc> sizeClasses = { 5,10,25,50,100,200 };

    double Pd = 0.05;
    double g = -0.0035;
    double b = - 4.5195;
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

    Population() { 
            normParams = nullptr; 
    };
    unsigned short int s;
    unsigned short int fg;
    size_t idx;

    double EC50s;
    double Growth_adj_tolerance;
    std::pair<double, double>* normParams;
    double size;
    size_t initPop;
    size_t ts_length;
    double epsilon;
    double growth;
    
    static const size_t potsize = 10000;
    
    std::array<double,potsize> c;
    // std::list<double> c;

    randGenerator* genie;
    static size_t cnt;


    void fillUntil(const size_t pos, const double& val)
    {
        // std::list<double> dummy(pos,val);
        // c = std::move(dummy);
        c.fill(0.);
        for (auto i = 0; i < pos; ++i) c[i] = val;
    }

    void initRandPop()
    {
        // std::list<double> dummy(initPop,size);
        // for (auto& i :  dummy) i = i * 2. / 3. * genie->get_rand_uni() + i / 3.;
        // c = std::move(dummy);
        c.fill(0.);
        for (auto i = 0; i < initPop; ++i) c[i] = size * 2. / 3. * genie->get_rand_uni() + size / 3.;        
    }

    double getMass()
    {
        double result = 0.;
        for (const auto& i : c) result +=i;
        return result;
    }

    size_t getNumAlive()
    {
        size_t result = 0;
        for (const auto& i : c) if (i != 0) result++;
        return result;
        // return c.size();
    }

    void doTheEvolution()
    {            
        std::vector<double> to_add;    
        for (auto it = c.begin(); it != c.end(); ++it)
        {
            *it += normParams->first / 2. * growth;
            if (*it > (normParams->second  * genie->get_rand_norm() + normParams->first))
            {
                *it /= 2.;
                // if (c.size() < potsize) c.insert(it,*it);
                to_add.push_back(*it);
            }            
        }
        auto it = to_add.begin();
        auto cit = std::find(c.begin(),c.end(),0.);
        while (it != to_add.end() && cit != c.end())
        {
            *cit = *it;
            ++it;
            cit = std::find(cit,c.end(),0.);
        }
        
        for (auto it = c.begin(); it != c.end();++it)
        {
            if (genie->get_rand_uni() < 0.05) *it = 0.;//it = c.erase(it);
            // else ++it;
        }
    }
};

size_t Population::cnt = 0;

using Container = boost::multi_index_container<
    Population*, // the data type stored
    boost::multi_index::indexed_by<
    boost::multi_index::hashed_non_unique<  
            boost::multi_index::tag<struct s>, 
            boost::multi_index::member<Population, unsigned short int, & Population::s> 
    >,
    boost::multi_index::hashed_non_unique<  
            boost::multi_index::tag<struct fg>, 
            boost::multi_index::member<Population, unsigned short int, & Population::fg >
    >,
    boost::multi_index::hashed_unique<
    boost::multi_index::tag<struct idx>,
    boost::multi_index::member<Population, size_t, & Population::idx >
    >
    >
>;

struct parB
{
    unsigned short int s;
    unsigned short int fg;
    size_t idx;
    size_t t;
    double mass;
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
    // for (auto i = 0; i < 100000000; ++i) double dummy = randGenerator::get_rand_norm();
    //----------------------------------------------------------------------------------------
    /*
    The performance bottleneck in this code is the random number generation. 
    /We attempt to mitigate that by generating them (uniform and normal distribution)
    in separate threads
    */
    D(std::cout << "Gonna start" << std::endl;)
    randGenerator::cancellation_token_norm = true;
    randGenerator::cancellation_token_uni = true;
    
    D(std::cout << "Starting threads..." << std::endl;)
    // randGenerator::thread_norm = std::thread(&randGenerator::async_norm);
    // randGenerator::thread_uni = std::thread(&randGenerator::async_uni);
    D(std::cout << "Threads away!" << std::endl;) 
    /*
    Note to self: when the threads are detached instead, the system takes a long
    time to liberate the resources. So long that if several instances of the program
    are launched after each other it bogs down the system. Why?
    */
    // randGenerator::thread_norm.detach();  
    // randGenerator::thread_uni.detach();

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
                    return std::make_pair(v * 0.92, v * 0.92 * 0.08);
            });

	D(for (auto& i : distParams) std::cout << i.first << " " << i.second << " / " ;)
	
    Container C;
    Population* c_array = new Population[par.sc * par.FGnum];
    auto& idxByS = C.get<s>();
    double kEC50 = 1. / (par.FG_EC50Orders + 1.);

    size_t cnt = 0;
    Population* head = c_array;
    for (auto sc = 0; sc < par.sc; ++sc)
    {		
        for (auto fg = 0; fg < par.FGnum; ++fg,++cnt)
        {			
            c_array[cnt].s = sc;
            c_array[cnt].fg = fg;
            c_array[cnt].idx = sc << 16 | fg;
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

    for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
    {
            head = *it;
            size_t idx = head->s << 16 | 0;
            auto it_first = C.get<2>().find(idx);
            double val = (*it_first)->EC50s;
            head->Growth_adj_tolerance = 1. -
                    pow(kEC50 * (val - c_array->EC50s), par.q);
    	D(std::cout << 
                    head->Growth_adj_tolerance << " " <<
                    c_array->EC50s << " " <<
                    val << " " <<
                    std::endl;)
    }

   

    //----------------------------------------------------------------------------------------
    //Allocating space for temporary variables and results
    const size_t ts = par.length / par.step;
    const size_t array_size = par.FGnum * par.sc * ts;

    parB* parB_array = new parB[array_size];
    for (auto i = 0; i < array_size; ++i) parB_array[i].parD.resize(par.reiterations);
    parBContainer parB_container;
    auto& idxByS_parB = parB_container.get<s_par>();
    

    double*** FQ_sc = new double**[par.reiterations];
    for (auto i = 0; i < par.reiterations; ++i)
    {
        FQ_sc[i] = new double*[ts];
        for (auto j =  0; j < ts; ++j) FQ_sc[i][j]=new double[par.sc];
    }

    // for (auto i = 0; i < par.sc; ++i) FQ_sc[i] = new double[ts];
    double* B = new double[ts];
    double** totD = new double* [par.reiterations];
    for (auto i = 0; i < par.reiterations; ++i) totD[i] = new double[ts];

    //----------------------------------------------------------------------------------------
    //Work work work work work!

    for (auto reit = 0; reit != par.reiterations; ++reit)
    {
        // parB_container.clear();
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
            if (reit == 0)
            {
                parB_array[parB_cnt].s = (unsigned short int) (**it).s; 
                parB_array[parB_cnt].fg = (unsigned short int) (**it).fg;
                parB_array[parB_cnt].idx = (**it).idx;
                parB_container.insert(&(parB_array[parB_cnt]));
            } 
        }

        // for (auto& i parB_container.get<IndexBys_parB()>().begin();

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
            
            for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it)
            {
                (**it).epsilon = (**it).Growth_adj_tolerance *
                        (1. - (1. / (1. + (pow(10, (**it).EC50s) / pow(pow(10, L), par.nlogit)))));

                (**it).growth = ((Unbiased_GR[(**it).s] -
                        (B[i - 1] * par.BM_inhibit)) * par.step) * (**it).epsilon;

                if ((**it).growth < 0.) (**it).growth = 0.;
                
                (**it).doTheEvolution();
            }

            for (auto it = C.get<idx>().begin(); it != C.get<idx>().end(); ++it, ++parB_cnt)
            {
                
                parB_array[parB_cnt].mass = (**it).getMass();
                parB_array[parB_cnt].parD[reit] = (**it).getNumAlive();
                parB_array[parB_cnt].t = i;
                if(reit == 0)
                {
                    parB_array[parB_cnt].s = (unsigned short int) (**it).s; 
                    parB_array[parB_cnt].fg = (unsigned short int) (**it).fg;
                    parB_array[parB_cnt].idx = (**it).idx;
                    parB_container.insert(&(parB_array[parB_cnt]));
                } 
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
            std::cout << std::fixed << std::setprecision(2) << B[i] << " " << totD[reit][i] << "\t";







        }
        std::cout << std::endl << std::endl << std::endl;

	}
	


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


	delete [] B;	
	for (auto i = 0; i < par.reiterations; ++i) delete[] totD[i];
	delete [] totD;

    D(std::cout << "Joining random number generation threads..." << std::endl;)
    randGenerator::cancellation_token_norm = false;
    randGenerator::cancellation_token_uni = false;
    // randGenerator::thread_norm.join();
    // randGenerator::thread_uni.join();
    D(std::cout << "Done!" << std::endl;)
        
	std::cout << "Time taken by peace4U: " << duration.count() << " milliseconds " << Population::cnt << std::endl;
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


