#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <thread>
#include <atomic>


enum class DISTRIBUTION
{
    NORMAL,
    UNIFORM
}; 

class RandGenerator
{

private: 
    static std::atomic<uint32_t> seed;

public:

    typedef boost::mt11213b base_generator_type;
    //typedef boost::mt19937 base_generator_type;

    
    
    //Factory method
    
    
    virtual void set_parameters(double a, double b) = 0;
    
    virtual void start_async_gen() = 0;
    virtual void do_async_gen() = 0;
    virtual void stop_async_gen() = 0;
    virtual double get_from_async() = 0;
    base_generator_type base_generator;
    // std::atomic<bool> cancellation_token;

    boost::lockfree::spsc_queue<double, boost::lockfree::capacity<10000>> queue;


    static RandGenerator* make_generator(DISTRIBUTION);
    virtual double get() = 0;
    RandGenerator();
};


template<typename T>
class Distribution: public RandGenerator
{
private:
   
    // typedef T distribution_type;
    T distribution;
    typedef boost::variate_generator<base_generator_type&,T> derived_generator;
    derived_generator gen;  
    std::thread rand_thread;
    std::atomic<bool> cancellation_token; 
    void do_async_gen() override;


public:
    Distribution();
    double get() override
    {
        return gen();
    }

    void set_parameters(double a, double b) override;
    
    void start_async_gen() override;
    void stop_async_gen() override; 

    //WARNING: get_from_async() shouldn't be used at the same time as get(). Horror
    //might ensue since both functions share the not-thread-safe base generator
    double get_from_async() override; 

};