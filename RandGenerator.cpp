#include "RandGenerator.h"

std::atomic<uint32_t> RandGenerator::seed{41};

RandGenerator::RandGenerator()
{
    base_generator.seed(seed++);
}


RandGenerator* RandGenerator::make_generator(DISTRIBUTION select_distribution)
{
    switch (select_distribution)
    {
    case DISTRIBUTION::NORMAL:
        return new Distribution<boost::normal_distribution<double>>();      
        break;
    
    case DISTRIBUTION::UNIFORM:
        return new Distribution<boost::uniform_real<double>>();   
        break;
    default:
        return new Distribution<boost::uniform_real<double>>(); 
        break;
    }
}


template<typename T> //Need to specialize the class if the distribution constructor takes a different number of arguments
Distribution<T>::Distribution() : distribution(0.0,1.0), gen(base_generator,distribution)
{
    // for (auto i = 0; i < 10; ++i) std::cout <<gen() << " ";
    // std::cout << std::endl;
}

template<typename T>
void  Distribution<T>::set_parameters(double a, double b)
{
    T distribution(a,b);
    gen.distribution() =  distribution;
}

template<typename T>
void Distribution<T>::do_async_gen()
{
    while (cancellation_token) queue.push(get());		
}

template<typename T>
void Distribution<T>::start_async_gen()
{
    cancellation_token = true;
    rand_thread = std::thread(&Distribution<T>::do_async_gen,this);
}

template<typename T>
void Distribution<T>::stop_async_gen()
{
    cancellation_token = false;
    rand_thread.join();
}

template<typename T>
double Distribution<T>::get_from_async()
{
    double rand;
    while (!queue.pop(rand)) {};
    return rand;
}