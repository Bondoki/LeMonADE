#ifndef RANDOMNUMBERGENERATORS_H
#define RANDOMNUMBERGENERATORS_H

#include <omp.h>
#include <chrono>
#include <random>
#include <iostream>
#include <stdint.h>

#include <omp.h>

class RandomNumberGenerators
{
public:
    typedef std::mt19937 Engine;
    typedef std::uniform_real_distribution<double> Distribution;

    RandomNumberGenerators() : engines(), distribution(0.0, 1.0)
    {
        int threads = std::max(1, omp_get_max_threads());
        for(int seed = 0; seed < threads; ++seed)
        {
            engines.push_back(Engine(seed));
        }
    }

    //! returns random unsignet 32 bit integer from R250Engine
    //inline uint32_t r250_rand32(){return r250Engine->r250_rand();} //range [0:2e31-1]
    //! returns random double from R250Engine
    //inline double r250_drand(){return distribution(engines[id]())} //range [0.0:1.0]

    inline uint32_t r250_rand32()
    {
    	int id = omp_get_thread_num();
    	return uint32_t(engines[id]()); //range [0:2e31-1]
    }

    inline double r250_drand()
    {
    	int id = omp_get_thread_num();
    	return distribution(engines[id]); //range [0:1]
    }

    std::vector<Engine> engines;
    Distribution distribution;
    
    void seedAll()
    {
        int threads = std::max(1, omp_get_max_threads());
        for(int id = 0; id < threads; ++id)
        {
            engines[id].seed(std::chrono::system_clock::now().time_since_epoch().count());
            
        }
       
    }
};

#endif

