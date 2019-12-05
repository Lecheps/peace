// peace.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



#include <iostream>
#include "peace_utilities.h"

int main()
{
    // auto start = std::chrono::high_resolution_clock::now();
    
    Parameters par;
    // par.reiterations = 30;  
    auto results = run(par,false);
    
    // auto stop = std::chrono::high_resolution_clock().now();
	// auto duration = std::chrono::duration_cast<std::chrono::milliseconds > (stop - start);
    // std::cout << "Time taken by peace4U: " << duration.count() << " milliseconds " << std::endl;
    
    return 0;

    
}





 










