#pragma once

#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/numpy.hpp>


class Parameters {
public:
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

    void setDeathProbabilities();

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



