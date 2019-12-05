#include "Parameters.h"


void Parameters::setDeathProbabilities()
    {
        Fi_1 = 0.3314 * exp((-1.001 * log(b + 5) - 0.0366) * g);
        Fi_1_Pd = Pd * exp((-1.001 * log(b_Pd + 5) - 0.0366) * g_Pd);
        Pd_C = sizeClasses;
        for (auto& i : Pd_C) i = Fi_1_Pd * pow(i + b_Pd, g_Pd);            
    }


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