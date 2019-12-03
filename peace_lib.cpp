#include <iostream>
#include "peace_utilities.h"



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