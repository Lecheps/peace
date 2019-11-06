#========================================================#
# Build the Boost dependencies for the project using a specific version of python #
#========================================================#

set(BoostVersion 1.71.0)
# set(BoostSHA256 96b34f7468f26a141f6020efb813f1a2f3dfb9797ecf76a7d7cbd843cc95f5bd)

string(REGEX REPLACE "beta\\.([0-9])$" "beta\\1" BoostFolderName ${BoostVersion})
string(REPLACE "." "_" BoostFolderName ${BoostFolderName})
set(BoostFolderName boost_${BoostFolderName})

#Boost jam will use the Python version active in the current terminal
set (libApp  ${Python_VERSION_MAJOR}${Python_VERSION_MINOR})


ExternalProject_Add(Boost
    PREFIX Boost
    URL  https://dl.bintray.com/boostorg/release/${BoostVersion}/source/${BoostFolderName}.tar.gz
    # URL_HASH SHA256 ${BoostSHA256}
    CONFIGURE_COMMAND ./bootstrap.sh
                                      --with-libraries=python
                                      --with-python=${Python_EXECUTABLE}
    BUILD_COMMAND  ./b2 install
                                                           include=${Python_INCLUDE_DIRS}
                                                           variant=release
                                                           optimization=speed
                                                           link=static
                                                           cxxflags='-fPIC'
                                                           --prefix=${CMAKE_BINARY_DIR}/extern
                                                           -d 0
                                                           -j8
                                                        
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    )

set(Boost_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib/ )
set(Boost_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include/boost/)

if(${Python_VERSION_MAJOR} EQUAL 3) #Note the numpy need to be installed for the library to have been created
  message(STATUS "Using Python3")
  set(Boost_LIBRARIES -lboost_python${libApp} -lboost_numpy${libApp})
else()
  message(STATUS "Using Python2")
  set(Boost_LIBRARIES -lboost_python -lboost_numpy)
endif()
