#include"ixs.angular.dat.hpp"
IXS_ANGULAR_DATA_SPEC( float );
IXS_ANGULAR_DATA_SPEC( double );
IXS_ANGULAR_DATA_SPEC( long double );
#include"mpreal.use.h"
#if defined (__MPREAL_USE)
#include"mpreal.h"
IXS_ANGULAR_DATA_SPEC( mpfr::mpreal );
#endif
