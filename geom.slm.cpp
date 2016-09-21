#include"geom.slm.hpp"

__GEOM_SLM_SPEC( float );
__GEOM_SLM_SPEC( double );
__GEOM_SLM_SPEC( long double );
#include"mpreal.use.h"
#ifdef  __MPREAL_USE
#include"mpreal.h"
__GEOM_SLM_SPEC( mpfr::mpreal );
#endif
