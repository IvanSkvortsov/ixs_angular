#ifndef __GEOM_SLM_H__
#define __GEOM_SLM_H__
#include<vector>
#include"mapping.t.h"
#include"matrix.slm.h"
#include"alpha.siz.h"

#define __GEOM_SLM_DEBUG
#ifdef  __GEOM_SLM_DEBUG
  #include<cassert>
#endif

// USAGE Interface
#define __GEOM_SLM_CX_INTERFACE( CX )\
inline reference CX##_x(int i){ return this->_##CX.x[i];}\
inline reference CX##_y(int i){ return  this->_##CX.y[i];}\
inline reference CX##_z(int i){ return this->_##CX.z[i];}\
inline const_reference CX##_x(int i)const{ return this->_##CX.x[i];}\
inline const_reference CX##_y(int i)const{ return this->_##CX.y[i];}\
inline const_reference CX##_z(int i)const{ return this->_##CX.z[i];}

#define __GEOM_SLM_SLMKX_INTERFACE( slm_kX )\
inline void slm_kX##_set_lx( int lx){ this->_##slm_kX##_it_lx = this->_##slm_kX + lx * lx + lx;}\
inline reference slm_kX(int m){ return this->_##slm_kX##_it_lx[m];}\
inline const_reference slm_kX(int m)const{ return this->_##slm_kX##_it_lx[m];}\
inline reference slm_kX(int lx, int m){ return this->_##slm_kX[ lx * lx + lx + m];}\
inline const_reference slm_kX(int lx, int m)const{ return this->_##slm_kX[ lx * lx + lx + m];}

template<typename T>
struct geom_slm : public std::vector<T>, public mapping_struct
{
public:
	typedef std::vector<T> base_class;
	typedef typename base_class::value_type value_type;
	typedef typename base_class::reference reference;
	typedef typename base_class::const_reference const_reference;
	typedef typename base_class::pointer pointer;
	typedef typename base_class::const_pointer const_pointer;
	typedef typename base_class::size_type size_type;

	typedef typename alpha_siz::_lmax_struct _lmax_struct;
	typedef struct
	{
		pointer x, y, z;
	} _3d_powers;
	typedef struct
	{
		int _la_size, _lb_size, _kA_lsize, _kB_lsize;
	} _lsize_struct;

protected:
	_lmax_struct _lmax;
	_3d_powers _CA, _CB;
	pointer _slm_kA, _slm_kB;
	_lsize_struct _lsize;

	void init_slmkX( T * __slm_kX, T * CX_x, T * CX_y, T * CX_z, const int lx_max, const int kx_lmax, T const * CX, matrix_slm<T> const & mx_slm);
	inline void error( const char * _method, const char _message[] = "nothing to do here" )const
	{ std::cerr << "Error: [" << this << "] geom_slm<T>::" << _method << ", " << _message << std::endl; }
private:
	pointer _slm_kA_it_lx, _slm_kB_it_lx;
public:
	geom_slm();
	geom_slm( geom_slm<T> const & v);
	geom_slm<T> & operator=( geom_slm<T> const & v);

	const size_type comp_size()const;
	void write();
	void init_geom(T const * CA, T const * CB, matrix_slm<T> const & __mx_slm );
	void init(T const * CA, T const * CB, matrix_slm<T> const & __mx_slm );
	void init(T const * CA, matrix_slm<T> const & __mx_slm );

	inline void set_lmax( _lmax_struct const & __lmax ){ this->_lmax = __lmax;}
	inline void set_l_max( const int __l_max ){ this->_lmax._l_max = __l_max;}
	inline void set_lso_max( const int __lso_max ){ this->_lmax._lso_max = __lso_max;}
	inline void set_la_max( const int __la_max ){ this->_lmax._la_max = __la_max;}
	inline void set_lb_max( const int __lb_max ){ this->_lmax._lb_max = __lb_max;}
	inline const int & l_max()const{ return this->_lmax._l_max;}
	inline const int & lso_max()const{ return this->_lmax._lso_max;}
	inline const int & la_max()const{ return this->_lmax._la_max;}
	inline const int & lb_max()const{ return this->_lmax._lb_max;}
	inline void init_lsize()
	{
		this->_lsize._la_size = this->la_max() + 1;
		this->_lsize._lb_size = this->lb_max() + 1;
		this->_lsize._kA_lsize = this->la_max() + this->l_max() + 1;
		if( this->is_mapping_mid() && this->l_max() < this->lb_max() )
			this->_lsize._kA_lsize = this->la_max() + this->lb_max() + 1;
		this->_lsize._kB_lsize = this->lb_max() + this->l_max() + 1;
	}
	inline const int & la_size()const{ return this->_lsize._la_size;}
	inline const int & lb_size()const{ return this->_lsize._lb_size;}
	inline const int & kA_lsize()const{ return this->_lsize._kA_lsize;}
	inline const int & kB_lsize()const{ return this->_lsize._kB_lsize;}

	// USAGE Interface
	// CX powers usage (e.i. CA):
	// in order to power 'n' of CA.x, CA.y, CA.z:
	//   %CA_x( n ); -- return value - is power 'n' of CA.x
	//   %CA_y( n );
	//   %CA_z( n );
	__GEOM_SLM_CX_INTERFACE( CA );
	__GEOM_SLM_CX_INTERFACE( CB );

	// slm_kX usage (e.i. slm_kA):
	// in order to get spherical values of ang.mom. 'lx', and projection 'm', there are 2 ways:
	//    1.1 %slm_kA_set_lx( lx ); -- set ang.mom. 'lx'
	//    1.2 %slm_kA( m ); -- get spherical with 'lx', 'm'
	//
	//    2.1 %slm_kA(lx, m); -- get it
	__GEOM_SLM_SLMKX_INTERFACE( slm_kA );
	__GEOM_SLM_SLMKX_INTERFACE( slm_kB );
	// Total usage of this class:
	// first (or zero'th), one should init mapping type
	// 0. mapping type initialization
	//   1. %set_mapping( mapping_struct cosnt & )
	//   2. %set_mapping( mapping_enum cosnt & )
	//   3. %set_mapping_min/mid/max() -- choose applicable one (min, mid, or max)
	//
	// second (or then) init input values
	// 1. input values initialization
	//   1.1 %set_lmax( _lmax_struct const & )
	//   1.2 %init_lsize()
	//
	//   2.1 %set_l_max( const int )
	//   2.2 %set_lso_max( const int )
	//   2.3 %set_la_max( const int )
	//   2.4 %set_la_max( const int )
	//   2.5 %init_lsize()
	//
	// then allocate memory and init pointers
	// 2. %write() - memory allocation and pointers initialization
	//
	// calculate powers and sphericals
	// 3. calculate powers of CA_x, CA_y, CA_z; CB_x, ...; and values of spherical harmonics for kA/|kA| = -CA/|CA|, and kB/|kB| = -CB/|CB|
	//   1. %init( CA, CB, mx_slm ); <-- for maximum mapping
	//   2. %init( CA, mx_slm ); <-- for middle mapping
	//   3. nothing for minimum mapping
	//   0. %init_geom( CA, CB, mx_slm ); <-- for all ocasions
	//
	// 4. use it - 'USAGE Interface'
};

#endif//__GEOM_SLM_H__
