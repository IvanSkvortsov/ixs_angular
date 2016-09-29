#ifndef __GEOM_SLM_HPP__
#define __GEOM_SLM_HPP__
#include"geom.slm.h"
#include<cstring>// memcpy, memset
#include"vector.3d.h"// vector_3d

template<typename T> geom_slm<T>::geom_slm(): std::vector<T>(), mapping_struct(maximum), _lmax(), _CA(), _CB(), _slm_kA(), _slm_kB()
	,_lsize(), _slm_kA_it_lx(), _slm_kB_it_lx(){}
template<typename T> geom_slm<T>::geom_slm(geom_slm<T> const & v):
	std::vector<T>(), mapping_struct(v), _lmax(v._lmax), _CA(v._CA), _CB(v._CB), _slm_kA(v._slm_kA), _slm_kB(v._slm_kB)
	,_lsize(v._lsize), _slm_kA_it_lx(v._slm_kA_it_lx), _slm_kB_it_lx(v._slm_kB_it_lx){}

template<typename T> geom_slm<T> & geom_slm<T>::operator=(geom_slm<T> const & v)
{
	if( this == &v )
		return *this;
	geom_slm<T> t1( v ), t2;
	memcpy( &t2, this, sizeof(*this) );
	memcpy( this, &t1, sizeof(*this) );
	memset( &t1, 0, sizeof(*this) );
	return *this;
}

template<typename T> const typename geom_slm<T>::size_type geom_slm<T>::comp_size()const
{
	size_type __size = 0;
	__size += this->la_size();
	__size += this->lb_size();
	__size *= 3;
	__size += (this->kA_lsize() * this->kA_lsize());
	__size += (this->kB_lsize() * this->kB_lsize());
	return __size;
}

template<typename T> void geom_slm<T>::write()
{
	this->resize( this->comp_size() );
	this->_CA.x = this->data();
	this->_CA.y = this->_CA.x + this->la_size();
	this->_CA.z = this->_CA.y + this->la_size();

	this->_CB.x = this->_CA.z + this->la_size();
	this->_CB.y = this->_CB.x + this->lb_size();
	this->_CB.z = this->_CB.y + this->lb_size();

	this->_slm_kA = this->_CB.z + this->lb_size();
	this->_slm_kB = this->_slm_kA + (this->kA_lsize() * this->kA_lsize());
}

template<typename T>
inline static void comp_powers_of( T const & value, const int max_power, T * storage, T & u )
{
	u = 1; *storage++ = 1;
	for(int i = 1; i <= max_power; ++i, ++storage)
	{
		u *= value;
		*storage = u;
	}
}

template<typename T>
inline static void slm_pown( T * __data, matrix_slm<T> const & mx_slm, T const * _vx, T const * _vy, T const * _vz, const int l_max)
{
	if( mx_slm.lmb_max() < l_max )
	{
		std::cerr << "Error: [slm_pown]" << std::endl;
		std::cerr << "mx_slm.lmax < lmax" << std::endl;
		std::cerr << "mx_slm.lmax : " << mx_slm.lmb_max() << std::endl;
		std::cerr << "lmax : " << l_max << std::endl;
		exit(1);
	}
	T u;
	typename matrix_slm<T>::const_pointer _slm = mx_slm.data();
	typename matrix_slm<T>::value_type::const_pointer _pol;
	for(int l = 0; l <= l_max; ++l)
	{
		for(int m = -l; m <= l; ++m, ++__data, ++_slm)
		{
			*__data = 0;
			_pol = _slm->data();
			for(int i = 0; i < _slm->size(); ++i, ++_pol )
			{
				u  = _vx[ _pol->x[0] ];
				u *= _vy[ _pol->x[1] ];
				u *= _vz[ _pol->x[2] ];
				u *= _pol->d;
				*__data += u;
			}
		}
	}
}

template<typename T>
void geom_slm<T>::init_slmkX( T * __slm_kX, T * CX_x, T * CX_y, T * CX_z, const int lx_max, const int kx_lmax,
		T const * CX, matrix_slm<T> const & mx_slm )
{
	vector_3d<T> kx( CX );

	T u;
	comp_powers_of<T>( kx.x(), lx_max, CX_x, u );
	comp_powers_of<T>( kx.y(), lx_max, CX_y, u );
	comp_powers_of<T>( kx.z(), lx_max, CX_z, u );

	T _kx_len;
	kx.neg();
	kx.len( _kx_len );
	kx /= _kx_len;

	T powk_x[kx_lmax+1], powk_y[kx_lmax+1], powk_z[kx_lmax+1];
	comp_powers_of<T>( kx.x(), kx_lmax, powk_x, u );
	comp_powers_of<T>( kx.y(), kx_lmax, powk_y, u );
	comp_powers_of<T>( kx.z(), kx_lmax, powk_z, u );

	slm_pown<T>( __slm_kX, mx_slm, powk_x, powk_y, powk_z, kx_lmax );
}

template<typename T>
void geom_slm<T>::init( T const * CA, matrix_slm<T> const & mx_slm )
{
	this->init_slmkX( this->_slm_kA, this->_CA.x, this->_CA.y, this->_CA.z, this->la_max(), this->kA_lsize()-1, CA, mx_slm );
}

template<typename T>
void geom_slm<T>::init( T const * CA, T const * CB, matrix_slm<T> const & mx_slm )
{
	this->init_slmkX( this->_slm_kA, this->_CA.x, this->_CA.y, this->_CA.z, this->la_max(), this->kA_lsize()-1, CA, mx_slm );
	this->init_slmkX( this->_slm_kB, this->_CB.x, this->_CB.y, this->_CB.z, this->lb_max(), this->kB_lsize()-1, CB, mx_slm );
}

template<typename T>
void geom_slm<T>::init_geom( T const * CA, T const * CB, matrix_slm<T> const & mx_slm )
{
	switch( this->get_mapping() )
	{
	case minimum:
		break;
	case middle :
		this->init( CA, mx_slm );
		break;
	case maximum:
		this->init( CA, CB, mx_slm );
		break;
	default:
		this->error("init_geom", "unknowen mapping type");
		std::cerr << "current mapping type : " << this->get_mapping() << std::endl;
		exit(1);
	}
}
#define __GEOM_SLM_SPEC( type )\
template class geom_slm<type>;\
void comp_powers_of( type const & value, const int max_power, type * storage, type & u );\
void slm_pown( type * __data, matrix_slm<type> const & mx_slm, type const * _vx, type const * _vy, type const * _vz, const int l_max);

#endif//__GEOM_SLM_HPP__
