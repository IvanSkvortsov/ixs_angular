#ifndef __IXS_ANGULAR_DATA_H__
#define __IXS_ANGULAR_DATA_H__
#include"ixs.angular.map.h"

#include"memory.map.h"
#include"memorystream.h"
#include"matrix.cursor.h"

#include"ixs.omega.h"
#include"geom.slm.h"
#include"alpha.slm.h"

//#define __IXS_ANGULAR_DATA_PRINT
#define __IXS_ANGULAR_DATA_SIZECHECK
//#define __IXS_ANGULAR_DATA_RADIALCHECK

#define __IXS_ANGULAR_DATA_DEBUG
#ifdef  __IXS_ANGULAR_DATA_DEBUG
  #include<cassert>
#endif


template<typename T, typename U>
struct ixs_angular_dat: public memory_map, public memorystream, public ixs_angular_map
{
public:
	__DATA_TYPEDEF( T );
	typedef T float_type;
	typedef int size_type;

	typedef typename memorystream::pos_type pos_type;
	typedef typename memorystream::off_type off_type;
	typedef typename memorystream::seek_dir seek_dir;
	typedef typename memory_map::mode_type mode_type;
	
	typedef struct
	{
		int la, ax, ay, az, lb, bx, by, bz, l;
	} _lxyz_struct;
	typedef struct
	{
		int na, nb, na_p_nb, lmbA, lmbB;
	} _nx_struct;
protected:
	matrix_cursor_1<T> * _M_mx1ang;
	pointer _mx1ang_it;

	inline void error(const char * _method, const char _message[] = "nothing to do here")const
	{
		std::cerr << "Error: [" << this << "] ixs_angular_dat<T,U>::" << _method << ", " << _message << std::endl;
	}
public:
	ixs_angular_dat();
	ixs_angular_dat(ixs_angular_dat<T,U> const & v);
	ixs_angular_dat<T,U> & operator=(ixs_angular_dat<T,U> const & v);

	void memory_create(const char * file, int __flags = O_RDWR |O_CREAT |O_TRUNC, mode_type __mode = S_IRUSR |S_IWUSR |S_IRGRP |S_IROTH);
	void memory_open(const char * file, int __flags = O_RDONLY, mode_type __mode = S_IRUSR |S_IRGRP |S_IROTH);
	
	void write_dat();
	void read_dat();

	inline void sync_stream(){this->memorystream::setbuf( this->memory_map::data(), this->memory_map::size() );}

	inline void set_angular_map( ixs_angular_map const & __angular_map){ this->ixs_angular_map::operator=(__angular_map);}

	// A != C != B  <==> max
	void comp_ang_max( geom_slm<U> & geom_s, alpha_slm<T,U> & alp_s, ixs_omega<T,U> const & ixs_omg );
	// max SemiLocal
	void comp_ang_max_SemiLocal( pointer __p_mx1ang, _lxyz_struct const & _lxyz, geom_slm<U> & geom_s, ixs_omega<T,U> const & ixs_omg );
	void comp_ang_max_SemiLocal( U & value, _lxyz_struct const & _lxyz, _nx_struct const & _nx,
			geom_slm<U> & geom_s, ixs_omega<T,U> & ixs_omg_a, ixs_omega<T,U> & ixs_omg_b );
	// max Local
	void comp_ang_max_Local( pointer __p_mx1ang, _lxyz_struct const & _lxyz,
			geom_slm<U> & geom_s, alpha_slm<T,U> & alp_s, ixs_omega<T,U> const & ixs_omg );
	void comp_ang_max_Local( U & value, _lxyz_struct const & _lxyz, _nx_struct const & _nx,
			geom_slm<U> & geom_s, alpha_slm<T,U> & alp_s, ixs_omega<T,U> & ixs_omg_x );

	// A != C == B  <==> mid
	void comp_ang_mid( geom_slm<U> & geom_s, ixs_omega<T,U> const & ixs_omg );
	// mid SemiLocal
	void comp_ang_mid_SemiLocal( pointer __p_mx1ang, _lxyz_struct const & _lxyz, geom_slm<U> & geom_s, ixs_omega<T,U> const & ixs_omg );
	void comp_ang_mid_SemiLocal( U & value, _lxyz_struct const & _lxyz, _nx_struct const & _nx,
			geom_slm<U> & geom_s, ixs_omega<T,U> & ixs_omg_a, ixs_omega<T,U> & ixs_omg_b );
	// mid Local
	void comp_ang_mid_Local( pointer __p_mx1ang, _lxyz_struct const & _lxyz, geom_slm<U> & geom_s, ixs_omega<T,U> const & ixs_omg );
	void comp_ang_mid_Local( U & value, _lxyz_struct const & _lxyz, _nx_struct const & _nx, geom_slm<U> & geom_s, ixs_omega<T,U> & ixs_omg_x );

	// A == C == B  <==> min
	void comp_ang_min( ixs_omega<T,U> const & ixs_omg );
	// min SemiLocal
	void comp_ang_min_SemiLocal( pointer __p_mx1ang, _lxyz_struct const & _lxyz, ixs_omega<T,U> const & ixs_omg );
	void comp_ang_min_SemiLocal( U & value, _lxyz_struct const & _lxyz, _nx_struct const & _nx, 
			ixs_omega<T,U> & ixs_omg_a, ixs_omega<T,U> & ixs_omg_b );
	// mid Local (not needed)
	// ...

	inline pointer M_mx1ang_data(){return this->_M_mx1ang->data();}
	inline const_pointer M_mx1ang_data()const{return this->_M_mx1ang->data();}
	inline pointer M_mx1ang_data(int i){return this->_M_mx1ang->data(i);}
	inline const_pointer M_mx1ang_data(int i)const{return this->_M_mx1ang->data(i);}
	inline size_type const & M_mx1ang_size()const{return this->_M_mx1ang->size();}

	inline void mx1ang_set_it( int __pos ){ this->_mx1ang_it = this->M_mx1ang_data() + this->ixs_angular_map::map3node_pos(); }
	inline reference mx1ang( int i ){ return this->_mx1ang_it[i]; }
	inline const_reference mx1ang( int i )const{ return this->_mx1ang_it[i]; }
};

#endif//__IXS_ANGULAR_DATA_H__
