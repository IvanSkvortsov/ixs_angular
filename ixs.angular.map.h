#ifndef __IXS_ANGULAR_MAP_H__
#define __IXS_ANGULAR_MAP_H__
#include"memorystream.h"
#include"matrix.cursor.h"
#include"mapping.t.h"
#include"alpha.map.h"

#define __IXS_ANGULAR_MAP_DEBUG
#ifdef  __IXS_ANGULAR_MAP_DEBUG
  #include<cassert>
#endif

#define __IXS_ANGULAR_MAP_MAP2LMB( lmbX )\
inline void map2##lmbX##_set_l( const size_type l ){ this->_map2##lmbX##_it_l = this->_M_map_lmb->data() + this->M_map_lmb_m() * l;}\
inline void map2##lmbX##_set_nx( const size_type nx ){ this->_map2##lmbX##_it_nx = this->_map2##lmbX##_it_l + nx;}\
inline int & map2##lmbX##_min (){ return this->_map2##lmbX##_it_nx->_min;}\
inline int & map2##lmbX##_size(){ return this->_map2##lmbX##_it_nx->_size;}\
inline int const & map2##lmbX##_min ()const{ return this->_map2##lmbX##_it_nx->_min;}\
inline int const & map2##lmbX##_size()const{ return this->_map2##lmbX##_it_nx->_size;}

#define __IXS_ANGULAR_MAP_MX3( mx3, mx3_value_type )\
inline const size_type & M_##mx3##_n()const{return this->_M_##mx3->n();}\
inline const size_type & M_##mx3##_m()const{return this->_M_##mx3->m();}\
inline const size_type & M_##mx3##_p()const{return this->_M_##mx3->p();}\
inline const size_type & M_##mx3##_size()const{return this->_M_##mx3->size();}\
inline const mx3_value_type * M_##mx3##_data()const{return this->_M_##mx3->data();}

#define __LMAX_STRUCT_INTERFACE( lx_max )\
inline const int & lx_max()const{ return this->_M_lmax->_##lx_max;}\
inline int & lx_max(){ return this->_M_lmax->_##lx_max;}

struct ixs_angular_map : public mapping_struct
{
public:
	typedef typename memorystream::pos_type pos_type;
	typedef typename memorystream::off_type off_type;
	typedef typename memorystream::seek_dir seek_dir;

	typedef typename size_struct<1>::size_type size_type;
	typedef typename alpha_siz::_lmax_struct _lmax_struct;
#pragma pack( push, 4 )
	typedef struct
	{
		int _pos, _size;
	} _pos1_struct;
	typedef struct
	{
		int _min, _size;
	} _min1_struct;
	typedef struct
	{
		int _lmb, _nx2, _node;
	} _mappos_struct;
#pragma pack( pop )
protected:
	size_type _mx1ang_size;
	_lmax_struct _lmax;
	_lmax_struct * _M_lmax;
	_mappos_struct * _M_mappos;
	mapping_struct * _M_mapping_t;

	matrix_cursor_2<_min1_struct> * _M_map_lmb;
	matrix_cursor_3<int> * _M_map_nx2;
	matrix_cursor_3<_pos1_struct> * _M_node;

	void init_node_min(size_type & __pos, const int & la, const int & lb);
	void init_node_mid(size_type & __pos, const int & la, const int & lb);
	void init_node_max(size_type & __pos, const int & la, const int & lb, alpha_map & alp_m);

	inline void error(char const * _method, const char _message[] = "nothing to do here" )const
	{ std::cerr << "Error: [" << this << "] ixs_angular_map::" << _method << ", " << _message << std::endl;}
private:
	_min1_struct * _map2lmbA_it_l, * _map2lmbA_it_nx;
	_min1_struct * _map2lmbB_it_l, * _map2lmbB_it_nx;
	int _map3nx2_it_l, _map3nx2_it_na, * _map3nx2_it_nb;
	int _map3node_it_la, _map3node_it_ia, _map3node_it_lb, _map3node_it_ib;
	_pos1_struct * _map3node_it_l;
public:
	// Usage:
	// creation
	//   1. you should init %_lmax; e.i. %set_lmax; it is used in %comp_size methods
	//   2. %write_map
	//   3. %init_map() or %init_map( alpha_map & )
	//   4. use it
	// opening
	//   1. %read_map
	//   2. use it
	ixs_angular_map();
	ixs_angular_map( ixs_angular_map const & v);
	ixs_angular_map & operator=(ixs_angular_map const & v);
	// read, write
	const size_type write_map( memorystream & ms );
	const size_type read_map( memorystream & ms );
	// init
	void init_map();// min, mid
	void init_map( alpha_map & alp_m );// max
	void init_map_lmb();
	void init_map_nx2();
	const size_type init_node_min();
	const size_type init_node_mid();
	const size_type init_node_max(alpha_map & alp_m);
	// comp_size
	const size_type comp_size()const;
	const size_type comp_lmb_size()const;
	const size_type comp_nx2_size()const;
	const size_type comp_node_size()const;
	// mx1ang_size
	inline const size_type & mx1ang_size()const{ return this->_mx1ang_size;}
	inline size_type & mx1ang_size(){ return this->_mx1ang_size;}
	// lmax interface
	inline void set_lmax( _lmax_struct const & __lmax){ this->_lmax = __lmax;}
	__LMAX_STRUCT_INTERFACE( l_max );
	__LMAX_STRUCT_INTERFACE( lso_max );
	__LMAX_STRUCT_INTERFACE( la_max );
	__LMAX_STRUCT_INTERFACE( lb_max );

	inline const size_type map_lmb_N()const{return (this->_lmax._l_max + this->_lmax._lso_max + 1);}
	inline const size_type map_lmb_M()const
	{
		return (this->_lmax._la_max > this->_lmax._lb_max ? (this->_lmax._la_max + 1) : (this->_lmax._lb_max + 1));
	}

	inline const size_type map_nx2_N()const{return (this->_lmax._l_max + 1 + this->_lmax._lso_max);}
	inline const size_type map_nx2_M()const{return (this->_lmax._la_max + 1);}
	inline const size_type map_nx2_P()const{return (this->_lmax._lb_max + 1);}

	inline const size_type node_N()const
	{
		size_type __size = this->_lmax._la_max + 1;
		__size *= (this->_lmax._la_max + 2);
		__size *= (this->_lmax._la_max + 3);
		__size /= 6;
		return __size;
	}
	inline const size_type node_M()const
	{
		size_type __size = this->_lmax._lb_max + 1;
		__size *= (this->_lmax._lb_max + 2);
		__size *= (this->_lmax._lb_max + 3);
		__size /= 6;
		return __size;
	}
	inline const size_type node_P()const{return (this->_lmax._l_max + 1 + this->_lmax._lso_max);}
	
	// USAGE Interface
	// lmb
	inline const size_type & M_map_lmb_n()const{ return this->_M_map_lmb->n();}
	inline const size_type & M_map_lmb_m()const{ return this->_M_map_lmb->m();}
	inline const size_type & M_map_lmb_size()const{ return this->_M_map_lmb->size();}
	inline const _min1_struct * M_map_lmb_data()const{ return this->_M_map_lmb->data();}
	__IXS_ANGULAR_MAP_MAP2LMB( lmbA );
	__IXS_ANGULAR_MAP_MAP2LMB( lmbB );
	// e.i.
	//    %map2lmbA_set_l( l );  -- l is angular momentum of pseudo potential
	//    %map2lmbA_set_nx( na ); -- na is sum of a,b,c, which are <= (ax,ay,az); and ax+ay+az = la
	// -> %map2lmbA_min(); -- minimum value of lambda A (l < na ? (l+na)%2 : (l-na))
	// -> %map2lmbA_size(); -- number of lambda A; size = ((l + na - lmb_min)/2 + 1); lmb = lmb_min, lmb_min+2, ..., l+na;

	// nx2
	__IXS_ANGULAR_MAP_MX3( map_nx2, int );
	inline void map3nx2_set_lmax(){ this->map3nx2_set_l( this->l_max() );}
	inline void map3nx2_set_l( const size_type l){ this->_map3nx2_it_l = l * this->M_map_nx2_m();}
	inline void map3nx2_set_na( const size_type na)
	{
		this->_map3nx2_it_na = this->_map3nx2_it_l +  na;
		this->_map3nx2_it_na *= this->M_map_nx2_p();
	}
	inline void map3nx2_set_nb( const size_type nb){ this->_map3nx2_it_nb = this->_M_map_nx2->data() + this->_map3nx2_it_na + nb;}
	inline int & map3nx2(){return *this->_map3nx2_it_nb;}
	inline int const & map3nx2()const{return *this->_map3nx2_it_nb;}
	// node
	__IXS_ANGULAR_MAP_MX3( node, _pos1_struct );
	inline void map3node_set_la( const size_type lx)
	{
		this->_map3node_it_la = lx * (lx + 1) * (lx + 2);
		this->_map3node_it_la /= 6;
	}
	inline void map3node_set_ia( const size_type ix)
	{
		this->_map3node_it_ia = this->_map3node_it_la + ix;
		this->_map3node_it_ia *= this->M_node_m();
	}
	inline void map3node_set_lb( const size_type lx)
	{
		this->_map3node_it_lb = lx * (lx + 1) * (lx + 2);
		this->_map3node_it_lb /= 6;
		this->_map3node_it_lb += this->_map3node_it_ia;
	}
	inline void map3node_set_ib( const size_type ix)
	{
		this->_map3node_it_ib = this->_map3node_it_lb + ix;
		this->_map3node_it_ib *= this->M_node_p();
	}
	inline void map3node_set_lmax(){ this->map3node_set_l( this->l_max() );}
	inline void map3node_set_l( const size_type l )
	{
		this->_map3node_it_l  = this->_M_node->data();
		this->_map3node_it_l += this->_map3node_it_ib;
		this->_map3node_it_l += l;
#ifdef  __IXS_ANGULAR_MAP_DEBUG
		assert( this->_map3node_it_l == this->_M_node->data( this->_map3node_it_ib + l ) );
#endif
	}
	inline int & map3node_pos (){ return this->_map3node_it_l->_pos;}
	inline int & map3node_size(){ return this->_map3node_it_l->_size;}
	inline int const & map3node_pos ()const{ return this->_map3node_it_l->_pos;}
	inline int const & map3node_size()const{ return this->_map3node_it_l->_size;}
};

#endif//__IXS_ANGULAR_MAP_H__
