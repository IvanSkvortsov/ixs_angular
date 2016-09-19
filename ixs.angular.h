#ifndef __IXS_ANGULAR_H__
#define __IXS_ANGULAR_H__
#include"memory.map.h"
#include"memorystream.h"
#include"matrix.cursor.h"
#include"mapping.t.h"

template<typename T, typename U>
struct ixs_angular : public memory_map, public memorystream, public mapping_struct
{
public:
	typedef typename memory_map::mode_type mode_type;
	typedef typename memorystream::pos_type pos_type;
	typedef typename memorystream::off_type off_type;
	typedef typename memorystream::seek_dir seek_dir;
	typedef typename alpha_map::mapping_type mapping_type;

	__DATA_TYPEDEF( T );
	typedef typename size_struct<1>::size_type size_type;
#pragma pack( push, 4 )
	typedef struct
	{
		int _l_max, _lso_max, _la_max, _lb_max;
	} _lmax_struct;
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
		int _lmb, _nx2, _node0;
	} _mappos_struct;
	typedef struct
	{
		int _lmb_size, _nx2_size, _node0_size, _mxang_size;
	} _mxsize_struct;
#pragma pack( pop )
protected:
	//_mxsize_struct _mxsize;
	size_type _mxang_size;
	_lmax_struct _lmax;
	_lmax_struct * _M_lmax;
	_mappos_struct * _M_mappos;
	mapping_struct * _M_mapping_t;

	matrix_cursor_2<_min1_struct> * _M_map_lmb;
	matrix_cursor_3<int> * _M_map_nx2;
	matrix_cursor_3<_pos1_struct> * _M_node0;
	matrix_cursor_1<T> * _M_mxang;
public:

	inline void set_lmax( _lmax_struct const & __lmax){ this->_lmax = __lmax;}
	const size_type comp_size()const;
	const size_type comp_lmb_size()const;
	const size_type comp_nx2_size()const;
	const size_type comp_node0_size()const;
	//const size_type comp_mxang_size()const;

	inline const size_type map_lmb_N()const{return (this->_lmax._l_max + this->_lmax._lso_max + 1);}
	inline const size_type map_lmb_M()const
	{
		return (this->_lmax._la_max > this->_lmax._lb_max ? (this->_lmax._la_max + 1) : (this->_lmax._lb_max + 1));
	}

	inline const size_type map_nx2_N()const{return (this->_lmax._l_max + 1 + this->_lmax._lso_max);}
	inline const size_type map_nx2_M()const{return (this->_lmax._la_max + 1);}
	inline const size_type map_nx2_P()const{return (this->_lmax._lb_max + 1);}

	inline const size_type node0_N()const
	{
		size_type __size = this->_lmax._la_max + 1;
		__size *= (this->_lmax._la_max + 2);
		__size *= (this->_lmax._la_max + 3);
		__size /= 6;
		return __size;
	}
	inline const size_type node0_M()const
	{
		size_type __size = this->_lmax._lb_max + 1;
		__size *= (this->_lmax._lb_max + 2);
		__size *= (this->_lmax._lb_max + 3);
		__size /= 6;
		return __size;
	}
	inline const size_type node0_P()const{return (this->_lmax._l_max + 1 + this->_lmax._lso_max);}
	// P.S. it should be like as if one had %lmax then size is (%lmax+1), l = 0, 1, ..., %lmax;
	// but here _l = 0, 1, ..., %_l_max; _lso = 1, 2, ..., %_lso_max;
	// that's why size = (%l_max + %lso_max)
};

#endif//__IXS_ANGULAR_H__
