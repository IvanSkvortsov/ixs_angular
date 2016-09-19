#ifndef __IXS_ANGULAR_HPP__
#define __IXS_ANGULAR_HPP__
#include"ixs.angular.h"
#include<cstdlib>// exit
#include<iomanip>// setw
#ifdef  __IXS_ANGULAR_DEBUG
  #include<cassert>
#endif

template<typename T, typename U> ixs_angular<T,U>::ixs_angular() : memory_map(), memorystream(), _lmax(),
	_M_lmax(), _M_map(), _M_map_lmb(), _M_map_nx2(), _M_node0(), _M_mxang(){}

template<typename T, typename U> ixs_angular<T,U>::ixs_angular(ixs_angular<T,U> const & v) : memory_map(), memorystream(v), _lmax(v._lmax),
	_M_lmax(v._M_lmax), _M_map(v._M_map), _M_map_lmb(v._M_map_lmb), _M_map_nx2(v._M_map_nx2), _M_node0(v._M_node0), _M_mxang(v._M_mxang){}

template<typename T, typename U>
ixs_angular<T,U> & ixs_angular<T,U>::operator=(ixs_angular<T,U> const & v)
{
	if( this == &v )
		return *this;
	if( this->data() || this->size() )
		this->close();
	ixs_angular<T,U> tmp( v );
	memcpy( this, &tmp, sizeof(*this) );
	return *this;
}

// memory_create
template<typename T, typename U>
void ixs_angular<T,U>::memory_create( const char * file, int __flags, mode_type __mode )
{
	static const int _AUG_BYTES = 1024;
	size_type __size = this->comp_size();
	int _st = this->memory_map::create( file, __size + _AUG_BYTES, __flags, __mode );
	if( _st )
	{
		this->error("memory_create", "from memory_map");
		std::cerr << "file      '" << file << "'" << std::endl;
		std::cerr << "  __size : " << std::setw(12) << __size << std::endl;
		std::cerr << "map.size : " << std::setw(12) << this->size() << std::endl;
		std::cerr << "map.data : [" << std::setw(12) << this->data() << "]" << std::endl;
		exit(1);
	}
	this->sync_stream();
	this->write_ang();
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type ixs_angular<T,U>::write_ang()
{
	memorystream & ms = *this;
	size_type _seek_start = ms.tell();

	this->_M_lmax = (_lmax_struct *)ms.getcur();
	*this->_M_lmax = this->_lmax;
	ms.seek( sizeof(_lmax_struct), seek_dir::cur );

	this->_M_mapping_t = (mapping_enum *)ms.getcur();
	*this->_M_mapping_t = this->_mapping_t;
	ms.seek( sizeof(mapping_enum), seek_dir::cur );

	this->_M_mappos = (_mappos_struct *)ms.getcur();
	ms.seek( sizeof(_mappos_struct), seek_dir::cur );

	union{
		void * _void;
		matrix_cursor_2<_min1_struct> * _map_lmb;
		matrix_cursor_3<_pos1_struct> * _node0;
		matrix_cursor_3<int> * _map_nx2;
		matrix_cursor_1<T> * _mxang;
	} __cnvrt;

	__cnvrt._void = ms.getcur();
	this->_M_map_lmb = __cnvrt._map_lmb;
	this->_M_map_lmb->n() = this->map_lmb_N();
	this->_M_map_lmb->m() = this->map_lmb_M();
	this->_M_map_lmb->init_size();
	ms.seek( sizeof(size_struct<2>) + sizeof(_min1_struct)*this->_M_map_lmb->size(), seek_dir::cur );

	__cnvrt._void = ms.getcur();
	this->_M_map_nx2 = __cnvrt._map_nx2;
	this->_M_map_nx2->n() = this->map_nx2_N();
	this->_M_map_nx2->m() = this->map_nx2_M();
	this->_M_map_nx2->p() = this->map_nx2_P();
	this->_M_map_nx2->init_size();
	ms.seek( sizeof(size_struct<3>) + sizeof(int)*this->_M_map_nx2->size(), seek_dir::cur );

	__cnvrt._void = ms.getcur();
	this->_M_node0 = __cnvrt._node0;
	this->_M_node0->n() = this->node0_N();
	this->_M_node0->m() = this->node0_M();
	this->_M_node0->p() = this->node0_P();
	this->_M_node0->init_size();
	ms.seek( sizeof(size_struct<3>) + sizeof(int)*this->_M_node0->size(), seek_dir::cur );

	/*
	__cnvrt._void = ms.getcur();
	this->_M_mxang = __cnvrt._mxang;
	this->_M_mxang->size() = this->_mxang_size();
	ms.seek( sizeof(size_struct<1>) + sizeof(T) * this->_M_mxang->size(), seek_dir::cur );
	*/

	return ms.tell() - _seek_start;
}

template<typename T, typename U> const typename ixs_angular<T,U>::size_type comp_size()const
{
	size_type __size = 0;
	__size += sizeof(_lmax_struct);
	__size += sizeof(mapping_enum);
	__size += sizeof(_mappos_struct);

	__size += this->comp_lmb_size();
	__size += this->comp_nx2_size();
	__size += this->comp_node0_size();
	//__size += this->comp_mxang_size();
	return __size;
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type comp_lmb_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<2>);
	__size += sizeof(_min1_struct) * this->map_lmb_N() * this->map_lmb_M();
	return __size;
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type comp_nx2_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<3>);
	__size += sizeof(int) * this->map_nx2_N() * this->map_nx2_M() * this->map_nx2_P();
	return __size;
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type comp_node0_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<3>);
	__size += sizeof(int) * this->node0_N() * this->node0_M() * this->node0_P();
	return __size;
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type comp_mxang_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<1>);
	this->_mxang_size = this->mxang_SIZE();
	__size += sizeof(T) * this->_mxang_size;
	return __size;
}

/*
template<typename T, typename U> const typename ixs_angular<T,U>::size_type mxang_SIZE()const
{
	size_type __size = 0;
}
*/
// init
template<typename T, typename U>
const typename ixs_angular<T,U>::size_type ixs_angular<T,U>::init()
{
	this->mx_lmb_init();
	this->mx_nxn_init();
#ifdef  __IXS_ANGULAR_ALPMAP
	this->mx_alpmap_init();
#endif
	return this->mx_big_init();
}
template<typename T, typename U> const typename ixs_angular<T,U>::size_type
ixs_angular<T,U>::init(const int * __alpmapA, const int __a_size, const int * __alpmapB, const int __b_size)
{
	this->mx_lmb_init();
	this->mx_nxn_init();
#ifdef  __IXS_ANGULAR_ALPMAP
	this->mx_alpmap_init( __alpmapA, __a_size, __alpmapB, __b_size );
#endif
	return this->mx_big_init();
}
// mx_lmb
template<typename T, typename U> void ixs_angular<T,U>::mx_lmb_init()
{
	// 2d mx_lmb
	matrix_cursor_2<mxlmb_type> & __mx_lmb = *this->_M_mx_lmb._matrix_lmb;
	mxlmb_type * __p_lmb = __mx_lmb.data(), * __p_lmb_i;
	for(int l = 0; l < this->_src_max._l_max; ++l, __p_lmb += __mx_lmb.m() )// n() = (l_max + lso_max + 1)
	{
		__p_lmb_i = __p_lmb;
		for(int n = 0; n < __mx_lmb.m(); ++n, ++__p_lmb_i )// m() = max( la_max + 1, lb_max + 1 )
		{
			__p_lmb_i->_min = (l - n < 0 ? (l + n)%2 : l - n);
			__p_lmb_i->_size = (l + n - __p_lmb_i->_min)/2 + 1;
		}
	}
}
// mx_nxn
template<typename T, typename U> void ixs_angular<T,U>::mx_nxn_init()
{
	// 2d mx_lmb
	matrix_cursor_2<mxlmb_type> & __mx_lmb = *this->_M_mx_lmb._matrix_lmb;
	mxlmb_type * __p_lmb = __mx_lmb.data(), * __p_lmb_i, * __p_lmb_j;
	// 3d mx_nxn
	matrix_cursor_3<mxnxn_type> & __mx_nxn = *this->_M_mx_nxn._matrix_nxn;
	mxnxn_type * __p_nxn = __mx_nxn.data(), * __p_nxn_i;
	__p_lmb = __mx_lmb.data();
	for(int l = 0; l < this->_src_max._l_max; ++l, __p_lmb += __mx_lmb.m(), __p_nxn += __mx_nxn.m() * __mx_nxn.p() )
	{
		__p_nxn_i = __p_nxn;
		__p_lmb_i = __p_lmb;
		for(int na = 0; na < __mx_nxn.m(); ++na, ++__p_lmb_i)
		{
			__p_lmb_j = __p_lmb;
			for(int nb = 0; nb < __mx_nxn.p(); ++nb, ++__p_lmb_j, ++__p_nxn_i)
			{
				*__p_nxn_i = __p_lmb_i->_size;
				*__p_nxn_i *= __p_lmb_j->_size;
			}
		}
	}
	// l = this->_src_max._l_max
	__p_nxn_i = __p_nxn;
	for(int na = 0; na < __mx_nxn.m(); ++na)
		for(int nb = 0; nb < __mx_nxn.p(); ++nb, ++__p_nxn_i)
			*__p_nxn_i = (na+nb)/2 + 1;
}
#ifdef  __IXS_ANGULAR_ALPMAP
// mx_alpmap
template<typename T, typename U> void ixs_angular<T,U>::mx_alpmap_init()
{
	for(int i = 0; i < this->mx_alpmapA().size(); ++i)
		this->alpmapA(i) = 1;
	for(int i = 0; i < this->mx_alpmapB().size(); ++i)
		this->alpmapB(i) = 1;
}
template<typename T, typename U>
void ixs_angular<T,U>::mx_alpmap_init(const int * __alpmapA, const int __a_size, const int * __alpmapB, const int __b_size)
{
	if( this->mx_alpmapA().size() != __a_size || this->mx_alpmapB().size() != __b_size )
	{
		this->error("mx_alpmap_init", "size conflict");
		std::cerr << "a_size : " << __a_size << std::endl;
		std::cerr << "b_size : " << __b_size << std::endl;
		std::cerr << "mx_alpmapA.size : " << this->mx_alpmapA().size() << std::endl;
		std::cerr << "mx_alpmapB.size : " << this->mx_alpmapB().size() << std::endl;
		exit(1);
	}
	int * p = &this->alpmapA(0);
	for(int i = 0; i < __a_size; ++i, ++p)
		*p = __alpmapA[i];
	p = &this->alpmapB(0);
	for(int i = 0; i < __b_size; ++i, ++p)
		*p = __alpmapB[i];
}
#endif
// mx_big
template<typename T, typename U> const typename ixs_angular<T,U>::size_type ixs_angular<T,U>::mx_big_init()
{
	// 3d mx_nxn
	matrix_cursor_3<mxnxn_type> & __mx_nxn = *this->_M_mx_nxn._matrix_nxn;
	mxnxn_type * __p_nxn = __mx_nxn.data(), * __p_nxn_i;
	// 3d mx_big
	matrix_cursor_3<mxbig_type> & __mx_big = *this->_M_mx_big._matrix_big;
	mxbig_type * __p_big = __mx_big.data(), * __p_big_i;
	size_type __size, __pos = 0;
	size_type __mx_big_mp = __mx_big.m() * __mx_big.p();
	size_type __mx_nxn_mp = __mx_nxn.m() * __mx_nxn.p();
#if defined (__IXS_ANGULAR_PRINT)
	std::cout << std::setw(8) << "iter" << std::setw(8) << "i" <<
		std::setw(4) << "la" << std::setw(4) << "ila" << ' ' <<
		std::setw(4) << "lb" << std::setw(4) << "ilb" << ' ' <<
		std::setw(4) << "l" <<
		std::setw(14) << "pos" << std::setw(8) << "size_i" << std::endl << std::endl;
	int iter = 0;
#endif
	for(int la = 0, ila_size = 1; la <= this->_src_cur._la_max; ++la, ila_size += (la + 1))
	{
		//for(int i_la = 0; i_la < ila_size; ++i_la, __p_big += __mx_big.m() * __mx_big.p() )
		for(int i_la = 0; i_la < ila_size; ++i_la, __p_big += __mx_big_mp )
		{
			__p_big_i = __p_big;
			for(int lb = 0, ilb_size = 1; lb <= this->_src_cur._lb_max; ++lb, ilb_size += (lb + 1))
			{
				for(int i_lb = 0; i_lb < ilb_size; ++i_lb, __p_big_i += __mx_big.p() )
				{
					__p_nxn = __mx_nxn.data();
					//for(int l = 0; l <= this->_src_cur._l_max; ++l, __p_nxn += __mx_nxn.m() * __mx_nxn.p())
					for(int l = 0; l < this->_src_cur._l_max; ++l, __p_nxn += __mx_nxn_mp )// not typo '%l <= %_l_max'
					//for(int l = 0; l <= this->_src_cur._l_max; ++l, __p_nxn += __mx_nxn_mp )// not typo '%l <= %_l_max'
					{
						__p_nxn_i = __p_nxn;
						__size = 0;
						for(int na = 0; na <= la; ++na, __p_nxn_i += __mx_nxn.p())
						{
							for(int nb = 0; nb <= lb; ++nb)
							{
								__size += __p_nxn_i[nb];
							}
						}
						__p_big_i[l]._pos = __pos;
#ifdef  __IXS_ANGULAR_ALPMAP
						if( l == this->_src_cur._l_max )
						{
							__size *= this->alpmapA(la);
							__size *= this->alpmapB(lb);
						}
#endif
						__p_big_i[l]._size = __size;
						__pos += __size;
#if defined (__IXS_ANGULAR_PRINT)
						std::cout << std::setw(8) << iter++ << std::setw(8) << &__p_big_i[l] - __mx_big.data() <<
							std::setw(4) << la << std::setw(4) << i_la << ' ' <<
							std::setw(4) << lb << std::setw(4) << i_lb << ' ' <<
							std::setw(4) << l <<
							std::setw(14) << __pos << std::setw(8) << __size << std::endl;
#endif
					}
#if defined (__IXS_ANGULAR_PRINT)
					std::cout << "---------------" << std::endl;
#endif
				}
			}
#if defined (__IXS_ANGULAR_PRINT)
			std::cout << std::endl;
#endif
		}
	}
	return __pos;
}
// memory_open
template<typename T, typename U> void ixs_angular<T,U>::memory_open(char const * file, int __flags, mode_type __mode)
{
	if( this->data() || this->size() )
		this->close();
	int _st = this->::memory_map::open( file, __flags, __mode );
	if( _st )
	{
		this->error("memory_open", "can't open");
		std::cerr << "file '" << file << '\'' << std::endl;
		exit(1);
	}
	this->sync_stream();
	this->seek( 0 );
	this->read_matrix();
}

#endif//__IXS_ANGULAR_HPP__
