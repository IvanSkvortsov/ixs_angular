#ifndef __IXS_ANGULAR_MAP_HPP__
#define __IXS_ANGULAR_MAP_HPP__
#include"ixs.angular.map.h"
#include<cstdlib>// exit
#include<iomanip>// setw
#ifdef  __IXS_ANGULAR_MAP_DEBUG
  #include<cassert>
#endif

ixs_angular_map::ixs_angular_map() : memory_map(), memorystream(), _lmax(),
	_M_lmax(), _M_map(), _M_map_lmb(), _M_map_nx2(), _M_node(), _M_mxang(){}

ixs_angular_map::ixs_angular_map(ixs_angular_map const & v) : memory_map(), memorystream(v), _lmax(v._lmax),
	_M_lmax(v._M_lmax), _M_map(v._M_map), _M_map_lmb(v._M_map_lmb), _M_map_nx2(v._M_map_nx2), _M_node(v._M_node), _M_mxang(v._M_mxang){}

template<typename T, typename U>
ixs_angular_map & ixs_angular_map::operator=(ixs_angular_map const & v)
{
	if( this == &v )
		return *this;
	if( this->data() || this->size() )
		this->close();
	ixs_angular_map tmp( v );
	memcpy( this, &tmp, sizeof(*this) );
	return *this;
}

// memory_create
void ixs_angular_map::memory_create( const char * file, int __flags, mode_type __mode )
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
	this->write_map();
}
// memory_open
void ixs_angular_map::memory_open( const char * file, int __flags, mode_type __mode )
{
	if( this->data() || this->size() )
		this->close();
	int _st = this->memory_map::open( file, __flags, __mode );
	if( _st )
	{
		this->error("memory_open", "from memory_map");
		std::cerr << "file      '" << file << "'" << std::endl;
		std::cerr << "map.size : " << std::setw(12) << this->size() << std::endl;
		std::cerr << "map.data : [" << std::setw(12) << this->data() << "]" << std::endl;
		exit(1);
	}
	this->sync_stream();
	this->read_map();
}
const typename ixs_angular_map::size_type ixs_angular_map::write_map()
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
		matrix_cursor_3<_pos1_struct> * _node;
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
	this->_M_node = __cnvrt._node;
	this->_M_node->n() = this->node_N();
	this->_M_node->m() = this->node_M();
	this->_M_node->p() = this->node_P();
	this->_M_node->init_size();
	ms.seek( sizeof(size_struct<3>) + sizeof(int)*this->_M_node->size(), seek_dir::cur );

	return ms.tell() - _seek_start;
}

const typename ixs_angular_map::size_type ixs_angular_map::read_map()
{
	memorystream & ms = *this;
	size_type _seek_start = ms.tell();

	this->_M_lmax = (_lmax_struct *)ms.getcur();
	this->_lmax = *this->_M_lmax;
	ms.seek( sizeof(_lmax_struct), seek_dir::cur );

	this->_M_mapping_t = (mapping_enum *)ms.getcur();
	this->_mapping_t = *this->_M_mapping_t;
	ms.seek( sizeof(mapping_enum), seek_dir::cur );

	this->_M_mappos = (_mappos_struct *)ms.getcur();
	ms.seek( sizeof(_mappos_struct), seek_dir::cur );

	union{
		void * _void;
		matrix_cursor_2<_min1_struct> * _map_lmb;
		matrix_cursor_3<_pos1_struct> * _node;
		matrix_cursor_3<int> * _map_nx2;
		matrix_cursor_1<T> * _mxang;
	} __cnvrt;

	__cnvrt._void = ms.getcur();
	this->_M_map_lmb = __cnvrt._map_lmb;
	ms.seek( sizeof(size_struct<2>) + sizeof(_min1_struct)*this->_M_map_lmb->size(), seek_dir::cur );

	__cnvrt._void = ms.getcur();
	this->_M_map_nx2 = __cnvrt._map_nx2;
	ms.seek( sizeof(size_struct<3>) + sizeof(int)*this->_M_map_nx2->size(), seek_dir::cur );

	__cnvrt._void = ms.getcur();
	this->_M_node = __cnvrt._node;
	ms.seek( sizeof(size_struct<3>) + sizeof(int)*this->_M_node->size(), seek_dir::cur );

	return ms.tell() - _seek_start;
}

const typename ixs_angular_map::size_type comp_size()const
{
	size_type __size = 0;
	__size += sizeof(_lmax_struct);
	__size += sizeof(mapping_enum);
	__size += sizeof(_mappos_struct);

	__size += this->comp_lmb_size();
	__size += this->comp_nx2_size();
	__size += this->comp_node_size();
	return __size;
}
const typename ixs_angular_map::size_type comp_lmb_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<2>);
	__size += sizeof(_min1_struct) * this->map_lmb_N() * this->map_lmb_M();
	return __size;
}
const typename ixs_angular_map::size_type comp_nx2_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<3>);
	__size += sizeof(int) * this->map_nx2_N() * this->map_nx2_M() * this->map_nx2_P();
	return __size;
}
const typename ixs_angular_map::size_type comp_node_size()const
{
	size_type __size = 0;
	__size += sizeof(size_struct<3>);
	__size += sizeof(int) * this->node_N() * this->node_M() * this->node_P();
	return __size;
}

// init
// mx_lmb
void ixs_angular_map::init_map_lmb()
{
	// 2d mx_lmb
	for(int l = 0; l < this->l_max(); ++l )// n() = (l_max + lso_max + 1)
	{
		this->map2lmbA_set_l( l );
		for(int n = 0; n < this->M_map_lmb_m(); ++n )// m() = max( la_max + 1, lb_max + 1 )
		{
			this->map2lmbA_set_nx( n );
			this->map2lmbA_min() = (l - n < 0 ? (l + n)%2 : l - n);
			this->map2lmbA_size() = (l + n - this->map2lmbA_min())/2 + 1;
		}
	}
}
// mx_nxn
void ixs_angular_map::init_map_nx2()
{
	// 3d mx_nxn
	for(int l = 0; l < this->l_max(); ++l )
	{
		this->map2lmbA_set_l( l );
		this->map2lmbB_set_l( l );
		this->map3nx2_set_l( l );
		for(int na = 0; na < this->M_map_nx2_m(); ++na)
		{
			this->map2lmbA_set_nx( na );
			this->map3nx2_set_na( na );
			for(int nb = 0; nb < this->M_map_nx2_p(); ++nb )
			{
				this->map2lmbB_set_nx( nb );

				this->map3nx2_set_nb( nb );
				this->map3nx2()  = this->map2lmbA_size();
				this->map3nx2() *= this->map2lmbB_size();
			}
		}
	}
	this->map3nx2_set_lmax();
	for(int na = 0; na < this->M_map_nx2_m(); ++na)
	{
		this->map3nx2_set_na( na );
		for(int nb = 0; nb < this->M_map_nx2_p(); ++nb )
		{
			this->map3nx2_set_nb( nb );
			this->map3nx2() = (na+nb)/2 + 1;
		}
	}
}

void ixs_angular_map::init_node_min( size_type & __pos, const int & la, const int & lb )
{
	if( la%2 != lb%2 )
		return;
	// semi-local
	for(int l = la%2; l < this->l_max(); l += 2 )
	{
		this->mx3node_set_l( l );
		this->mx3node_pos() = __pos++;
		this->mx3node_size() = 1;
	}
	// local
	// not need
}
void ixs_angular_map::init_node_mid( size_type & __pos, const int & la, const int & lb )
{
	size_type __size;
	// semi-local
	for(int l = lb%2; l < this->l_max(); l+=2 )
	{
		this->mx3node_set_l( l );
		this->map3nx2_set_l( l );
		__size = 0;
		for(int na = 0; na <= la; ++na )
		{
			this->map3nx2_set_na( na );
			this->map3nx2_set_nb( lb );
			__size += this->map3nx2();
		}
		this->mx3node_pos() = __pos;
		this->mx3node_size() = __size;
		__pos += __size;
	}
	this->mx3node_set_lmax();
	this->map3nx2_set_lmax();
	__size = 0;
	// local
	for(int na = 0; na <= la; ++na)
	{
		this->map3nx2_set_na( na );
		this->map3nx2_set_nb( lb );
		__size += this->map3nx2();
	}
	__size *= alp_m.map2AB_size();
	this->mx3node_pos() = __pos;
	this->mx3node_size() = __size;
	__pos += __size;
}

void ixs_angular_map::init_node_max( size_type & __pos, const int & la, const int & lb, alpha_map & alp_m )
{
	size_type __size;
	// semi-local
	for(int l = 0; l < this->l_max(); ++l )
	{
		this->mx3node_set_l( l );
		this->map3nx2_set_l( l );
		__size = 0;
		for(int na = 0; na <= la; ++na )
		{
			this->map3nx2_set_na( na );
			for(int nb = 0; nb <= lb; ++nb)
			{
				this->map3nx2_set_nb( nb );
				__size += this->map3nx2();
			}
		}
		this->mx3node_pos() = __pos;
		this->mx3node_size() = __size;
		__pos += __size;
	}
	this->mx3node_set_lmax();
	this->map3nx2_set_lmax();
	__size = 0;
	// local
	for(int na = 0; na <= la; ++na)
	{
		this->map3nx2_set_na( na );
		for(int nb = 0; nb <= lb; ++nb)
		{
			this->map3nx2_set_nb( nb );
			__size += this->map3nx2();
		}
	}
	__size *= alp_m.map2AB_size();
	this->mx3node_pos() = __pos;
	this->mx3node_size() = __size;
	__pos += __size;
}

// mx_big
const typename ixs_angular_map::size_type ixs_angular_map::init_node_max( alpha_map & alp_m )
{
	// 3d mx_big
	size_type __size, __pos = 0;
	for(int la = 0, ila_size = 1; la <= this->la_max(); ++la, ila_size += (la + 1))
	{
		this->mx3node_set_la( la );
		alp_m.map2AB_set_la( la );
		for(int i_la = 0; i_la < ila_size; ++i_la )
		{
			this->mx3node_set_ia( i_la );
			for(int lb = 0, ilb_size = 1; lb <= this->lb_max(); ++lb, ilb_size += (lb + 1))
			{
				this->mx3node_set_lb( lb );
				alp_m.map2AB_set_lb( lb );
				for(int i_lb = 0; i_lb < ilb_size; ++i_lb )
				{
					this->mx3node_set_ib( i_lb );
					this->init_node_max( __pos, la, lb, alp_m );
				}
			}
		}
	}
	this->_mxang_size = __pos;
	return __pos;
}

const typename ixs_angular_map::size_type ixs_angular_map::init_node_mid()
{
	// 3d mx_big
	size_type __size, __pos = 0;
	for(int la = 0, ila_size = 1; la <= this->la_max(); ++la, ila_size += (la + 1))
	{
		this->mx3node_set_la( la );
		for(int i_la = 0; i_la < ila_size; ++i_la )
		{
			this->mx3node_set_ia( i_la );
			for(int lb = 0, ilb_size = 1; lb <= this->lb_max(); ++lb, ilb_size += (lb + 1))
			{
				this->mx3node_set_lb( lb );
				for(int i_lb = 0; i_lb < ilb_size; ++i_lb )
				{
					this->mx3node_set_ib( i_lb );
					this->init_node_mid( __pos, la, lb );
				}
			}
		}
	}
	this->_mxang_size = __pos;
	return __pos;
}

const typename ixs_angular_map::size_type ixs_angular_map::init_node_min()
{
	size_type __size, __pos = 0;
	for(int la = 0, ila_size = 1; la <= this->la_max(); ++la, ila_size += (la + 1))
	{
		this->mx3node_set_la( la );
		for(int i_la = 0; i_la < ila_size; ++i_la )
		{
			this->mx3node_set_ia( i_la );
			for(int lb = 0, ilb_size = 1; lb <= this->lb_max(); ++lb, ilb_size += (lb + 1))
			{
				if( lb%2 != la%2 )
					continue;
				this->mx3node_set_lb( lb );
				for(int i_lb = 0; i_lb < ilb_size; ++i_lb )
				{
					this->mx3node_set_ib( i_lb );
					this->init_node_min( __pos, la, lb );
				}
			}
		}
	}
	this->_mxang_size = __pos;
	return __pos;
}

void ixs_angular_map::init_map(alpha_map & alp_m)
{
	this->init_map_lmb();
	this->init_map_nx2();
	switch( this->get_mapping() )
	{
	case minimum :
		this->init_node_min();
		break;
	case middle  :
		this->init_node_mid();
		break;
	case maximum :
		this->init_node_max(alp_m);
		break;
	default :
		this->error("init_map", "unknowen mapping type");
		std::cerr << "current mapping type : " << this->get_mapping() << std::endl;
		exit(1);
	}
}

#endif//__IXS_ANGULAR_MAP_HPP__
