#ifndef __IXS_ANGULAR_MEM_H__
#define __IXS_ANGULAR_MEM_H__
#include"memory.map.h"
#include"ixs.angular.map.h"

#define __IXS_ANGULAR_MEM_DEBUG
#ifdef  __IXS_ANGULAR_MEM_DEBUG
  #include<cassert>
#endif

struct ixs_angular_mem : public memory_map, public memorystream, public ixs_angular_map
{
protected:
	inline void error( const char * _method, const char _message[] = "nothing to do here" )const
	{ std::cerr << "Error: [" << this << "] ixs_angular_mem::" << _method << ", " << _message << std::endl; }
public:
	typedef typename memorystream::pos_type pos_type;
	typedef typename memorystream::off_type off_type;
	typedef typename memorystream::seek_dir seek_dir;

	typedef typename size_struct<1>::size_type size_type;

	inline ixs_angular_mem() : memory_map(), memorystream(), ixs_angular_map(){}
	inline ixs_angular_mem( ixs_angular_mem const & v) : memory_map(), memorystream( v ), ixs_angular_map( v ){}
	inline ixs_angular_mem & operator=(ixs_angular_mem const & v)
	{
		if( this == &v )
			return *this;
		if( this->data() || this->size() )
			this->close();
		this->memorystream::operator=( v );
		this->ixs_angular_map::operator=( v );
		return *this;
	}

	inline void sync_stream(){ this->memorystream::setbuf( this->memory_map::data(), this->memory_map::size() ); }
	inline void memory_create( const char * file, int __flags = O_RDWR| O_TRUNC| O_CREAT, mode_type __mode = memory_map::MODE_644 )
	{
		static const int _AUG_BYTES = 1024;
		size_type __size = this->comp_size();
		int _st = this->memory_map::create( file, __size + _AUG_BYTES, __flags, __mode );
		if( _st )
		{
			this->error("memory_create", "from memory_map");
			std::cerr << "file      '" << file << "'" << std::endl;
			std::cerr << "  __size : " << __size << std::endl;
			std::cerr << "map.size : " << this->size() << std::endl;
			std::cerr << "map.data : [" << this->data() << "]" << std::endl;
			exit(1);
		}
		this->sync_stream();
		memorystream & ms = *this;
		size_type __written_size = this->write_map( ms );
#ifdef  __IXS_ANGULAR_MEM_DEBUG
		assert( __size == __written_size );
#endif
	}
	inline void memory_open( const char * file, int __flags = O_RDONLY, mode_type __mode = memory_map::MODE_444 )
	{
		if( this->data() || this->size() )
				this->close();
		int _st = this->memory_map::open( file, __flags, __mode );
		if( _st )
		{
			this->error("memory_open", "from memory_map");
			std::cerr << "file      '" << file << "'" << std::endl;
			std::cerr << "map.size : " << this->size() << std::endl;
			std::cerr << "map.data : [" << this->data() << "]" << std::endl;
			exit(1);
		}
		this->sync_stream();
		memorystream & ms = *this;
		size_type __read_size = this->read_map( ms );
#ifdef  __IXS_ANGULAR_MEM_DEBUG
		assert( __read_size == this->comp_size()  );
#endif
	}
};

#endif//__IXS_ANGULAR_MEM_H__
