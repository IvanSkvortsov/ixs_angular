#include"ixs.angular.mem.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include"get.work.time.h"

using namespace std;

void print_i( const int value, const char * name )
{
	cout << setw(14) << name << " : " << setw(14) << value << endl;
}

void print_p( void const * _pointer, const char * name )
{
	cout << setw(14) << name << " : [" << _pointer << "]" << endl;
}

inline void timer_start( get_work_time & t ){ t.counting_start();}
inline void timer_stop( get_work_time & t, const char name[] = "" )
{
	t.counting_end();
	fprintf(stdout, "%24s time : %14.4f sec, %12llu nsec\n", name, t.get_realseconds(), t.get_nanoseconds() );
}

int demo_ang_map(const char * file )
{
	union{
		typename ixs_angular_map::_lmax_struct _lmax_ang;
		typename alpha_map::_lmax_struct _lmax_alp;
		typename ixs_angular_map::_lmax_struct _lmax;
	};
	get_work_time t0;
	timer_start( t0 );
	//_lmax_struct _lmax;
	_lmax._l_max = 4;
	_lmax._lso_max = 4;
	_lmax._la_max = 4;
	_lmax._lb_max = 4;

	const int Asize[] = {15, 10,  5,  3, 1, 1, 1, 1, 1, 1, 1};
	const int Bsize[] = {17, 13,  8,  5, 3, 1, 1, 1, 1, 1, 1};
	const int Csize[] = { 8,  5,  3,  2, 7, 5, 4, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1};

	mapping_struct map_t;
	map_t.set_mapping_max();

	ixs_angular_mem ixs_am;
	ixs_am.set_lmax( _lmax );
	ixs_am.set_mapping( map_t );

	memory_map alp_mmap;
	memorystream ms;
	alpha_map alp_m;
	alp_m.set_mapping( map_t );

	get_work_time t1;
	timer_start( t1 );
	alp_mmap.create( "alpha.map.tmp", alp_m.comp_size(_lmax_alp) );
	timer_stop( t1, "alpha map create" );

	timer_start( t1 );
	ms.setbuf( alp_mmap.data(), alp_mmap.size() );
	alp_m.write_map( ms, _lmax_alp );
	alp_m.init_map( Asize, Bsize, Csize );
	timer_stop( t1, "alpha map init" );

	print_i( ixs_am.comp_size(), "comp_size" );

	timer_start( t1 );
	ixs_am.memory_create( file );
	timer_stop( t1, "ixs_ang_map create" );

	print_i( ixs_am.size(), "map_size" );
	print_p( ixs_am.data(), "map_data" );
	if( ixs_am.data() )
		print_p( (char *)ixs_am.data() + ixs_am.size(), "map_end" );

	timer_start( t1 );
	ixs_am.init_map(alp_m);
	timer_stop( t1, "ixs_ang_map init" );
	timer_stop( t0, "total" );
	return 0;
}

int main(int argc, char ** argv )
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './main.exe file.inp'" << endl;
		return 1;
	}
	return demo_ang_map( argv[1] );
}
