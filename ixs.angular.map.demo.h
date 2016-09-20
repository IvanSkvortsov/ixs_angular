#include"ixs.angular.map.h"
#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

void print_i( const int value, const char * name )
{
	cout << setw(14) << name << " : " << setw(14) << value << endl;
}

void print_p( void const * _pointer, const char * name )
{
	cout << setw(14) << name << " : [" << _pointer << "]" << endl;
}

int demo_ang_map(const char * file )
{
	typedef typename ixs_angular_map::_lmax_struct _lmax_struct;
	_lmax_struct _lmax;
	_lmax._l_max = 4;
	_lmax._lso_max = 4;
	_lmax._la_max = 4;
	_lmax._lb_max = 4;

	ixs_angular_map ixs_am;
	ixs_am.set_lmax( _lmax );
	print_i( ixs_am.comp_size(), "comp_size" );
	ixs_am.memory_create( file );
	print_i( ixs_am.size(), "map_size" );
	print_p( ixs_am.data(), "map_data" );
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
