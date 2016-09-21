#include"geom.slm.h"
#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

template<typename T>
int demo_geom_slm(char const * file )
{
	ifstream inp( file );
	if( !inp.is_open() )
	{
		cerr << "Error: [demo_geom_slm] can't open file '" << file << "'" << endl;
		return 1;
	}
	int l_max, lso_max, la_max, lb_max;
	inp >> l_max >> lso_max >> la_max >> lb_max;

	T CA[3], CB[3];
	for(int i = 0; i < 3; ++i) inp >> CA[i];
	for(int i = 0; i < 3; ++i) inp >> CB[i];

	geom_slm<T> gs;
	gs.set_mapping_max();

	gs.set_l_max( l_max );
	gs.set_lso_max( lso_max );
	gs.set_la_max( la_max );
	gs.set_lb_max( lb_max );
	gs.init_lsize();

	gs.write();

	matrix_slm<T> mx_slm;
	const char file_slm[] = "vector.slm.32.src";
	mx_slm.read( file_slm );

	gs.init( CA, CB, mx_slm );
	return 0;
}

int main(int argc, char ** argv)
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './main.exe file.inp'" << endl;
		return 1;
	}
	return demo_geom_slm<double>( argv[1] );
}
