#include"geom.slm.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include"vector.3d.h"
#include<cassert>
#include"mpreal.h"

using namespace std;

template<typename T>
T pown_f( T const & v, int n )
{
	static T u, res;
	if( n < 0 )
		return 1/pown_f<T>( v, -n );
	res = ( n&1 ? v : 1 );
	n >>= 1;
	u = v;
	for(int i = 1; i < 31 && n; ++i, n >>= 1 )
	{
		u *= u;
		if( n & 1 )
			res *= u;
	}
	return res;
}

template<typename T>
T pown_f_( T const & v, int n )
{
	static T res;
	if( n < 0 )
		return 1/pown_f_<T>( v, -n );
	res = 1;
	for(int i = 0; i < n; ++i, res *= v);
	return res;
}

template<typename T>
void unit_test_CX( T const & gs_powCX, T const & powCX, const int i, const char * name )
{
	if( gs_powCX == powCX )
		return;
	static T err_f;
	static int it = 0;
	err_f = gs_powCX - powCX;
	err_f = ( err_f < 0 ? -err_f : err_f );
	if( it < i )
	{
		cerr << endl;
		++it;
	}
	cerr.setf( ios::scientific );
	static const int prec = 17, w = prec + 9;
	cerr << "Error " << setw(4) << i << " : " <<
		setw(w) << setprecision(prec) << gs_powCX <<
		" [" << &gs_powCX << "] " <<
		setw(w) << setprecision(prec) << powCX <<
		setw(14) << setprecision(2) << err_f <<
		setw(8) << name << endl;
	//fprintf(stderr, "Error %4d : %26.17e [%p] %26.17e%14.2e %8s\n", i, gs_powCX, &gs_powCX, powCX, err_f, name);
}
template<typename T>
void gs_test_powCX( geom_slm<T> const & gs, T const * __CA, T const * __CB )
{
	vector_3d<T> CA( __CA ), CB( __CB );
	for(int i = 0; i < gs.la_size(); ++i)
	{
		unit_test_CX( gs.CA_x(i), pown_f( CA.x(), i ), i, "CA_x");
		unit_test_CX( gs.CA_y(i), pown_f( CA.y(), i ), i, "CA_y");
		unit_test_CX( gs.CA_z(i), pown_f( CA.z(), i ), i, "CA_z");

		unit_test_CX( gs.CB_x(i), pown_f( CB.x(), i ), i, "CB_x");
		unit_test_CX( gs.CB_y(i), pown_f( CB.y(), i ), i, "CB_y");
		unit_test_CX( gs.CB_z(i), pown_f( CB.z(), i ), i, "CB_z");
	}
}

template<typename T>
inline void print_3d( T const * CX, const char * name )
{
	cout << setw(10) << name << " : ";
	for(int i = 0; i < 3; ++i)
		cout << setw(20) << CX[i];
	cout << endl;
}

inline void print_i(const int i, const char * name )
{
	cout << setw(10) << name << " : " << setw(6) << i << endl;
}
inline void print_lmax(const int l_max, const int lso_max, const int la_max, const int lb_max, const char * name )
{
	cout << setw(10) << name << " : " << endl;
	print_i( l_max, "l_max" );
	print_i( lso_max, "lso_max" );
	print_i( la_max, "la_max" );
	print_i( lb_max, "lb_max" );
	cout << endl;
}

template<typename T>
inline void print_lsize( const geom_slm<T> & gs )
{
	cout << setw(10) << "gs lsize" << " : " << endl;
	print_i( gs.la_size(), "la_size" );
	print_i( gs.lb_size(), "lb_size" );
	print_i( gs.kA_lsize(), "kA_lsize" );
	print_i( gs.kB_lsize(), "kB_lsize" );
	cout << endl;
}

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
	print_lmax( l_max, lso_max, la_max, lb_max, "input lmax");

	T CA[3], CB[3];
	for(int i = 0; i < 3; ++i) inp >> CA[i];
	for(int i = 0; i < 3; ++i) inp >> CB[i];
	print_3d( CA, "CA" );
	print_3d( CB, "CB" );
	cout << endl;

	geom_slm<T> gs;
	gs.set_mapping_max();

	gs.set_l_max( l_max );
	gs.set_lso_max( lso_max );
	gs.set_la_max( la_max );
	gs.set_lb_max( lb_max );
	gs.init_lsize();

	print_lmax( gs.l_max(), gs.lso_max(), gs.la_max(), gs.lb_max(), "gs lmax");
	assert( gs.l_max() == l_max );
	assert( gs.lso_max() == lso_max );
	assert( gs.la_max() == la_max );
	assert( gs.lb_max() == lb_max );
	print_lsize( gs );

	gs.write();
	cout << setw(10) << "gs.size" << " : " << setw(6) << gs.size() << endl;
	cout << setw(10) << "gs.data" << " : [" << setw(6) << gs.data() << "]" << endl;
	if( gs.data() )
		cout << setw(10) << "gs.end" << " : [" << setw(6) << gs.data() + gs.size() << "]" << endl;
	cout << endl;

	matrix_slm<T> mx_slm;
	const char file_slm[] = "vector.slm.32.src";
	mx_slm.read( file_slm );

	gs.init( CA, CB, mx_slm );

	assert( gs._CA.x == gs.data() );
	gs_test_powCX( gs, CA, CB );

	return 0;
}

int main(int argc, char ** argv)
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './main.exe file.inp'" << endl;
		return 1;
	}
	mpfr::mpreal::set_default_prec( mpfr::digits2bits(30) );
	return demo_geom_slm<mpfr::mpreal>( argv[1] );
	//return demo_geom_slm<double>( argv[1] );
}
