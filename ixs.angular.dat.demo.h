#include"ixs.angular.dat.h"
#include"ixs.angular.mem.h"
#include"alpha.slm.h"
#include"alpha.mem.h"
//#include"geom.slm.h"
#include<cmath>
#include<iostream>
#include<fstream>
#include<iomanip>

#include"mpreal.h"

using namespace std;


template<typename T>
void print_lmax( std::ostream & out, T const & v )
{
	out << "l_max    : " << setw(4) << v.get_l_max() << endl;
	out << "lso_max  : " << setw(4) << v.get_lso_max() << endl;
	out << "la_max   : " << setw(4) << v.get_la_max() << endl;
	out << "lb_max   : " << setw(4) << v.get_lb_max() << endl;
	out << endl;
}

void print_ixs( std::ostream & out, ixs_angular_mem const & ixs_mem )
{
	out << "map.data : [" << ixs_mem.data() << "]" << endl;
	out << "map.size : " << ixs_mem.size() << endl;
	out << "map.end  : [" << (void *)((char *)ixs_mem.data() + ixs_mem.size()) << "]" << endl;
	out << endl;
	out << "map_lmb.n : " << setw(8) << ixs_mem.M_map_lmb_n() << endl;
	out << "map_lmb.m : " << setw(8) << ixs_mem.M_map_lmb_m() << endl;
	out << "map_lmb.size" << setw(8) << ixs_mem.M_map_lmb_size() << endl;
	out << endl;
	out << "map_nx2.n : " << setw(8) << ixs_mem.M_map_nx2_n() << endl;
	out << "map_nx2.m : " << setw(8) << ixs_mem.M_map_nx2_m() << endl;
	out << "map_nx2.p : " << setw(8) << ixs_mem.M_map_nx2_p() << endl;
	out << "map_nx2.size" << setw(8) << ixs_mem.M_map_nx2_size() << endl;
	out << endl;
	out << "node.n : " << setw(8) << ixs_mem.M_node_n() << endl;
	out << "node.m : " << setw(8) << ixs_mem.M_node_m() << endl;
	out << "node.p : " << setw(8) << ixs_mem.M_node_p() << endl;
	out << "node.size" << setw(8) << ixs_mem.M_node_size() << endl;
	out << endl;
}

template<typename T, typename U>
int demo_angular(const char * file_inp)
{
	ifstream inp( file_inp );
	if( !inp.is_open() )
	{
		cerr << "Error: [demo_angular] can't open file '" << file_inp << '\'' << endl;
		return 1;
	}

	typename alpha_siz::_lmax_struct _lmax;
	inp >> _lmax._l_max >> _lmax._lso_max >> _lmax._la_max >> _lmax._lb_max;

	mapping_struct map_t;
	map_t.set_mapping_min();

	alpha_mem<T,U> alp_mem;
	alp_mem.set_lmax( _lmax );
	alp_mem.set_mapping( map_t );

	std::ostream & out = std::clog;
	print_lmax( out, alp_mem );

	const char file[] = "ixs.angular.mem.mmap";
	out << "ixs_map.file '" << file << "'" << endl;
	ixs_angular_mem ixs_mem;
	ixs_mem.ixs_angular_map::set_lmax( _lmax );
	ixs_mem.set_mapping( map_t );

	ixs_mem.memory_create( file );
	if( ixs_mem.is_mapping_max() )
		ixs_mem.init_map( alp_mem );
	else
		ixs_mem.init_map();

	print_ixs( out, ixs_mem );

	ixs_mem.close();
	ixs_angular_mem ixs_mem2;
	cout << "opening..." << endl;
	ixs_mem2.memory_open( file );
	print_ixs( cout, ixs_mem2 );

	ixs_angular_dat<T,U> ixs_dat;
	ixs_dat.set_angular_map( ixs_mem2 );
	//ixs_dat.ixs_angular_map::operator=( ixs_mem2 );// - same as %ixs_angular_dat<T,U>::set_angular_map( ixs_angular_map const & )

	ixs_dat.memory_create("ixs.angular.dat.mmap");
	cout << "map.data   : [" << ixs_dat.data() << "]" << endl;
	cout << "map.size   : " << setw(8) << ixs_dat.size() << endl;
	cout << "array.size : " << setw(8) << ixs_dat.M_mx1ang_size() << endl;
	cout << "array.data : [" << ixs_dat.M_mx1ang_data() << "]" << endl;

	ixs_omega<T,U> ixs_omg;
	ixs_omg.memory_open( "ixs.omega.06.mmap" );

	cout << "angular integral values computation:" << endl << flush;
	cout << "start..." << endl << flush;
	cerr << "start..." << endl;

	if( ixs_dat.is_mapping_max() )
		ixs_dat.comp_ang_max( geom_s, alp_s, ixs_omg );
	else if( ixs_dat.is_mapping_mid() )
		ixs_dat.comp_ang_mid( geom_s, ixs_omg );
	else if( ixs_dat.is_mapping_min() )
		ixs_dat.comp_ang_min( ixs_omg );
	cout << "end" << endl << endl << flush;
	cerr << "end" << endl;
	return 0;
}

int demo_ixsang_d(const char * file)
{
	return demo_ixsang<double,long double>(file);
}
int main(int argc, char ** argv)
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './main.exe file.inp'" << endl;
		return 1;
	}
	return demo_ixsang_d(argv[1]);
}

#include"alpha.mem.h"

#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cmath>
#include"sf.math.h"
#include"mpreal.h"
#include"vector.3d.h"

using namespace std;

inline void print_i( char const * name, const int value)
{
	cout << setw(14) << name << " : " << setw(10) << value << endl;
}

inline const int read_lxsize( const char * name, istream & inp, int * _lxsize, const int _size )
{
	int __size = 0, _num;
	inp >> _num;
	cout << setw(10) << name << " : ";
	for(int i = 0; i < _size || i < _num; ++i)
	{
		inp >> _lxsize[i];
		if( i < _size )
		{
			__size += _lxsize[i];
			cout << setw(4) << _lxsize[i];
		}
	}
	cout << endl;
	cout << setw(4) << __size << endl << endl;
	return __size;
}
template<typename T>
inline void read_alp( const char * name, istream & inp, T * _vec, const int _size )
{
	cout << setw(10) << name << " :" << endl;
	for(int i = 0; i < _size; ++i)
	{
		inp >> _vec[i];
		if( i%8 == 1 )
			cout << setw(4) << i << " : " << setw(24) << setprecision(14) << scientific << _vec[i] << endl;
	}
	cout << endl;
}
template<typename T, typename U>
int demo_mem(const char * file)
{
	typedef typename alpha_map::_lmax_struct _lmax_struct;
	typedef typename alpha_map::_alpsize_struct _alpsize_struct;
	ifstream inp( file );
	if( !inp.is_open() )
	{
		cerr << "Error: [demo_val] can't open file '" << file << "'" << endl;
		return 1;
	}

	_lmax_struct _lmax;
	inp >> _lmax._la_max >> _lmax._lb_max >> _lmax._l_max >> _lmax._lso_max;
	print_i( "la_max", _lmax._la_max );
	print_i( "lb_max", _lmax._lb_max );
	print_i( "l_max", _lmax._l_max );
	print_i( "lso_max", _lmax._lso_max );
	cout << endl;

	int _Asize[20], _Bsize[20], _Csize[20];
	int _a_size, _b_size, _c_size;
	_a_size = read_lxsize( "alpha A", inp, _Asize, _lmax._la_max + 1 );
	_b_size = read_lxsize( "alpha B", inp, _Bsize, _lmax._lb_max + 1 );
	_c_size = read_lxsize( "alpha C", inp, _Csize, _lmax._l_max + _lmax._lso_max + 1 );

	vector_3d<U> CA, CB;
	U CACB;
	for(int i = 0; i < 3; ++i)
	{
		inp >> CA[i];
		cout << setw(3) << i << " : " << setw(24) << setprecision(16) << scientific << CA[i] << endl;
	}
	for(int i = 0; i < 3; ++i)
	{
		inp >> CB[i];
		cout << setw(3) << i << " : " << setw(24) << setprecision(16) << scientific << CB[i] << endl;
	}
	CACB = CA.scalar( CB );

	vector<U> _Avec( _a_size ), _Bvec( _b_size ), _Cvec( _c_size );
	read_alp( "alpha A", inp, _Avec.data(), _a_size );
	read_alp( "alpha B", inp, _Bvec.data(), _b_size );
	read_alp( "alpha C", inp, _Cvec.data(), _c_size );

	alpha_mem<T,U> m_alp;
	m_alp.set_mapping_max();
	_alpsize_struct _alpsize = {._a_size = _a_size, _b_size, _c_size};
	m_alp.set_lmax( _lmax );
	m_alp.set_alpsize( _alpsize );
	m_alp.set_A_size( _Asize );
	m_alp.set_B_size( _Bsize );
	m_alp.set_C_size( _Csize );
	cout << "m_alp.size : " << setw(10) << m_alp.comp_size( _lmax, _alpsize, _Asize, _Bsize ) << endl << flush;
	if( m_alp.is_mapping_min() )
		cout << "mapping : " << setw(10) << "minimum" << endl;
	else if( m_alp.is_mapping_mid() )
		cout << "mapping : " << setw(10) << "middle" << endl;
	else if( m_alp.is_mapping_max() )
		cout << "mapping : " << setw(10) << "maximum" << endl;
	else
		cout << "mapping : " << setw(10) << "undefined" << endl;

	cout << "sizeof class : " << setw(4) << sizeof(m_alp) << endl;
	m_alp.memory_create( "alpha.mem.map" );

	m_alp.init_map();// init_map
	// init_mxalp
	if( m_alp.is_mapping_min() )
		m_alp.init_mxalp( _Avec.data(), _Bvec.data(), _Cvec.data() );
	else if( m_alp.is_mapping_mid() )
		m_alp.init_mxalp( _Avec.data(), _Bvec.data(), _Cvec.data(), CA.len() );
	else if( m_alp.is_mapping_max() )
		m_alp.init_mxalp( _Avec.data(), _Bvec.data(), _Cvec.data(), CA.len(), CB.len(), CACB );
	else
	{
		cerr << "Error: [demo_val] mapping undefined" << endl;
		return 2;
	}
	m_alp.init_mxpow(CA.data(), _Avec.data(), CB.data(), _Bvec.data() );// init_mxpow

	cout << endl;
	print_i( "mx1A.size", m_alp.mx1A_size() );
	print_i( "mx1B.size", m_alp.mx1B_size() );
	print_i( "mx1C.size", m_alp.mx1C_size() );
	if( m_alp.is_mapping_mid() || m_alp.is_mapping_max() )
	{
		print_i( "mx1kA.size", m_alp.mx1kA_size() );
		//print_i( "mx1powkA.size", m_alp.mx1powkA_size() );
	}
	if( m_alp.is_mapping_max() )
	{
		print_i( "mx1kB.size", m_alp.mx1kB_size() );
		print_i( "mx2k.size", m_alp.mx2k_size() );
		print_i( "mx2slm.size", m_alp.mx2slm_size() );
		//print_i( "mx1powkB.size", m_alp.mx1powkB_size() );
		//print_i( "mx2powk.size", m_alp.mx2powk_size() );
	}
	if( m_alp.is_mapping_min() )
		return 0;
	for(int i = 0; i < 52; ++i) cout << '-'; cout << endl;
	cout << setw(4) << "i" << setw(24) << "alpA" << setw(24) << "kA" << endl;
	U ta, tb, tc, _k;
	for(int la = 0, lb = 2; la < m_alp.la_size(); ++la)
	{
		m_alp.map1A_set_l( la );
		m_alp.map2AB_set_la( la );
		m_alp.map2AB_set_lb( lb );
		m_alp.map1B_set_l( lb );
		m_alp.map2powk_set_la( la );
		m_alp.map2powk_set_lb( lb );
		tb = CB.len() * _Bvec[m_alp.map1B_pos()] * 2;
		tb *= tb;
		for(int it = 0, i = m_alp.map1A_pos(); it < m_alp.map1A_size(); ++i, ++it )
		{
			if( i % 6 != 1 && i > 2) continue;
			ta = CA.len() * _Avec[i] * 2;
			ta *= ta;
			tc = 8 * _Avec[i] * _Bvec[m_alp.map1B_pos()] * CACB;
			_k = ta + tc + tb;
			_k = sqrt( _k );
			cout << setw(4) << i <<
				//setw(24) << setprecision(16) << scientific << m_alp.mx1A( i ) <<
				setw(24) << setprecision(16) << scientific << m_alp.mx1kA( la, it ) <<
				setw(24) << setprecision(16) << scientific << m_alp.mx1powkA( la, it)[2] <<
				setw(24) << setprecision(16) << scientific << m_alp.mx1kA( la, it ) * m_alp.mx1kA( la, it );
			if( m_alp.is_mapping_max() )
			{
				cout <<
				//setw(24) << setprecision(16) << scientific << m_alp.mx2k( la, it, lb, 0 ) <<
				setw(24) << setprecision(16) << scientific << m_alp.mx2k( it, 0 ) <<

				//setw(24) << setprecision(16) << scientific << m_alp.mx2powk( it, 0 )[2] <<
				setw(24) << setprecision(16) << scientific << m_alp.mx2powk( la, it, lb, 0 )[2] <<
				//" [" << m_alp.mx2powk( la, it, lb, 0 ) + 1 << "]";
				setw(24) << setprecision(16) << scientific << m_alp.mx2k( it, 0) * m_alp.mx2k( it, 0);
			}
			cout << endl;
		}
	}
	cout << endl;
	return 0;
}

int demo_mem_d(const char * file)
{
	//mpfr::mpreal::set_default_prec( mpfr::digits2bits(40) );
	//return demo_val<double, mpfr::mpreal>( file );
	return demo_mem<double, long double>( file );
}

int main(int argc, char ** argv )
{
	if( argc != 2 )
	{
		cerr << "Error: [main] usage './main.exe file.inp'" << endl;
		return 1;
	}
	return demo_mem_d( argv[1] );
}
