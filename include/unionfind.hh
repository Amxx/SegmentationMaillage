#ifndef UNIONFIND_HH
#define UNIONFIND_HH



#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<double>						Local_kernel;
typedef typename Local_kernel::Point_3		Local_point;
typedef typename Local_kernel::Vector_3		Local_vector;




class UnionFind
{
	private:
		int						_id;
		UnionFind*		_ft;
	
		Local_vector	_norm;
		double				_norm_w;
	
	public:
		UnionFind();
		int&					id();
		Local_vector&	normal();
	
		UnionFind& 		root();
		static void 	merge(UnionFind& uf1, UnionFind& uf2, bool = false);
};






UnionFind::UnionFind()
{
	static int cnt = 0;
	_id			= cnt++;
	_ft			= this;
	_norm		= Local_vector(0., 0., 0.);
	_norm_w	= 1.0;
}
int& UnionFind::id()
{
	return _id;
}
Local_vector&	UnionFind::normal()
{
	return _norm;
}
UnionFind& UnionFind::root()
{
	if (_ft != this)
		_ft = &(_ft->root());
	return *_ft;
}
void UnionFind::merge(UnionFind& uf1, UnionFind& uf2, bool norm)
{
	UnionFind& r1 = uf1.root();
	UnionFind& r2 = uf2.root();
	
	if (&r1 == &r2) return;
	
	r2._ft	= r1._ft;
	
	if (norm)
	{
		Local_vector med = Local_vector(r1._norm_w * r1._norm.x() + r2._norm_w * r2._norm.x(),
																		r1._norm_w * r1._norm.y() + r2._norm_w * r2._norm.y(),
																		r1._norm_w * r1._norm.z() + r2._norm_w * r2._norm.z()	);
		
		r1._norm		= med / CGAL::sqrt(med*med);
		r1._norm_w	= r1._norm_w + r2._norm_w;
	}
}


#endif