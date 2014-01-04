#ifndef UNIONFIND_HH
#define UNIONFIND_HH

class UnionFind
{
	private:
		int					_id;
		UnionFind*	_ft;
	
	public:
		UnionFind();
		int&				id();
		UnionFind& 	root();
		static void merge(UnionFind& uf1, UnionFind& uf2);	
};






UnionFind::UnionFind()
{
	static int cnt = 0;
	_id	= cnt++;
	_ft	= this;
}
int& UnionFind::id()
{
	return _id;
}
UnionFind& UnionFind::root()
{
	if (_ft != this)
		_ft = &(_ft->root());
	return *_ft;
}
void UnionFind::merge(UnionFind& uf1, UnionFind& uf2)
{
	uf1.root()._ft = uf2.root()._ft;
}



#endif