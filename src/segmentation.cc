#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/File_header_OFF.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/Cartesian.h>
#include <iostream>
#include <fstream>

/* If you do not want to use a viewer, you can comment the following file. */
#include "linear_cell_complex_3_viewer_qt.hh"



#include "unionfind.hh"



struct Sum_functor
{
  template<class Cell_attribute>
  void operator() (Cell_attribute& ca1, Cell_attribute& ca2) { ca1.info() = ca1.info() + ca2.info(); }
};
struct Divide_by_two_functor
{
  template<class Cell_attribute>
  void operator() (Cell_attribute& ca1, Cell_attribute& ca2) { ca1.info() = (ca1.info()/2); ca2.info() = (ca1.info()); }
};
struct Average_functor
{
  template<class Cell_attribute>
  void operator() (Cell_attribute& ca1, Cell_attribute& ca2) { ca1.info() = (ca1.info() + ca2.info()) / 2; }
};


struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart											< 3, Refs >																															Dart;
    typedef CGAL::Cell_attribute_with_point	< Refs, Color, CGAL::Tag_true, Average_functor > 												Vertex_attribute;
    typedef CGAL::Cell_attribute						< Refs, Color, CGAL::Tag_true, Sum_functor, Divide_by_two_functor >			Color_attribute;
    typedef CGAL::Cell_attribute						< Refs, UnionFind, CGAL::Tag_true, Sum_functor, Divide_by_two_functor >	UnionFind_attribute;
	//typedef CGAL::cpp11::tuple							< Vertex_attribute, Color_attribute, Color_attribute >									Attributes;
    typedef CGAL::cpp11::tuple							< Vertex_attribute, Color_attribute, UnionFind_attribute >							Attributes;
  };
};

/* ==========================================================================
   ========================================================================== */
typedef CGAL::Cartesian									< double >             		Kernel;
typedef CGAL::Linear_cell_complex_traits< 3, Kernel >         		Traits;
typedef CGAL::Linear_cell_complex				< 3, 3, Traits, Myitem >	LCC_3;
/* ==========================================================================
   ========================================================================== */





typedef enum { LOCAL, GLOBAL } Sampler;


int main(int argc, char* argv[])
{
	srand(time(NULL));
	
	float		threshold = 0.75f;
	Sampler	method		= GLOBAL;
	
	if ( argc < 2 || !strcmp(argv[1],"-h") || !strcmp(argv[1],"-?") )
	{
		std::cout	<< "Usage : segmentation filename [options]"															<< std::endl
							<< "   filename being an '.off' file."																		<< std::endl
							<< "Available options :"																									<< std::endl
							<< " -s <float>        select the threshold for variation [-1, 1]"				<< std::endl
							<< "                   (default 0.75)"																		<< std::endl
							<< "                   the higher the threshold the smaller the segments"	<< std::endl
							<< " -m [local|global] select segmentation method"												<< std::endl
							<< "                   (default global)"																	<< std::endl;
		return EXIT_FAILURE;
	}
	
	for (int i = 2; i<argc; ++i)
	{
		if (argc > i+1 && !strcmp(argv[i], "-s"))
			threshold = atof(argv[++i]);
		else if (argc > i+1 && !strcmp(argv[i], "-m"))
		{
			++i;
			if (!strcmp(argv[i], "local"))
				method = LOCAL;
			else if (!strcmp(argv[i], "global"))
				method = GLOBAL;
			else
				printf("Warning : Unknown segmentation method : %s\n", argv[i]);
		}
		else
			printf("Warning : Unknown option : %s\n", argv[i]);
	}
			
	std::ifstream ifs(argv[1]);
	if ( ifs.fail() )
	{
		std::cout << "Error : impossible to open file " << argv[1] << std::endl;
		return EXIT_FAILURE;
	}
	
	
	LCC_3							lcc;
	Geom_utils<LCC_3>	geomutils;
	CGAL::load_off(lcc, ifs);

	
	LCC_3::Dart_range::iterator it_begin = lcc.darts().begin();
	LCC_3::Dart_range::iterator it_end   = lcc.darts().end();
	
	for ( LCC_3::Dart_range::iterator it = it_begin; it != it_end; ++it )
	{
		it->attribute<0>()->info().r = 1;
		it->attribute<0>()->info().g = 0;
		it->attribute<0>()->info().b = 0;
		
		if (it->attribute<1>() == NULL)
		{
			lcc.set_attribute<1>(it, lcc.create_attribute<1>());
			it->attribute<1>()->info().r = 0.5;
			it->attribute<1>()->info().g = 0;
			it->attribute<1>()->info().b = 0.5;
		}
		
		if (it->attribute<2>() == NULL)
		{
			lcc.set_attribute<2>(it, lcc.create_attribute<2>());
			it->attribute<2>()->info().normal() = geomutils.get_facet_normal(lcc, it);
		}
	}

	int nb_cells = lcc.one_dart_per_cell<2>().size();
	printf("Nb Cells %d\n", nb_cells);
	
	
	/* ======================================================================== *
	 * SEGMENTATION																															*
	 * ======================================================================== */
	std::cout << "Computing Segmentation ... ";
	fflush(stdout);
	
	switch (method)

	{
		case LOCAL:
		{
			for ( LCC_3::Dart_range::iterator it1 = it_begin; it1 != it_end; ++it1 )
			{
				Local_vector n1 = geomutils.get_facet_normal(lcc, it1);
				for ( LCC_3::Dart_of_orbit_range<2>::iterator it2 = lcc.darts_of_orbit<2>(it1).begin(); it2.cont(); ++it2 )
				{
					Local_vector n2 = geomutils.get_facet_normal(lcc, it2);
					if (n1 * n2 >= threshold)
						UnionFind::merge(it1->attribute<2>()->info(), it2->attribute<2>()->info());
				}
			}
			break;
		}
		
		case GLOBAL:
		{
			int treated = lcc.get_new_mark();
			for ( LCC_3::Dart_range::iterator it = it_begin; it != it_end; ++it )
				{
				if ( !lcc.is_marked(it, treated) )
				{
					lcc.mark(it, treated);		
				
					std::list<LCC_3::Dart_of_orbit_range<2>::iterator> queue;
					for (	LCC_3::Dart_of_orbit_range<2>::iterator it1 = lcc.darts_of_orbit<2>(it).begin(); it1.cont(); ++it1 )
						queue.push_back(it1);
			
					while (!queue.empty())
					{
						LCC_3::Dart_of_orbit_range<2>::iterator it1 = queue.front();
						queue.pop_front();
				
						Local_vector n1 = it->attribute<2>()->info().root().normal();
						Local_vector n2 = it1->attribute<2>()->info().root().normal();
				
						if (n1 * n2 >= threshold)
						{
							UnionFind::merge(it->attribute<2>()->info(), it1->attribute<2>()->info(), true);
							if (!lcc.is_marked(it1, treated))
							{
								for (	LCC_3::Dart_of_orbit_range<2>::iterator it2 = lcc.darts_of_orbit<2>(it1).begin(); it2.cont(); ++it2 )
									queue.push_back(it2);
								lcc.mark(it1, treated);
							}
						}
					}
				}
			}
			CGAL_assertion(lcc.is_whole_map_marked(treated));
			lcc.free_mark(treated);  
			break;
		}
	}
	
	std::cout << "done" << std::endl;
	/* ======================================================================== */
    
	lcc.display_characteristics(std::cout) << ", valid=" << lcc.is_valid() << std::endl;
	display_lcc(lcc);

	return EXIT_SUCCESS;
}