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

struct Sum_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  { ca1.info()=ca1.info()+ca2.info(); }
};
struct Divide_by_two_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  {
    ca1.info()=(ca1.info()/2);
    ca2.info()=(ca1.info());
  }
};
struct Average_functor
{
  template<class CellAttribute>
  void operator()(CellAttribute& ca1, CellAttribute& ca2)
  { ca1.info()=(ca1.info()+ ca2.info())/2; }
};

struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, Refs> Dart;
    
    typedef CGAL::Cell_attribute_with_point< Refs, Color, CGAL::Tag_true,
                                             Average_functor >
    Vertex_attribute;
    typedef CGAL::Cell_attribute<Refs, Color, CGAL::Tag_true,
                                 Sum_functor, Divide_by_two_functor>
    Color_attribute;

    typedef CGAL::cpp11::tuple<Vertex_attribute,Color_attribute,Color_attribute> Attributes;
  };
};

/* ==========================================================================
   ========================================================================== */
typedef CGAL::Cartesian<double>                             Kernel;
typedef CGAL::Linear_cell_complex_traits<3, Kernel>         Traits;
typedef CGAL::Linear_cell_complex<3, 3, Traits, Myitem>     LCC_3;

/* ==========================================================================
   ========================================================================== */

int main(int narg, char** argv)
{
    srand(time(NULL));
    
    if (narg<2 || (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
    {
        std::cout<<"Usage : load_off filename"<<std::endl
                    <<"   filename being an '.off' file."<<std::endl;
        return EXIT_FAILURE;
    }
    std::ifstream ifs(argv[1]);
    if ( ifs.fail() )
    {
        std::cout<<"Error : impossible to open file "<<argv[1]<<std::endl;
        return EXIT_FAILURE;
    }

    LCC_3 lcc;
    CGAL::load_off(lcc, ifs);


    //INIT
    for(typename LCC_3::Dart_range::iterator it=lcc.darts().begin(),
        itend=lcc.darts().end(); it!=itend; ++it)
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
            it->attribute<2>()->info().r = 0;
            it->attribute<2>()->info().g = 0;
            it->attribute<2>()->info().b = 0;
        }
    }

    int nb_cells = lcc.one_dart_per_cell<2>().size();
    printf("Nb Cells %d\n", nb_cells);
    int n_seed = 3;
    //RANDOM SEED
    /*for(int i=0; i<n_seed; i++)
    {
        int r = rand()%nb_cells;
        printf("R %d\n", r);
        typename LCC_3::Dart_range::iterator it = lcc.one_dart_per_cell<2>().begin();
        for(int j=0; j<r; j++)
            it++;
        it->attribute<2>()->info().r = 0.2;
        it->attribute<2>()->info().g = 0.4;
        it->attribute<2>()->info().b = 0.6;
    }*/
    //Voir iterer sur les attributs
    typename LCC_3::Dart_range::iterator it = lcc.one_dart_per_cell<2>().begin();
    for(int j=0; j<nb_cells; j++)
    {
        it->attribute<2>()->info().r = 0.2;
        it->attribute<2>()->info().g = 0.4;
        it->attribute<2>()->info().b = 0.6;
        it++;
    }
    
    lcc.display_characteristics(std::cout) << ", valid=" << lcc.is_valid() << std::endl;
    display_lcc(lcc);

    return EXIT_SUCCESS;
}