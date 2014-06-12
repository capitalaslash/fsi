#include "biquad.hpp"

using namespace libMesh;

void generate_biquad( Mesh & mesh, GetPot & param )
{
    uint const flag_f = param("flag_f", 11);
    uint const flag_s = param("flag_s", 12);
    uint const nx = param("nx", 16);
    uint const ny = param("ny", 8);
    Real const ox = param("ox", 0.0);
    Real const lx = param("lx", 2.0);
    Real const oy = param("oy", 0.0);
    Real const ly = param("ly", 1.0);
    Real const li = param("li", 1.0);

    Real const toll = 1.e-8;

    bool const vert = param("mesh_vert", false);

    MeshTools::Generation::build_square (mesh,
                                         nx, ny,
                                         ox, lx,
                                         oy, ly,
                                         QUAD9);

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    mesh.boundary_info->clear();

    for ( ; el != end_el; ++el)
    {
        Elem* elem = *el;

        if (vert)
        {
            subdomain_id_type region_flag = flag_f;
            if (elem->point(8)(1) > li)
            {
                region_flag = flag_s;
            }
            elem->subdomain_id() = region_flag;

            for (uint s=0; s<elem->n_sides(); s++)
            {
                // check if side in on boundary
                AutoPtr<Elem> side (elem->build_side(s));
                Real side_x = side->point(2)(0);
                Real side_y = side->point(2)(1);

                subdomain_id_type side_flag = 0;

                if (elem->neighbor(s) == NULL)
                {
                    if ( side_y > oy && side_y < li &&
                         std::fabs( side_x - lx ) < toll )
                        side_flag = 1;
                    else if ( side_y > li && side_y < ly &&
                              std::fabs( side_x - lx ) < toll )
                        side_flag = 2;
                    else if ( std::fabs( side_y - ly) < toll &&
                              side_x > ox && side_x < lx )
                        side_flag = 3;
                    else if ( side_y > li && side_y < ly &&
                              std::fabs( side_x - ox ) < toll )
                        side_flag = 4;
                    else if ( side_y > oy && side_y < li &&
                              std::fabs( side_x - ox ) < toll )
                        side_flag = 5;
                    else if ( std::fabs( side_y - oy) < toll &&
                              side_x > ox && side_x < lx )
                        side_flag = 6;
                    else
                        abort();
                }
                // check if facet is on interface
                else if( std::fabs(side_y - li) < toll)
                {
                    side_flag = 10;
                }

                if(side_flag != 0)
                    mesh.boundary_info->add_side(elem,s,side_flag);
            }
        }
        else
        {
            subdomain_id_type region_flag = flag_f;
            if (elem->point(8)(0) > li)
            {
                region_flag = flag_s;
            }
            elem->subdomain_id() = region_flag;

            for (uint s=0; s<elem->n_sides(); s++)
            {
                // check if side in on boundary
                AutoPtr<Elem> side (elem->build_side(s));
                Real side_x = side->point(2)(0);
                Real side_y = side->point(2)(1);

                subdomain_id_type side_flag = 0;

                if (elem->neighbor(s) == NULL)
                {
                    if ( side_x > ox && side_x < li &&
                         std::fabs( side_y - oy ) < toll )
                        side_flag = 1;
                    else if ( side_x > li && side_x < lx &&
                              std::fabs( side_y - oy ) < toll )
                        side_flag = 2;
                    else if ( std::fabs( side_x - lx) < toll &&
                              side_y > oy && side_y < ly )
                        side_flag = 3;
                    else if ( side_x > li && side_x < lx &&
                              std::fabs( side_y - ly ) < toll )
                        side_flag = 4;
                    else if ( side_x > ox && side_x < li &&
                              std::fabs( side_y - ly ) < toll )
                        side_flag = 5;
                    else if ( std::fabs( side_x - ox) < toll &&
                              side_y > oy && side_y < ly )
                        side_flag = 6;
                    else
                        abort();
                }
                // check if facet is on interface
                else if( std::fabs(side_x - li) < toll)
                {
                    side_flag = 10;
                }

                if(side_flag != 0)
                    mesh.boundary_info->add_side(elem,s,side_flag);
            }
        }
    }

    if(nx*ny < 20)
        mesh.boundary_info->print_info();
}

void generate_barrel( Mesh & mesh, GetPot & param )
{
    uint const flag_f = param("flag_f", 11);
    uint const flag_s = param("flag_s", 12);
    uint const nx = param("nx", 10);
    uint const ny = param("ny", 10);
    Point const o( param("ox", 0.0), param("oy", 0.0) );
    Point const l( param("lx", 2.0), param("ly", 2.0) );
    Point const i( param("ix", 1.0), param("iy", 1.0) );

    Real const toll = 1.e-8;

    MeshTools::Generation::build_square (mesh,
                                         nx, ny,
                                         o(0), l(0),
                                         o(1), l(1),
                                         QUAD9);

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    mesh.boundary_info->clear();

    for ( ; el != end_el; ++el)
    {
        Elem* elem = *el;

        subdomain_id_type region_flag = flag_s;
        if (elem->point(8)(0) < i(0) && elem->point(8)(1) > i(1))
        {
            region_flag = flag_f;
        }
        elem->subdomain_id() = region_flag;

        for (uint s=0; s<elem->n_sides(); s++)
        {
            // check if side in on boundary
            AutoPtr<Elem> side (elem->build_side(s));
            Point const& centroid = side->centroid();

            subdomain_id_type side_flag = 0;

            if (elem->neighbor(s) == NULL)
            {
                if ( centroid(0) > o(0) && centroid(0) < l(0) &&
                     std::fabs( centroid(1) - o(1) ) < toll )
                    side_flag = 1;
                else if ( std::fabs( centroid(0) - l(0) ) < toll &&
                          centroid(1) > o(1) && centroid(1) < l(1) )
                    side_flag = 2;
                else if ( centroid(0) > i(0) && centroid(0) < l(0) &&
                          std::fabs( centroid(1) - l(1)) < toll)
                    side_flag = 3;
                else if ( centroid(0) > o(0) && centroid(0) < i(0) &&
                          std::fabs( centroid(1) - l(1) ) < toll )
                    side_flag = 4;
                else if ( std::fabs( centroid(0) - o(0) ) < toll &&
                          centroid(1) > i(1) && centroid(1) < l(1) )
                    side_flag = 5;
                else if ( std::fabs( centroid(0) - o(0)) < toll &&
                          centroid(1) > o(1) && centroid(1) < i(1) )
                    side_flag = 6;
                else
                    abort();
            }
            // check if facet is on interface
            else if( (centroid(0) > o(0) && centroid(0) < i(0) && std::fabs(centroid(1) - i(1)) < toll) ||
                     (centroid(1) > i(1) && centroid(1) < l(1) && std::fabs(centroid(0) - i(0)) < toll) )
            {
                side_flag = 10;
            }

            if(side_flag != 0)
                mesh.boundary_info->add_side(elem,s,side_flag);
        }
    }

    if(nx*ny < 20)
        mesh.boundary_info->print_info();
}
