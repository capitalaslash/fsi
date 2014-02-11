#include "biquad.hpp"

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
                    else if ( side_y > 0. && side_y < li &&
                              std::fabs( side_x - ox ) < toll )
                        side_flag = 5;
                    else if ( std::fabs( side_y - oy) < toll &&
                              side_x > ox && side_x < lx )
                        side_flag = 6;
                    else
                        abort();
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

                if(side_flag != 0)
                    mesh.boundary_info->add_side(elem,s,side_flag);
            }
        }
    }

    //mesh.boundary_info->print_info();
}
