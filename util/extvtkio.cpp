#include <iomanip>

#include <libmesh/equation_systems.h>

#include "extvtkio.hpp"

ExtVTKIO::ExtVTKIO (MeshBase const & mesh, Parameters const & par):
    VTKIO (mesh),
    _dirname (par.get<std::string>("output_dir")),
    _basename (par.get<std::string>("basename")),
    _print_step(par.get<uint>("print_step"))
{
    system (("mkdir -p " + _dirname).c_str());
    if( _dirname.substr( _dirname.size()-1,1) != "/" )
    {
        _dirname.append("/");
    }

#ifdef FSI_HAS_LIBXML2
    write_pvd( par.get<Real>("t_in"), par.get<Real>("t_out"), par.get<Real>("dt") );
#endif
}

void ExtVTKIO::write_solution(const EquationSystems & es,
                              const std::set<std::string> *system_names)
{
    std::stringstream file_name;
    file_name << _dirname << _basename << '_';
    file_name << std::setw(6) << std::setfill('0') << es.parameters.get<uint>("timestep");
    file_name << ".pvtu";
    VTKIO::write_equation_systems(file_name.str(), es, system_names);
}

#ifdef FSI_HAS_LIBXML2
void ExtVTKIO::write_pvd( Real const t_in, Real const t_out, Real const dt )
{
    xmlTextWriterPtr writer = xmlNewTextWriterFilename((_dirname + _basename + ".pvd").c_str(), 0);

    if (writer == NULL)
    {
      std::cerr << "testXmlwriterFilename: Error creating the xml writer" << std::endl;
      exit(1);
    }
    int rc = xmlTextWriterSetIndent(writer, 2);

    rc = xmlTextWriterStartDocument(writer, NULL, "UTF-8", NULL);
    if (rc < 0)
    {
      std::cerr << "testXmlwriterFilename: Error at xmlTextWriterStartDocument" << std::endl;
      exit(rc);
    }

    rc = xmlTextWriterStartElement(writer, BAD_CAST "VTKFile");
    xmlTextWriterWriteAttribute(writer, BAD_CAST "type",       BAD_CAST "Collection");
    xmlTextWriterWriteAttribute(writer, BAD_CAST "version",    BAD_CAST "0.1");
    xmlTextWriterWriteAttribute(writer, BAD_CAST "byte_order", BAD_CAST "LittleEndian");
    rc = xmlTextWriterStartElement(writer, BAD_CAST "Collection");

    Real time = t_in;
    uint timestep = 0;
    while( time <= t_out )
    {
      if ((timestep)%_print_step == 0)
      {
          std::stringstream ss;
          ss << time;
          xmlTextWriterStartElement(writer, BAD_CAST "DataSet");
          xmlTextWriterWriteAttribute(writer, BAD_CAST "timestep", BAD_CAST ss.str().c_str());

          std::stringstream file_name;
          file_name << _basename << '_';
          file_name << std::setw(6) << std::setfill('0') << timestep;
          file_name << ".pvtu";

          xmlTextWriterWriteAttribute(writer, BAD_CAST "file", BAD_CAST file_name.str().c_str());
          xmlTextWriterEndElement(writer);
      }

      time += dt;
      timestep++;
    }

    rc = xmlTextWriterEndDocument(writer);

    xmlFreeTextWriter(writer);
}
#endif
