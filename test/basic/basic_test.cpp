#include <iostream>

#include <libmesh/libmesh.h>
#include <libmesh/mesh.h>

using namespace libMesh;

int main(int argc, char *argv[])
{
    LibMeshInit init(argc,argv);

    std::cout << "hello world!" << std::endl;
    return 0;
}
