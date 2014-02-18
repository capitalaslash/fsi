Fluid Structure Interaction tests
=================================

the tests are built using the libMesh library (https://github.com/libMesh/libmesh).

Configuration
-------------

libMesh (git version):
````
$LIBMESH_SRC/configure --prefix=$LIBMESH_INSTALL \
    --disable-default-comm-world \
    --enable-parmesh \
    --with-methods="opt dbg" \
    --with-vtk-include=/usr/include/vtk-5.10 --with-vtk-lib=/usr/lib/vtk-5.10 \
    --enable-hdf5 \
    --disable-trilinos
````

FSItest:
````
cmake $FSItest_SRC \
    -DCMAKE_INSTALL_PREFIX:PATH=$FSItest_INSTALL
    -DCMAKE_BUILD_TYPE:STRING=Debug \
    -DLIBMESH_DIR:PATH=$LIBMESH_INSTALL
````
