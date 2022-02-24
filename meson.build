project('TetMesh', ['cpp'], version : '1.2.1')

smallvector = subproject('SmallVector').get_variable('smallvector_dep')

tetmesh_dep = declare_dependency(include_directories : include_directories('.'),
  dependencies : smallvector)

