project('TetMesh', ['cpp'], version : '1.2.1')

smallvector = subproject('SmallVector').get_variable('smallvector')

tetmesh = declare_dependency(include_directories : include_directories('.'),
  dependencies : smallvector)

