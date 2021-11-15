lib_types = library('types', ['types.hpp'])
lib_odeint = library('odeint', ['odeint.hpp'])
lib_util = library('util', ['util.hpp'])
lib_collisions = library('collisions', ['collisions.hpp'])

lib_dep = declare_dependency(link_with: [lib_types, lib_odeint, lib_util, lib_collisions], include_directories: include_directories('.'))