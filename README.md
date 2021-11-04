This repository uses meson [meson](https://mesonbuild.com/) and [ninja](https://ninja-build.org/) as a build system. In order to build the project:

```bash
meson build
cd build
ninja
```

### Stucture

#### `types.hpp` 

Contains and implements types used to store vectors and states (position, velocity) and defines the operations with them

#### `odeint.hpp`

Defines (and implements) the motion equation for the particles and the methods for integrating them.

#### `util.hpp`

Defines (and implements) methods that are mos specific (like, for instance, the ASDEX magnetic field fucntion `B_ASDEX`)

