# gap (1.6)

* GitHub: <https://github.com/jinghuazhao/R>
* Email: <mailto:jinghuazhao@hotmail.com>
* GitHub mirror: <https://github.com/cran/gap>

Run `revdepcheck::revdep_details(, "gap")` for more info

## In both

*   checking whether package ‘gap’ can be installed ... ERROR
     ```
     Installation failed.
     See ‘/Users/jayhesselberth/devel/rnabioco/valr/revdep/checks.noindex/gap/new/gap.Rcheck/00install.out’ for details.
     ```

## Installation

### Devel

```
* installing *source* package ‘gap’ ...
** this is package ‘gap’ version ‘1.6’
** package ‘gap’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
using C compiler: ‘Apple clang version 17.0.0 (clang-1700.6.3.2)’
using Fortran compiler: ‘GNU Fortran (GCC) 12.2.0’
using SDK: ‘MacOSX26.2.sdk’
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2k.c -o 2k.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2ld.c -o 2ld.o
...
clang: warning: overriding deployment version from '16.0' to '26.0' [-Woverriding-deployment-version]
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c pgc_c.c -o pgc_c.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c whscore.c -o whscore.o
clang -arch arm64 -std=gnu2x -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o gap.so 2k.o 2ld.o cline.o gcontrol_c.o gcx.o gif_c.o hap_c.o hwe.hardy.o kin.morgan.o makeped_c.o mia.o muvar.o package_native_routine_registration_skeleton.o pfc.o pfc.sim.o pgc_c.o whscore.o -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lemutls_w -lheapt_w -lgfortran -lquadmath -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0' not found
ld: library 'emutls_w' not found
clang: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [gap.so] Error 1
ERROR: compilation failed for package ‘gap’
* removing ‘/Users/jayhesselberth/devel/rnabioco/valr/revdep/checks.noindex/gap/new/gap.Rcheck/gap’


```
### CRAN

```
* installing *source* package ‘gap’ ...
** this is package ‘gap’ version ‘1.6’
** package ‘gap’ successfully unpacked and MD5 sums checked
** using staged installation
** libs
using C compiler: ‘Apple clang version 17.0.0 (clang-1700.6.3.2)’
using Fortran compiler: ‘GNU Fortran (GCC) 12.2.0’
using SDK: ‘MacOSX26.2.sdk’
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2k.c -o 2k.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2ld.c -o 2ld.o
...
clang: warning: overriding deployment version from '16.0' to '26.0' [-Woverriding-deployment-version]
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c pgc_c.c -o pgc_c.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c whscore.c -o whscore.o
clang -arch arm64 -std=gnu2x -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -L/Library/Frameworks/R.framework/Resources/lib -L/opt/R/arm64/lib -o gap.so 2k.o 2ld.o cline.o gcontrol_c.o gcx.o gif_c.o hap_c.o hwe.hardy.o kin.morgan.o makeped_c.o mia.o muvar.o package_native_routine_registration_skeleton.o pfc.o pfc.sim.o pgc_c.o whscore.o -L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0 -L/opt/gfortran/lib -lemutls_w -lheapt_w -lgfortran -lquadmath -F/Library/Frameworks/R.framework/.. -framework R
ld: warning: search path '/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/14.2.0' not found
ld: library 'emutls_w' not found
clang: error: linker command failed with exit code 1 (use -v to see invocation)
make: *** [gap.so] Error 1
ERROR: compilation failed for package ‘gap’
* removing ‘/Users/jayhesselberth/devel/rnabioco/valr/revdep/checks.noindex/gap/old/gap.Rcheck/gap’


```
