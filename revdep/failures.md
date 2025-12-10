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
using C compiler: ‘Apple clang version 17.0.0 (clang-1700.4.4.1)’
sh: /opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran: No such file or directory
using SDK: ‘MacOSX26.1.sdk’
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2k.c -o 2k.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2ld.c -o 2ld.o
...
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c kin.morgan.c -o kin.morgan.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c makeped_c.c -o makeped_c.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c mia.c -o mia.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c muvar.c -o muvar.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c package_native_routine_registration_skeleton.c -o package_native_routine_registration_skeleton.o
/opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran  -fPIC  -Wall -g -O2  -c pfc.f -o pfc.o
make: /opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran: No such file or directory
make: *** [pfc.o] Error 1
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
using C compiler: ‘Apple clang version 17.0.0 (clang-1700.4.4.1)’
sh: /opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran: No such file or directory
using SDK: ‘MacOSX26.1.sdk’
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2k.c -o 2k.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c 2ld.c -o 2ld.o
...
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c kin.morgan.c -o kin.morgan.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c makeped_c.c -o makeped_c.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c mia.c -o mia.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c muvar.c -o muvar.o
clang -arch arm64 -std=gnu2x -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/opt/R/arm64/include    -fPIC  -falign-functions=64 -Wall -g -O2  -c package_native_routine_registration_skeleton.c -o package_native_routine_registration_skeleton.o
/opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran  -fPIC  -Wall -g -O2  -c pfc.f -o pfc.o
make: /opt/homebrew/Cellar/gcc/15.1.0/bin/gfortran: No such file or directory
make: *** [pfc.o] Error 1
ERROR: compilation failed for package ‘gap’
* removing ‘/Users/jayhesselberth/devel/rnabioco/valr/revdep/checks.noindex/gap/old/gap.Rcheck/gap’


```
