project: roots-fortran
output_dir: ./doc
media_dir: ./media
src_dir: ./src
project_github: https://github.com/jacobwilliams/roots-fortran
summary: Root solvers for modern Fortran
author: Jacob Williams
github: https://github.com/jacobwilliams
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
source: true
graph: true
search: true
preprocessor: gfortran -E
exclude: **/pyplot_module.f90
         **/face.F90
exclude_dir: ./test
extra_mods: pyplot_module:https://github.com/jacobwilliams/pyplot-fortran
            iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html

{!README.md!}