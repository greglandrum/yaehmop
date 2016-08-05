cmake_minimum_required(VERSION 2.6)

project(yaehmop C)

# Source files
set(YAEHMOP_SRCS
  avg_props.c
  bands.c
  charge_mat.c
  chg_it.c
  COOP_stuff.c
  distance_mat.c
  DOS_stuff.c
  electrostat.c
  fileio.c
  FMO_stuff.c
  genutil.c
  geom_frags.c
  globals.c
  K_hamil.c
  K_overlap_mat.c
  kpoints.c
  main.c
  matrices.c
  memory.c
  mod_mulliken.c
  mov.c
  muller.c
  mulliken.c
  netCDF_support.c
  new3_fileio.c
  postprocess.c
  princ_axes.c
  R_hamil.c
  R_overlap_mat.c
  recip_space.c
  solid_symmetry.c
  symmetry.c
  transforms.c
  walsh.c
  xtal_coords.c
  zetas.c
  Zmat.c
)

# Should we use F2C?
option(USE_F2C
  "Whether to compile using only C compilers with the F2C library
   Note that you must be able to link to the f2c library for this to work."
  ON)

if(USE_F2C)
  # We link to the static library of F2C to avoid an error
  set(F2C_LIB "f2c.a")

  message("-- Building with F2C. Fortran compilers are NOT required.")
  message("-- The F2C library IS required, however.")
  message("-- F2C_LIB is currently: ${F2C_LIB}")
  message("-- Set the F2C library with -DF2C_LIB=<f2c_loc>. Static library is preferred.")
  message("-- Turn F2C off with -DUSE_F2C=OFF")

  set(YAEHMOP_SRCS
    ${YAEHMOP_SRCS}
    f2c_files/abfns.c
    f2c_files/cboris.c
    f2c_files/diag.c
    f2c_files/lovlap.c
  )
else(USE_F2C)
  message("-- Building without F2C. Fortran compilers ARE required.")
  # We have to find a Fortran compiler if we are not using F2C
  enable_language(Fortran)
  set(YAEHMOP_SRCS
    ${YAEHMOP_SRCS}
    abfns.f
    cboris.f
    diag.f
    lovlap.f
  )
endif(USE_F2C)

# We link to the static library to avoid an error

# The location of the parameters file needs to be defined
set(PARM_FILE_LOC ${yaehmop_SOURCE_DIR}/eht_parms.dat)
message("-- Parameter file is currently set to ${PARM_FILE_LOC}")
message("-- Change the parameter file with -DPARM_FILE_LOC=<location>")
add_definitions(-DEHT_PARM_FILE=\"${PARM_FILE_LOC}\")

# This just adds an underscore after some function names in the source...
add_definitions(-DUNDERSCORE_FORTRAN)

# Add the executable and link to the C math library
add_executable(bind ${YAEHMOP_SRCS})
target_link_libraries(bind m)

# If we are using F2C, we have to link to the library
if(USE_F2C)
  target_link_libraries(bind ${F2C_LIB})
endif(USE_F2C)
