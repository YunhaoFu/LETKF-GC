cmake_minimum_required (VERSION 3.2)
enable_language(Fortran)

project (LETKF-GRAPES_CoLM)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "Build type, options are: Debug Release RelWithDebInfo" FORCE)
endif()

# set(CMAKE_Fortran_COMPILER,"/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/ifort")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# MPI
#find_package(MPI REQUIRED)

# NetCDF
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

# LAPACK / BLAS
find_package(LAPACK REQUIRED)

set( da_lib_srcs
    precision.F90
    ini_fnl.F90
    io_control.F90
    prepost.F90
    kdtree.F90
    localization.F90
    letkf.F90
)

add_library( da STATIC ${da_lib_srcs} )

# LAPACK
target_link_libraries(da  PRIVATE ${LAPACK_LIBRARIES})
#NetCDF
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries(da  PRIVATE  ${NETCDF_LIBRARIES})

add_executable(letkf_run.exe letkf_run.F90)

set_target_properties(letkf_run.exe PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(letkf_run.exe da)
