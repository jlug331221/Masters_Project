cmake_minimum_required(VERSION 3.9)
project(Masters_Project C)

set(CMAKE_C_STANDARD 99)

add_executable(Masters_Project main.c Source_Files/LSH.c Header_Files/LSH.h Header_Files/Headers.h Header_Files/Defs.h Header_Files/kdtree.h Source_Files/kdtree.c Header_Files/bkmeans_z.h Source_Files/bkmeans_z.c Header_Files/helping_procedures.h Source_Files/helping_procedures.c Header_Files/kdtree_median.h Source_Files/kdtree_median.c Header_Files/bkmeans_j.h Source_Files/bkmeans_j.c Source_Files/LSH_cluster_ADT.c Header_Files/LSH_cluster_ADT.h)