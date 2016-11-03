cmake_minimum_required (VERSION 2.8.7)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)

project (NeuroBayes)
enable_language (Fortran)

if(DEFINED ENV{MAXNODE_PREPRO})
    set(MAXNODE_PREPRO $ENV{MAXNODE_PREPRO})
else()
    set(MAXNODE_PREPRO 100)
endif()

set(CMAKE_Fortran_FLAGS " -fno-automatic -fcheck=bound -Wno-unused -cpp ")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++0x -Wall -Wextra")

configure_file (
  "${PROJECT_SOURCE_DIR}/cpp_interface/nb_param.hh.in"
  "${PROJECT_SOURCE_DIR}/cpp_interface/nb_param.hh"
)

add_library(NeuroBayesCore_shared SHARED IMPORTED)
set_target_properties(NeuroBayesCore_shared PROPERTIES IMPORTED_LOCATION "${PROJECT_SOURCE_DIR}/core/libNeuroBayesCore_shared.so")

set(NeuroBayesCPPInterface_HEADERS
    "${PROJECT_SOURCE_DIR}/cpp_interface/nb_cpp_utils.h"
    "${PROJECT_SOURCE_DIR}/cpp_interface/NeuroBayesExpert.hh"
    "${PROJECT_SOURCE_DIR}/cpp_interface/NeuroBayesTeacher.hh"
    "${PROJECT_SOURCE_DIR}/cpp_interface/nb_param.hh"
    "${PROJECT_SOURCE_DIR}/cpp_interface/nb_c_utils.h"
)

set(NeuroBayesCPPInterface_SOURCES
    "${PROJECT_SOURCE_DIR}/cpp_interface/NeuroBayesExpert.cc"
    "${PROJECT_SOURCE_DIR}/cpp_interface/NeuroBayesTeacher.cc"
    "${PROJECT_SOURCE_DIR}/cpp_interface/nb_cpp_utils.cc"
)

add_library(NeuroBayesCPPInterface_shared SHARED ${NeuroBayesCPPInterface_SOURCES} ${NeuroBayesCPPInterface_HEADERS})
target_link_libraries(NeuroBayesCPPInterface_shared NeuroBayesCore_shared)

include_directories("${PROJECT_SOURCE_DIR}/cpp_interface/")

add_executable(crosstest "${PROJECT_SOURCE_DIR}/cpp_interface_tests/crosstest.cc")
target_link_libraries(crosstest NeuroBayesCPPInterface_shared)

add_executable(interfacetest "${PROJECT_SOURCE_DIR}/cpp_interface_tests/interfacetest.cc" ${NeuroBayesCPPInterface_HEADERS})
target_link_libraries(interfacetest NeuroBayesCPPInterface_shared)

add_executable(minimaltest "${PROJECT_SOURCE_DIR}/cpp_interface_tests/minimaltest.cc" ${NeuroBayesCPPInterface_HEADERS})
target_link_libraries(minimaltest NeuroBayesCPPInterface_shared)


install(TARGETS NeuroBayesCPPInterface_shared minimaltest
     LIBRARY DESTINATION lib
     ARCHIVE DESTINATION lib
     RUNTIME DESTINATION bin   
)

get_target_property(NeuroBayesCore_LOCATION NeuroBayesCore_shared LOCATION)
install(FILES ${NeuroBayesCore_LOCATION} DESTINATION lib)
install(FILES ${NeuroBayesCPPInterface_HEADERS} DESTINATION include)

set(CPACK_PACKAGE_VERSION 4.3)
set(CPACK_GENERATOR "RPM;DEB;TGZ")
set(CPACK_PACKAGE_NAME "NeuroBayes")
set(CPACK_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "neurobayes@blue-yonder.com")
set(CPACK_PACKAGE_VENDOR "Blue Yonder")
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/neurobayes.spec.in" "${CMAKE_CURRENT_BINARY_DIR}/neurobayes.spec" @ONLY IMMEDIATE)
set(CPACK_RPM_USER_BINARY_SPECFILE "${CMAKE_CURRENT_BINARY_DIR}/neurobayes.spec")

SET(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
SET(CPACK_DEBIAN_PACKAGE_SECTION "libs")
SET(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})


include(CPack)
