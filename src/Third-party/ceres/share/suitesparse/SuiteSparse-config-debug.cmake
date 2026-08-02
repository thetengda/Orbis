#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "suitesparseconfig" for configuration "Debug"
set_property(TARGET suitesparseconfig APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(suitesparseconfig PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/suitesparseconfig.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS suitesparseconfig )
list(APPEND _IMPORT_CHECK_FILES_FOR_suitesparseconfig "${_IMPORT_PREFIX}/debug/lib/suitesparseconfig.lib" )

# Import target "amd" for configuration "Debug"
set_property(TARGET amd APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(amd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS amd )
list(APPEND _IMPORT_CHECK_FILES_FOR_amd "${_IMPORT_PREFIX}/debug/lib/libamd.lib" )

# Import target "btf" for configuration "Debug"
set_property(TARGET btf APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(btf PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libbtf.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS btf )
list(APPEND _IMPORT_CHECK_FILES_FOR_btf "${_IMPORT_PREFIX}/debug/lib/libbtf.lib" )

# Import target "camd" for configuration "Debug"
set_property(TARGET camd APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(camd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libcamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS camd )
list(APPEND _IMPORT_CHECK_FILES_FOR_camd "${_IMPORT_PREFIX}/debug/lib/libcamd.lib" )

# Import target "ccolamd" for configuration "Debug"
set_property(TARGET ccolamd APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ccolamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libccolamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS ccolamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_ccolamd "${_IMPORT_PREFIX}/debug/lib/libccolamd.lib" )

# Import target "colamd" for configuration "Debug"
set_property(TARGET colamd APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(colamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libcolamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS colamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_colamd "${_IMPORT_PREFIX}/debug/lib/libcolamd.lib" )

# Import target "cholmod" for configuration "Debug"
set_property(TARGET cholmod APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(cholmod PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "amd;camd;colamd;ccolamd;suitesparseconfig;${_IMPORT_PREFIX}/debug/lib/openblas.lib;${_IMPORT_PREFIX}/debug/lib/lapack.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libcholmod.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cholmod )
list(APPEND _IMPORT_CHECK_FILES_FOR_cholmod "${_IMPORT_PREFIX}/debug/lib/libcholmod.lib" )

# Import target "cxsparse" for configuration "Debug"
set_property(TARGET cxsparse APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(cxsparse PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libcxsparse.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cxsparse )
list(APPEND _IMPORT_CHECK_FILES_FOR_cxsparse "${_IMPORT_PREFIX}/debug/lib/libcxsparse.lib" )

# Import target "klu" for configuration "Debug"
set_property(TARGET klu APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(klu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "amd;colamd;btf;suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libklu.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS klu )
list(APPEND _IMPORT_CHECK_FILES_FOR_klu "${_IMPORT_PREFIX}/debug/lib/libklu.lib" )

# Import target "ldl" for configuration "Debug"
set_property(TARGET ldl APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(ldl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "suitesparseconfig"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libldl.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS ldl )
list(APPEND _IMPORT_CHECK_FILES_FOR_ldl "${_IMPORT_PREFIX}/debug/lib/libldl.lib" )

# Import target "umfpack" for configuration "Debug"
set_property(TARGET umfpack APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(umfpack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "amd;camd;colamd;ccolamd;cholmod;suitesparseconfig;${_IMPORT_PREFIX}/debug/lib/openblas.lib;${_IMPORT_PREFIX}/debug/lib/lapack.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libumfpack.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS umfpack )
list(APPEND _IMPORT_CHECK_FILES_FOR_umfpack "${_IMPORT_PREFIX}/debug/lib/libumfpack.lib" )

# Import target "spqr" for configuration "Debug"
set_property(TARGET spqr APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(spqr PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_DEBUG "camd;ccolamd;cholmod;suitesparseconfig;${_IMPORT_PREFIX}/debug/lib/openblas.lib;${_IMPORT_PREFIX}/debug/lib/lapack.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/debug/lib/libspqr.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS spqr )
list(APPEND _IMPORT_CHECK_FILES_FOR_spqr "${_IMPORT_PREFIX}/debug/lib/libspqr.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
