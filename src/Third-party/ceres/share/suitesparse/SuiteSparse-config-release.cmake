#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "suitesparseconfig" for configuration "Release"
set_property(TARGET suitesparseconfig APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(suitesparseconfig PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/suitesparseconfig.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS suitesparseconfig )
list(APPEND _IMPORT_CHECK_FILES_FOR_suitesparseconfig "${_IMPORT_PREFIX}/lib/suitesparseconfig.lib" )

# Import target "amd" for configuration "Release"
set_property(TARGET amd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(amd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS amd )
list(APPEND _IMPORT_CHECK_FILES_FOR_amd "${_IMPORT_PREFIX}/lib/libamd.lib" )

# Import target "btf" for configuration "Release"
set_property(TARGET btf APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(btf PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libbtf.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS btf )
list(APPEND _IMPORT_CHECK_FILES_FOR_btf "${_IMPORT_PREFIX}/lib/libbtf.lib" )

# Import target "camd" for configuration "Release"
set_property(TARGET camd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(camd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS camd )
list(APPEND _IMPORT_CHECK_FILES_FOR_camd "${_IMPORT_PREFIX}/lib/libcamd.lib" )

# Import target "ccolamd" for configuration "Release"
set_property(TARGET ccolamd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ccolamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libccolamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS ccolamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_ccolamd "${_IMPORT_PREFIX}/lib/libccolamd.lib" )

# Import target "colamd" for configuration "Release"
set_property(TARGET colamd APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(colamd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcolamd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS colamd )
list(APPEND _IMPORT_CHECK_FILES_FOR_colamd "${_IMPORT_PREFIX}/lib/libcolamd.lib" )

# Import target "cholmod" for configuration "Release"
set_property(TARGET cholmod APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cholmod PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "amd;camd;colamd;ccolamd;suitesparseconfig;${_IMPORT_PREFIX}/lib/openblas.lib;${_IMPORT_PREFIX}/lib/lapack.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcholmod.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cholmod )
list(APPEND _IMPORT_CHECK_FILES_FOR_cholmod "${_IMPORT_PREFIX}/lib/libcholmod.lib" )

# Import target "cxsparse" for configuration "Release"
set_property(TARGET cxsparse APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(cxsparse PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libcxsparse.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS cxsparse )
list(APPEND _IMPORT_CHECK_FILES_FOR_cxsparse "${_IMPORT_PREFIX}/lib/libcxsparse.lib" )

# Import target "klu" for configuration "Release"
set_property(TARGET klu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(klu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "amd;colamd;btf;suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libklu.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS klu )
list(APPEND _IMPORT_CHECK_FILES_FOR_klu "${_IMPORT_PREFIX}/lib/libklu.lib" )

# Import target "ldl" for configuration "Release"
set_property(TARGET ldl APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ldl PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "suitesparseconfig"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libldl.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS ldl )
list(APPEND _IMPORT_CHECK_FILES_FOR_ldl "${_IMPORT_PREFIX}/lib/libldl.lib" )

# Import target "umfpack" for configuration "Release"
set_property(TARGET umfpack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(umfpack PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "amd;camd;colamd;ccolamd;cholmod;suitesparseconfig;${_IMPORT_PREFIX}/lib/openblas.lib;${_IMPORT_PREFIX}/lib/lapack.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libumfpack.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS umfpack )
list(APPEND _IMPORT_CHECK_FILES_FOR_umfpack "${_IMPORT_PREFIX}/lib/libumfpack.lib" )

# Import target "spqr" for configuration "Release"
set_property(TARGET spqr APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(spqr PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE "camd;ccolamd;cholmod;suitesparseconfig;${_IMPORT_PREFIX}/lib/openblas.lib;${_IMPORT_PREFIX}/lib/lapack.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libspqr.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS spqr )
list(APPEND _IMPORT_CHECK_FILES_FOR_spqr "${_IMPORT_PREFIX}/lib/libspqr.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
