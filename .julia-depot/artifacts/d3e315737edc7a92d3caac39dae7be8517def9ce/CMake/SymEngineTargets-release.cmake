#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "symengine" for configuration "Release"
set_property(TARGET symengine APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(symengine PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libsymengine.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libsymengine-0.11.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS symengine )
list(APPEND _IMPORT_CHECK_FILES_FOR_symengine "${_IMPORT_PREFIX}/lib/libsymengine.dll.a" "${_IMPORT_PREFIX}/bin/libsymengine-0.11.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
