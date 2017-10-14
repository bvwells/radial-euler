########################################
# Set up how to compile the source files
########################################

# Define the executable name
SET(RADIALEULEREXE radialeuler.exe)

# Add the source files
SET(RADIALEULER_Sources ${SRC}/main.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${RADIALEULEREXE} ${RADIALEULER_Sources})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${RADIALEULEREXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

set_target_properties( ${RADIALEULEREXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
