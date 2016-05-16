cmake_minimum_required(VERSION 3.0)

# Find the pulsar core
list(APPEND CMAKE_PREFIX_PATH "${PULSAR_PATH}")
find_package(pulsar REQUIRED)


if("${SIMINT_PATH}" STREQUAL "")
    message(FATAL_ERROR "SIMINT_PATH not set")
endif()

include_directories("${SIMINT_PATH}/include")
include_directories(${CMAKE_CURRENT_LIST_DIR})


# Common flags for standard modules
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-std=c++11")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-w3")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-restrict;-DRESTRICT=restrict")

    # some useless remarks
    #  1418  : external function definition with no prior declaration
    #  1419  : external declaration in primary source file (happens in boost files)
    #   383  : value copied to temporary, reference to temporary used (generally ok in c++)
    #   981  : operands are evaluated in unspecified order (generally ok in c++)
    #   869  : parameter "XXX" was never referenced.
    #          This is due to a bug in Intel before Intel 16 that causes this
    #          remark to be printed for variadic templates
    # 11074  : Inlining inhibited by limit max-size (etc)
    # 11076  : To get full report use -qopt-report=4 -qopt-report-phase ipo
    #   444  : destructor for base class .. not virtual (even warns for classes without virtual functions)
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd1418")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd1419")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd383")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd981")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd869")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd11074")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd11076")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-wd444")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-std=c++11")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-DRESTRICT=__restrict__")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-Wall;-Wextra;-pedantic")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-Wfloat-equal;-Wshadow")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-Wswitch-default")
    list(APPEND PULSAR_CXX_STRICT_FLAGS "-Wconversion")
endif()



set(SIMINT_FILES creator.cpp
                 SimintERI.cpp
   )

set (SIMINT_PYFILES
    )


add_library(SimintERI MODULE ${SIMINT_FILES})

target_link_libraries(SimintERI pulsar ${SIMINT_PATH}/lib/libsimint.a)
target_compile_options(SimintERI PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

install(TARGETS SimintERI LIBRARY DESTINATION SimintERI
                          RUNTIME DESTINATION SimintERI)


install(FILES "__init__.py" DESTINATION SimintERI)
install(FILES "modinfo.py"  DESTINATION SimintERI)
install(FILES ${SIMINT_PYFILES} DESTINATION SimintERI)
