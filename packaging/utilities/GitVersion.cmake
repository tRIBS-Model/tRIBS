################################################################################
# GitVersion.cmake -- stamps git provenance into the generated tribsVersion.h
#
# Run in script mode (cmake -P) from a build-time custom target, so the recorded
# hash and dirty state reflect the tree at build time rather than at configure
# time. Resolving this at configure time was misleading: editing a .cpp and
# rebuilding does not re-run cmake, so a binary could report a clean commit hash
# while containing uncommitted changes.
#
# Script mode starts a fresh cmake process with none of the project's variables,
# so everything needed is passed in with -D:
#
#   SRC_DIR        repository root, used as the git working directory
#   TEMPLATE_FILE  path to tribsVersion.h.in, the template alongside this script
#   OUTPUT_FILE    path to the generated tribsVersion.h
#   TRIBS_VERSION  semantic version, from CMakeLists.txt
#   TRIBS_RELEASE  release label, from CMakeLists.txt
#
# configure_file() leaves OUTPUT_FILE untouched when the generated content is
# identical, so a normal rebuild costs one short cmake run and recompiles
# nothing. Only an actual change in version or git state triggers a rebuild of
# the single translation unit that includes the header.
################################################################################

foreach(required SRC_DIR TEMPLATE_FILE OUTPUT_FILE TRIBS_VERSION TRIBS_RELEASE)
    if(NOT DEFINED ${required})
        message(FATAL_ERROR "GitVersion.cmake: ${required} was not passed in with -D")
    endif()
endforeach()

find_program(GIT_EXECUTABLE git)

# EXISTS covers worktrees too, where .git is a file rather than a directory.
if(GIT_EXECUTABLE AND EXISTS "${SRC_DIR}/.git")
    execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY "${SRC_DIR}"
            OUTPUT_VARIABLE TRIBS_GIT_SHA
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)

    # Non-zero exit means tracked files differ from HEAD, staged or not.
    # Untracked files are deliberately not counted as dirty.
    execute_process(
            COMMAND ${GIT_EXECUTABLE} diff --quiet HEAD
            WORKING_DIRECTORY "${SRC_DIR}"
            RESULT_VARIABLE TRIBS_GIT_DIRTY
            ERROR_QUIET)

    if(TRIBS_GIT_SHA AND NOT TRIBS_GIT_DIRTY EQUAL 0)
        set(TRIBS_GIT_SHA "${TRIBS_GIT_SHA}-dirty")
    endif()
endif()

# Deliberately left empty for source archives and git-less environments.
if(NOT TRIBS_GIT_SHA)
    set(TRIBS_GIT_SHA "")
endif()

configure_file("${TEMPLATE_FILE}" "${OUTPUT_FILE}" @ONLY)

################################################################################
#                          End of GitVersion.cmake
################################################################################
