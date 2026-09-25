####################################################################################################
# Copyright (c) 2016 - 2021
# Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
#
# Author(s)
#       Marwan Abdellah <marwan.abdellah@epfl.ch>
#
# For complete list of authors, please see AUTHORS.md
#
# This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
#
# This library is free software; you can redistribute it and/or modify it under the terms of the
# GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301 USA.
# You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
####################################################################################################

set(_glm_HEADER_SEARCH_DIRS
    "/usr/include"
    "/usr/local/include"
    "${CMAKE_SOURCE_DIR}/includes"
    "C:/Program Files (x86)/glm"
)

# Check environment variable
set(_glm_ENV_ROOT_DIR "$ENV{GLM_ROOT_DIR}")

if(NOT GLM_ROOT_DIR AND _glm_ENV_ROOT_DIR)
    set(GLM_ROOT_DIR "${_glm_ENV_ROOT_DIR}")
endif(NOT GLM_ROOT_DIR AND _glm_ENV_ROOT_DIR)

# Put user specified location at beginning of search
if(GLM_ROOT_DIR)
    set(_glm_HEADER_SEARCH_DIRS "${GLM_ROOT_DIR}"
        "${GLM_ROOT_DIR}/include"
        ${_glm_HEADER_SEARCH_DIRS})
endif(GLM_ROOT_DIR)

# Locate header
find_path(GLM_INCLUDE_DIR "glm/glm.hpp"
    PATHS ${_glm_HEADER_SEARCH_DIRS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLM DEFAULT_MSG GLM_INCLUDE_DIR)

if(GLM_FOUND)
    set(GLM_INCLUDE_DIRS "${GLM_INCLUDE_DIR}")
    message(STATUS "GLM: ${GLM_INCLUDE_DIR}")
    link_libraries(${GLM_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_GLM")
    set(ULTRALISER_USE_GLM TRUE)
else(GLM_FOUND)
    set(ULTRALISER_USE_GLM FALSE)
endif(GLM_FOUND)
