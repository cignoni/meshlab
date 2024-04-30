# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES filter_ransac.cpp)

set(HEADERS filter_ransac.h ransac_plane_fitting.h)

add_meshlab_plugin(filter_ransac ${SOURCES} ${HEADERS})
