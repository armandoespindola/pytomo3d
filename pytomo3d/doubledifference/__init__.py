#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Double difference module

:copyright:
   Ridvan Orsvuran (orsvuran@geoazur.unice.fr), 2017
:license:
    GNU Lesser General Public License, version 3 (LGPLv3)
    (http://www.gnu.org/licenses/lgpl-3.0.en.html)
"""

from .adjoint import (calculate_adjoint_pair,  # NOQA
                      add_adjoint_sources,
                      multiply_adjoint_source)

from .pairing import find_pairs  # NOQA

from .utils import (deconstruct_winname,  # NOQA
                    get_stanames_of_pair)

from .windows import (component_based_windows_data,  # NOQA
                      convert_to_sta_based_windows,
                      filter_paired_windows)
