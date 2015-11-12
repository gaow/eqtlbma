/** \file grid.hpp
 *
 *  `Grid' is a class
 *  Copyright (C) 2013 Timothee Flutre
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef QUANTGEN_GRID_HPP
#define QUANTGEN_GRID_HPP

#include <cstdlib>

#include <vector>
#include <string>
#include <iostream>
#include <gsl/gsl_matrix.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"

namespace quantgen {

  class Grid {
  public:
    std::vector<double> phi2s;
    std::vector<double> oma2s;
    std::vector<double> phi2s_fix;
    std::vector<double> oma2s_fix;
    std::vector<double> phi2s_maxh;
    std::vector<double> oma2s_maxh;
    Grid();
    Grid(const std::string & gridFile, const bool & makeFixMaxh, const int & verbose);
    size_t size (void) const { return phi2s.size(); }
  };

  class PriorMatrices {
  public:
    std::vector<gsl_matrix*> Wgs;
    std::vector<double> Wg_scalars;
    std::vector<std::string> Wg_names;
    PriorMatrices();
    PriorMatrices(const std::string & file_pattern,
                  const std::string & scalar_file,
                  const size_t & matrix_dimension,
                  const int & verbose);
    size_t size (void) const { return Wgs.size(); }
    ~PriorMatrices() {
      for (size_t m = 0; m < Wgs.size(); ++m)
        if (Wgs[m])
          gsl_matrix_free(Wgs[m]);
    }
  };

} // namespace quantgen

#endif // QUANTGEN_GRID_HPP
