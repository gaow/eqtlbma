/** \file grid.cpp
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

#include "quantgen/grid.hpp"

using namespace std;

using namespace utils;

namespace quantgen {

  Grid::Grid(const string & gridFile, const bool & makeFixMaxh,
	     const int & verbose)
  {
    if (! gridFile.empty()) {
      if (verbose > 0)
	cout << "load grid in " << gridFile << " ..." << endl << flush;

      gzFile gridStream;
      vector<string> tokens;
      string line;
      openFile(gridFile, gridStream, "rb");
      while (getline(gridStream, line)) {
	split(line, " \t", tokens);
	if (tokens.size() != 2) {
	  cerr << "ERROR: format of file " << gridFile
	       << " should be phi2<space/tab>oma2" << endl;
	  exit(EXIT_FAILURE);
	}
	phi2s.push_back(atof(tokens[0].c_str()));
	oma2s.push_back(atof(tokens[1].c_str()));
	if (makeFixMaxh) {
	  phi2s_fix.push_back(0.0);
	  oma2s_fix.push_back(atof(tokens[0].c_str())
			      + atof(tokens[1].c_str()));
	  phi2s_maxh.push_back(atof(tokens[0].c_str())
			       + atof(tokens[1].c_str()));
	  oma2s_maxh.push_back(0.0);
	}
      }
      if (! gzeof(gridStream)) {
	cerr << "ERROR: can't read successfully file " << gridFile
	     << " up to the end" << endl;
	exit(EXIT_FAILURE);
      }
      closeFile(gridFile, gridStream);

      if (verbose > 0)
	cout << "grid size: " << phi2s.size() << endl;
    }
  }

  PriorMatrices::PriorMatrices(const std::string & file_pattern,
                               const std::string & scalar_file,
                               const size_t & matrix_dimension,
                               const int & verbose) {
    if (file_pattern.empty() || scalar_file.empty()) {
      return;
    }
    // load prior matrix by file patterns
    Wg_names = glob(file_pattern);
    if(Wg_names.size() == 0){
      cerr << "ERROR: no input prior matrix was found from pattern " <<
        file_pattern << endl;
      exit(EXIT_FAILURE);
    }
    for (size_t m = 0; m < Wg_names.size(); ++m) {
      gsl_matrix * priorM = gsl_matrix_alloc(matrix_dimension,
                                             matrix_dimension);
      gzFile gridStream;
      vector<string> tokens;
      string line;
      size_t nb_line = 0;
      openFile(Wg_names[m], gridStream, "rb");
      while (getline(gridStream, line)) {
        split(line, " \t", tokens);
        if (tokens.size() != matrix_dimension) {
        cerr << "ERROR: file " << Wg_names[m]
             << " should have " << matrix_dimension << "columns"
             << endl;
        gsl_matrix_free(priorM);
        exit(EXIT_FAILURE);
        }
        for (size_t i = 0; i < tokens.size(); ++i) {
          gsl_matrix_set(priorM, nb_line, i, atof(tokens[0].c_str()));
        }
        nb_line += 1;
      }
      if (nb_line > matrix_dimension) {
        cerr << "ERROR: file " << Wg_names[m]
             << " should have " << matrix_dimension << "rows"
             << endl;
        gsl_matrix_free(priorM);
        exit(EXIT_FAILURE);
      }
      if (! gzeof(gridStream)) {
        cerr << "ERROR: can't read successfully file " << Wg_names[m]
             << " up to the end" << endl;
        gsl_matrix_free(priorM);
        exit(EXIT_FAILURE);
      }
      closeFile(Wg_names[m], gridStream);
      Wgs.push_back(priorM);
    }
    if (verbose > 0) {
        cout << "number of customized prior matrices: "
             << Wgs.size() << endl;
    }

    // load prior matrix scalars
    gzFile gridStream;
    vector<string> tokens;
    string line;
    openFile(scalar_file, gridStream, "rb");
    while (getline(gridStream, line)) {
      split(line, " \t", tokens);
      if (tokens.size() != 1) {
        cerr << "ERROR: file " << scalar_file
             << " should have a single column of values and no empty rows"
             << endl;
        exit(EXIT_FAILURE);
      }
      Wg_scalars.push_back(atof(tokens[0].c_str()));
    }
    if (! gzeof(gridStream)) {
      cerr << "ERROR: can't read successfully file " << scalar_file
           << " up to the end" << endl;
      exit(EXIT_FAILURE);
    }
    closeFile(scalar_file, gridStream);

    if (verbose > 0) {
        cout << "customized prior scalar grid size: " << Wg_scalars.size() << endl;
    }
  }

} //namespace quantgen
