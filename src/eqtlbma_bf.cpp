/** \file eqtlbma_bf.cpp
 *
 *  `eqtlbma_bf' performs eQTL mapping in multiple subgroups via a Bayesian model.
 *  Copyright (C) 2012-2015 Timothée Flutre
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

#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <libgen.h>
#include <sys/stat.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <sstream>
#include <numeric>
using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>

#include <omp.h>

#include "utils/utils_io.hpp"
#include "utils/utils_math.hpp"
using namespace utils;

#include "quantgen/data_loader.hpp"
#include "quantgen/gene.hpp"
#include "quantgen/snp.hpp"
#include "quantgen/gene_snp_pair.hpp"
using namespace quantgen;

#include "tabix/bgzf.h"
#include "tabix/tabix.h"

/** \brief Display the help on stdout.
 */
void help(char ** argv)
{
  cout << "`" << argv[0] << "'"
       << " performs eQTL mapping in multiple subgroups via a Bayesian model." << endl
       << endl
       << "Usage: " << argv[0] << " [OPTIONS] ..." << endl
       << endl
       << "Options:" << endl
       << "  -h, --help\tdisplay the help and exit" << endl
       << "  -V, --version\toutput version information and exit" << endl
       << "  -v, --verbose\tverbosity level (0/default=1/2/3)" << endl
       << "      --geno\tfile with absolute paths to genotype files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: can be in three formats (VCF/IMPUTE/custom)" << endl
       << "\t\tVCF: see specifications on 1kG website" << endl
       << "\t\tIMPUTE: row 1 is a header chr<del>name<del>coord<del>a1<del>a2" << endl
       << "\t\t followed by >sample1_a1a1<del>sample1_a1a2<del>sample1_a2a2<del>..." << endl
       << "\t\tcustom: genotypes as allele dose, same as for MatrixEQTL" << endl
       << "\t\t and missing data can be NA or -1 (as used by vcftools --012)" << endl
       << "      --scoord\tfile with the SNP coordinates" << endl
       << "\t\tcompulsory if custom genotype format; forbidden otherwise" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tSNPs in the genotype files without coordinate are skipped (see also --snp)" << endl
       << "\t\tif a tabix-indexed file is also present, it will be used" << endl
       << "      --exp\tfile with absolute paths to expression level files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: custom format, same as for MatrixEQTL" << endl
       << "\t\t row 1 for sample names, column 1 for gene names" << endl
       << "\t\tsubgroups can have different genes" << endl
       << "\t\tall genes should be in the --gcoord file" << endl
       << "      --gcoord\tfile with the gene coordinates" << endl
       << "\t\tshould be in the BED format (delimiter: tab)" << endl
       << "\t\tgenes in the exp level files without coordinates are skipped" << endl
       << "      --anchor\tgene boundary(ies) for the cis region" << endl
       << "\t\tdefault=TSS (assumed to be start in BED file)" << endl
       // << "\t\tdefault=TSS, can also be TSS+TES" << endl
       << "      --cis\tlength of half of the cis region (radius, in bp)" << endl
       << "\t\tapart from the anchor(s), default=100000" << endl
       << "      --inss\tfile with absolute paths to files with summary statistics" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsstats file: custom format, similar to the one from --outss (see below)" << endl
       << "\t\t header should have gene, snp, n, sigmahat, betahat.geno and sebetahat.geno" << endl
       << "\t\t order doesn't matter" << endl
       << "      --out\tprefix for the output files" << endl
       << "\t\tall output files are gzipped and have a header line" << endl
       << "      --lik\tlikelihood to use" << endl
       << "\t\t'normal' (default)" << endl
       << "\t\t'poisson' or 'quasipoisson'" << endl
       << "      --analys\tanalysis to perform" << endl
       << "\t\t'sep': separate analysis of each subgroup" << endl
       << "\t\t'join': joint analysis of all subgroups" << endl
       << "      --outss\twrite the output file with all summary statistics" << endl
       << "      --outw\twrite the output file with the ABFs averaged over the grid" << endl
       << "\t\tgrid weights are uniformly equal" << endl
       << "      --qnorm\tquantile-normalize the exp levels to a N(0,1)" << endl
       << "      --maf\tminimum minor allele frequency (default=0.0)" << endl
       << "      --covar\tfile with absolute paths to covariate files" << endl
       << "\t\ttwo columns: subgroup identifier<space/tab>path to file" << endl
       << "\t\tcan be a single line (single subgroup)" << endl
       << "\t\tadd '#' at the beginning of a line to comment it" << endl
       << "\t\tsubgroup file: row 1 is a header sample<space/tab>covariate1 ..." << endl
       << "\t\tall sample names should be in the respective genotype and exp level files" << endl
       << "\t\tthe covariates should be numbers, no missing value is allowed" << endl
       << "\t\tsubgroups can have different covariates" << endl
       << "\t\tthe order of rows is not important" << endl
       << "      --gridL\tfile with a 'large' grid for prior variances in standardized effect sizes" << endl
       << "\t\tfirst column is phi^2 and second column is omega^2, no header" << endl
       << "\t\tthis grid is used with model 1 ('general alternative') trying to capture" << endl
       << "\t\t all sorts of heterogeneity" << endl
       << "\t\trequired with --analys join" << endl
       << "      --gridS\tfile with a 'small' grid of values for phi^2 and omega^2" << endl
       << "\t\tsame format as --gridL" << endl
       << "\t\trequired with --analyis join if --bfs is 'sin' or 'all' or 'customized'" << endl
       << "      --priorM\tfiles of matrices to be used as prior values" << endl
       << "\t\teach file contains one RxR matrix where R is number of subgroups" << endl
       << "\t\trequired with --analyis join if --bfs is 'customized'" << endl
       << "\t\twhen used, --error must be 'mvlr' or 'hybrid'" << endl
       << "\t\twildcard character (*) for file names is supported (file pattern string must be double quoted!)" << endl
       << "      --gridM\tfile with grid of scalars for customized prior matrices" << endl
       << "\t\tthis is a one column file containing all scalar values" << endl
       << "\t\tthis grid should be 'smaller' than gridL" << endl
       << "\t\trequired when --priorM is specified" << endl
       << "      --bfs\twhich Bayes Factors to compute for the joint analysis" << endl
       << "\t\tonly the Laplace-approximated BF from Wen and Stephens (AoAS 2013) is implemented" << endl
       << "\t\tif --outw, each BF for a given configuration is the average of the BFs over one of the grids, with equal weights" << endl
       << "\t\t'gen' (default): general way to capture any level of heterogeneity" << endl
       << "\t\t correspond to the consistent configuration with the large grid" << endl
       << "\t\t fixed-effect and maximum-heterogeneity BFs are also calculated" << endl
       << "\t\t'sin': compute also the BF for each singleton (subgroup-specific configuration)" << endl
       << "\t\t they use the small grid (BF_BMAlite is also reported)" << endl
       << "\t\t'customized': compute also the BF for customized prior specifications" << endl
       << "\t\t customized BFs use the grid provided by --gridM (BF_BMACustomized is also reported)" << endl
       << "\t\t'all': compute the BFs for all configurations (costly if many subgroups)" << endl
       << "\t\t all BFs use the small grid (BF_BMA is also reported)" << endl
       << "      --error\tmodel for the errors (if --analys join)" << endl
       << "\t\t'uvlr': default, errors are not correlated between subgroups (different individuals)" << endl
       << "\t\t'mvlr': errors can be correlated between subgroups (same individuals)" << endl
       << "\t\t'hybrid': errors can be correlated between pairs of subgroups (common individuals)" << endl
       << "      --fiterr\tparam used when estimating the variance of the errors (if --analys join, only with 'mvlr' or 'hybrid')" << endl
       << "\t\tdefault=0.5 but can be between 0 (null model) and 1 (full model)" << endl
       << "      --nperm\tnumber of permutations" << endl
       << "\t\tdefault=0, otherwise 10000 is recommended" << endl
       << "      --seed\tseed for the two random number generators" << endl
       << "\t\tone for the permutations, another for the trick" << endl
       << "\t\tby default, both are initialized via microseconds from epoch" << endl
       << "\t\tthe RNGs are re-seeded before each subgroup and before the joint analysis" << endl
       << "\t\tthis, along with --trick 2, allows for proper comparison of separate and joint analyzes" << endl
       << "      --trick\tapply trick to speed-up permutations" << endl
       << "\t\tstop after the tenth permutation for which the test statistic" << endl
       << "\t\t is better than or equal to the true value, and sample from" << endl
       << "\t\t a uniform between 11/(nbPermsSoFar+2) and 11/(nbPermsSoFar+1)" << endl
       << "\t\tif '1', the permutations really stops" << endl
       << "\t\tif '2', all permutations are done but the test statistics are not computed" << endl
       << "\t\tallows to compare different test statistics on the same permutations" << endl
       << "      --tricut\tcutoff for the trick (default=10)" << endl
       << "\t\tstop permutations once the nb of permutations for which permTestStat is more extreme" << endl
       << "\t\t than trueTestStat equals this cutoff" << endl
       << "      --permsep\twhich permutation procedure for the separate analysis" << endl
       << "\t\t0 (default): no permutations are done for the separate analysis" << endl
       << "\t\t1: use the minimum P-value over SNPs and subgroups as a test statistic (keeps correlations)" << endl
       << "\t\t2: use the minimum P-value over SNPs but in each subgroup separately (breaks correlations)" << endl
       << "      --pbf\twhich BF to use as the test statistic for the joint-analysis permutations" << endl
       << "\t\t'none' (default): no permutations are done for the joint analysis" << endl
       << "\t\t'gen': general BF (see --bfs above)" << endl
       << "\t\t'gen-sin': 0.5 BFgen + 0.5 BFsin (also called BF_BMAlite)" << endl
       << "\t\t'all': average over all configurations (also called BF_BMA)" << endl
       << "      --maxbf\tuse the maximum ABF over SNPs as test statistic for permutations" << endl
       << "\t\totherwise the average ABF over SNPs is used (more Bayesian)" << endl
       << "      --thread\tnumber of threads (default=1, parallelize over SNPs)" << endl
       << "      --snp\tfile with a list of SNPs to analyze" << endl
       << "\t\tone SNP name per line, useful when launched in parallel" << endl
       << "\t\tprogram exits if an empty file is given" << endl
       << "      --sbgrp\tidentifier of the subgroup to analyze" << endl
       << "\t\tuseful for quick analysis and debugging" << endl
       << "\t\tcan be 'sbgrp1+sbgrp3' for instance" << endl
       << "      --wrtsize\tnumber of genes which results are written at once (default=10)" << endl
       << "\t\tto prevent excessive memory usage" << endl
       << "\t\ttune it depending on the average number of cis SNPs per gene" << endl
       << endl;
}

/** \brief Display version and license information on stdout.
 */
void version(char ** argv)
{
  cout << argv[0] << " " << VERSION << endl
       << endl
       << "Copyright (C) 2012-2014 Timothée Flutre and Xiaoquan Wen." << endl
       << "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>" << endl
       << "This is free software; see the source for copying conditions.  There is NO" << endl
       << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << endl
       << endl
       << "Written by Timothée Flutre and Xiaoquan Wen." << endl;
}

/** \brief Parse the command-line arguments and check the values of the
 *  compulsory ones.
 */
void
parseCmdLine(
  int argc,
  char ** argv,
  string & file_genopaths,
  string & snpCoordsFile,
  string & file_exppaths,
  string & file_genecoords,
  string & anchor,
  size_t & radius,
  string & file_sstats,
  string & out_prefix,
  bool & save_sstats,
  bool & save_weighted_abfs,
  string & likelihood,
  string & analysis,
  bool & need_qnorm,
  float & min_maf,
  string & covarFile,
  string & file_largegrid,
  string & file_smallgrid,
  string & file_prior,
  string & file_prior_scalar,
  string & bfs,
  string & error_model,
  float & prop_cov_errors,
  size_t & nb_types,
  size_t & nb_permutations,
  size_t & seed,
  int & trick,
  size_t & trick_cutoff,
  int & perm_sep,
  string & permbf,
  bool & use_max_bf,
  int & nb_threads,
  string & file_snpstokeep,
  vector<string> & subgroups_tokeep,
  int & write_group_size,
  int & verbose)
{
  int c = 0;
  static struct option long_options[] = {
      {"help", no_argument, 0, 'h'},
      {"version", no_argument, 0, 'V'},
      {"verbose", required_argument, 0, 'v'},
      {"geno", required_argument, 0, 0},
      {"scoord", required_argument, 0, 0},
      {"exp", required_argument, 0, 0},
      {"gcoord", required_argument, 0, 0},
      {"anchor", required_argument, 0, 0},
      {"cis", required_argument, 0, 0},
      {"inss", required_argument, 0, 0},
      {"out", required_argument, 0, 0},
      {"outss", no_argument, 0, 0},
      {"outm", no_argument, 0, 0},
      {"outw", no_argument, 0, 0},
      {"lik", required_argument, 0, 0},
      {"analys", required_argument, 0, 0},
      {"qnorm", no_argument, 0, 0},
      {"maf", required_argument, 0, 0},
      {"covar", required_argument, 0, 0},
      {"gridL", required_argument, 0, 0},
      {"gridS", required_argument, 0, 0},
      {"priorM", required_argument, 0, 0},
      {"gridM", required_argument, 0, 0},
      {"bfs", required_argument, 0, 0},
      {"error", required_argument, 0, 0},
      {"fiterr", required_argument, 0, 0},
      {"nperm", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"trick", required_argument, 0, 0},
      {"tricut", required_argument, 0, 0},
      {"permsep", required_argument, 0, 0},
      {"pbf", required_argument, 0, 0},
      {"maxbf", no_argument, 0, 0},
      {"thread", required_argument, 0, 0},
      {"snp", required_argument, 0, 0},
      {"sbgrp", required_argument, 0, 0},
      {"wrtsize", required_argument, 0, 0},
      {0, 0, 0, 0}
    };
  while(true){
    int option_index = 0;
    c = getopt_long(argc, argv, "hVv:",
		    long_options, &option_index);
    if(c == -1)
      break;
    switch (c){
    case 0:
      if(long_options[option_index].flag != 0)
	break;
      if(strcmp(long_options[option_index].name, "help") == 0){
	help(argv);
	exit(0);
      }
      if(strcmp(long_options[option_index].name, "version") == 0){
	version(argv);
	exit(0);
      }
      if(strcmp(long_options[option_index].name, "verbose") == 0){
	verbose = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "geno") == 0){
	file_genopaths = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "scoord") == 0){
	snpCoordsFile = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "exp") == 0){
	file_exppaths = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gcoord") == 0){
	file_genecoords = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "anchor") == 0){
	anchor = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "cis") == 0){
	radius = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "inss") == 0){
	file_sstats = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "out") == 0){
	out_prefix = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "outss") == 0){
	save_sstats = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "outw") == 0){
	save_weighted_abfs = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "lik") == 0){
	likelihood = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "analys") == 0){
	analysis = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "qnorm") == 0){
	need_qnorm = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "maf") == 0){
	min_maf = atof(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "covar") == 0){
	covarFile = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gridL") == 0){
	file_largegrid = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gridS") == 0){
	file_smallgrid = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "priorM") == 0){
	file_prior = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "gridM") == 0){
	file_prior_scalar = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "bfs") == 0){
	bfs = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "error") == 0){
	error_model = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "fiterr") == 0){
	prop_cov_errors = atof(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "nperm") == 0){
	nb_permutations = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "seed") == 0){
	seed = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "trick") == 0){
	trick = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "tricut") == 0){
	trick_cutoff = atol(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "permsep") == 0){
	perm_sep = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "pbf") == 0){
	permbf = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "maxbf") == 0){
	use_max_bf = true;
	break;
      }
      if(strcmp(long_options[option_index].name, "thread") == 0){
	nb_threads = atoi(optarg);
	break;
      }
      if(strcmp(long_options[option_index].name, "snp") == 0){
	file_snpstokeep = optarg;
	break;
      }
      if(strcmp(long_options[option_index].name, "sbgrp") == 0){
	split(optarg, "+", subgroups_tokeep);
	break;
      }
      if(strcmp(long_options[option_index].name, "wrtsize") == 0){
	write_group_size = atoi(optarg);
	break;
      }
    case 'h':
      help(argv);
      exit(0);
    case 'V':
      version(argv);
      exit(0);
    case 'v':
      verbose = atoi(optarg);
      break;
    case '?':
      printf("\n"); help(argv);
      abort();
    default:
      printf("\n"); help(argv);
      abort();
    }
  }

  if(file_sstats.empty()){
    if(file_genopaths.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --geno" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_genopaths)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_genopaths << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(! snpCoordsFile.empty() && ! doesFileExist(snpCoordsFile)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << snpCoordsFile << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(file_exppaths.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --exp" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_exppaths)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_exppaths << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(file_genecoords.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: missing compulsory option --gcoord" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(! doesFileExist(file_genecoords)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_genecoords << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(anchor.empty()){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: SNPs in trans not yet implemented, see --anchor and --cis" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(anchor != "TSS" && anchor != "TSS+TES"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --anchor should be TSS or TSS+TES" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
  } else{ // ! file_stats.empty()
    if(! doesFileExist(file_sstats)){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: can't find " << file_sstats << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(analysis != "join"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --inss requires --analys join" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
    if(error_model != "uvlr"){
      cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	   << "ERROR: --inss requires --error uvlr" << endl << endl;
      // help(argv);
      exit(EXIT_FAILURE);
    }
  }
  if(out_prefix.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --out" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if (write_group_size < 1){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --wrtsize should be greater than 1" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }

  // check writing files is possible
  stringstream ssOutFile;
  ssOutFile << out_prefix << "_test.txt.gz";
  gzFile outStream;
  openFile(ssOutFile.str(), outStream, "wb");
  closeFile(ssOutFile.str(), outStream);
  remove(ssOutFile.str().c_str());
  if(likelihood.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --lik" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(likelihood != "normal" && likelihood != "poisson"
     && likelihood != "quasipoisson"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --lik " << likelihood << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(likelihood != "normal" && analysis != "join"
     && error_model != "uvlr"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --lik " << likelihood << " is not valid with --analys join"
	 << " and --error " << error_model << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --analys" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis != "sep" && analysis != "join"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --analys " << analysis << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis != "join" && file_largegrid.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --gridL with --analys join" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_largegrid.empty() && ! doesFileExist(file_largegrid)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_largegrid << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_smallgrid.empty() && ! doesFileExist(file_smallgrid)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_smallgrid << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_prior_scalar.empty() && ! doesFileExist(file_prior_scalar)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_prior_scalar << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_prior.empty() && file_prior_scalar.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --gridM with --priorM" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if( file_prior.empty() && ! file_prior_scalar.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --priorM with --gridM" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(file_prior.empty() && bfs == "customized"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: missing compulsory option --priorM when --bfs customized is used" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_prior.empty() && bfs != "customized"){
    cerr << "WARNING: --priorM will be ignored because --bfs customized is not used" << endl;
  }
  if(bfs!= "gen" && bfs != "sin" && bfs != "all" && bfs != "customized"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --bfs " << bfs << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && bfs != "gen" &&
     file_smallgrid.empty()){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --gridS is required with --analys join and --bfs " << bfs << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis != "join" && bfs == "customized") {
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --analys join is required when --bfs customized is used" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(error_model == "uvlr" && bfs == "customized") {
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --error_model uvlr cannot be used with --bfs customized" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(error_model != "uvlr" && error_model != "mvlr" &&
     error_model != "hybrid"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --error " << error_model << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && error_model == "mvlr")
    cerr << "WARNING: summary statistics per subgroup won't be saved with --error mvlr" << endl;
  if(trick != 0 && trick != 1 && trick != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --trick " << trick << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(trick != 0 && trick_cutoff > nb_permutations){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --tricut " << trick_cutoff << " is larger than --nperm "
	 << nb_permutations << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(perm_sep != 0 && perm_sep != 1 && perm_sep != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --permsep " << perm_sep << " is not valid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "sep" && nb_permutations > 0 &&
     perm_sep != 1 && perm_sep != 2){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --type sep --nperm > 0, --permsep should be '1' or '2'" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(permbf != "none" && permbf != "gen" &&
     permbf != "gen-sin" && permbf != "all"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --pbf " << permbf << " is unvalid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf == "none"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --nperm > 0, --pbf should be different than 'none'" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf != "none" &&
     bfs == "gen" && permbf != "gen"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --bfs gen --nperm > 0, --pbf should be 'gen'" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(analysis == "join" && nb_permutations > 0 && permbf != "none" &&
     bfs == "sin" && permbf == "all"){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: if --analys join --bfs sin --nperm > 0, --pbf should be 'gen' or 'gen-sin'" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(nb_threads <= 0){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: --thread " << nb_threads << " is invalid" << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(! file_snpstokeep.empty() && ! doesFileExist(file_snpstokeep)){
    cerr << "cmd-line: " << getCmdLine(argc, argv) << endl << endl
	 << "ERROR: can't find " << file_snpstokeep << endl << endl;
    // help(argv);
    exit(EXIT_FAILURE);
  }
  if(seed == string::npos)
    seed = getSeed();
}

vector<Snp*> getSnpsForChr(
  const map<string, vector<Snp*> > & mChr2VecPtSnps,
  const string & chromosome)
{
  vector<Snp*> vPtSnps;
  map<string, vector<Snp*> >::const_iterator itVecPtSnps =
    mChr2VecPtSnps.find(chromosome);
  if(itVecPtSnps != mChr2VecPtSnps.end())
    vPtSnps = itVecPtSnps->second;
  return(vPtSnps);
}

void testForAssociations(
  const bool & hasDataNotSstats,
  const map<string,vector<Snp*> > & mChr2VecPtSnps,
  const string & anchor,
  const size_t & radius,
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const string & analysis,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const PriorMatrices & iPriorM,
  const string & bfs,
  const string & error_model,
  const float & prop_cov_errors,
  const int & verbose,
  map<string,Gene>::iterator & itG_begin,
  map<string,Gene>::iterator & itG_end,
  size_t & nbAnalyzedGenes,
  size_t & nbAnalyzedPairs)
{
  for(map<string,Gene>::iterator itG = itG_begin;
      itG != itG_end; ++itG){
    if(verbose > 1)
      cout << "gene " << itG->first << " (chr "
           << itG->second.GetChromosome() << ")" << endl;

    if(hasDataNotSstats){
      if(verbose > 1){
        vector<Snp*> vPtSnps = getSnpsForChr(mChr2VecPtSnps,
                                             itG->second.GetChromosome());
        size_t nbSnps = vPtSnps.size();
        cout << "chr " << itG->second.GetChromosome() << ": "
             << nbSnps << " SNP" << (nbSnps > 0 ? "s" : "") << endl;
      }
      itG->second.SetCisSnps(mChr2VecPtSnps, anchor, radius);
      if(verbose > 1)
        cout << "gene '" << itG->first << "': " << itG->second.GetNbCisSnps()
             << " cis SNP(s)" << endl;
      if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup()){
        if(verbose > 1)
          cerr << "WARNING: skip gene " << itG->second.GetName()
               << " because it has no SNP in cis" << endl;
        continue;
      }
    }

    if(analysis == "join" && error_model != "uvlr"
       && ! itG->second.HasExplevelsInAllSubgroups(subgroups)){
      if(verbose > 1)
	cerr << "WARNING: skip gene " << itG->second.GetName()
	     << " because option --error " << error_model
	     << " requires expression levels in all subgroups" << endl;
      continue;
    }
    itG->second.TestForAssociations(hasDataNotSstats, subgroups, samples,
				    likelihood, analysis, need_qnorm,
            covariates, iGridL, iGridS, iPriorM, bfs,
				    error_model, prop_cov_errors, verbose-1);
    ++nbAnalyzedGenes;
    nbAnalyzedPairs += itG->second.GetNbGeneSnpPairs();
  }

}

void makePermutationsSep(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const bool & need_qnorm,
  const Covariates & covariates,
  size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const int & perm_sep,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  map<string, Gene>::iterator & itG_begin,
  map<string, Gene>::iterator & itG_end)
{
  if(perm_sep == 1){
    gsl_rng_set(rngPerm, seed);
    if(trick != 0)
      gsl_rng_set(rngTrick, seed);
    for(map<string,Gene>::iterator itG = itG_begin;
	itG != itG_end; ++itG){
      if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
	continue;
      itG->second.MakePermutationsSepAllSubgroups(subgroups, samples,
						  likelihood,
						  need_qnorm, covariates,
						  nb_permutations, trick,
						  trick_cutoff,
						  rngPerm, rngTrick);
    }
  }
  else if(perm_sep == 2){
    for(vector<string>::const_iterator it_sbgrp = subgroups.begin();
	it_sbgrp != subgroups.end(); ++it_sbgrp){
      gsl_rng_set(rngPerm, seed);
      if(trick != 0)
	gsl_rng_set(rngTrick, seed);
      for(map<string,Gene>::iterator itG = itG_begin;
	  itG != itG_end; ++itG){
	if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
	  continue;
	itG->second.MakePermutationsSepPerSubgroup(*it_sbgrp, samples,
						   likelihood,
						   need_qnorm, covariates,
						   nb_permutations, trick,
						   trick_cutoff,
						   rngPerm, rngTrick);
      }
    }
  }
}

void makePermutationsJoin(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const PriorMatrices & iPriorM,
  const string & error_model,
  const float & prop_cov_errors,
  size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const string & permbf,
  const bool & use_max_bf,
  const gsl_rng * rngPerm,
  const gsl_rng * rngTrick,
  map<string, Gene>::iterator & itG_begin,
  map<string, Gene>::iterator & itG_end)
{
  gsl_rng_set(rngPerm, seed);
  if(trick != 0)
    gsl_rng_set(rngTrick, seed);

  for(map<string,Gene>::iterator itG = itG_begin;
       itG != itG_end; ++itG){
    if(! itG->second.HasAtLeastOneCisSnpInAtLeastOneSubgroup())
      continue;
    if(error_model != "uvlr" &&
       ! itG->second.HasExplevelsInAllSubgroups(subgroups))
      continue;
    itG->second.MakePermutationsJoin(subgroups, samples, likelihood, need_qnorm,
                                     covariates, iGridL, iGridS, iPriorM,
                                     error_model, prop_cov_errors,
                                     nb_permutations, trick,
                                     trick_cutoff, permbf, use_max_bf,
                                     rngPerm, rngTrick);
  }
}

void
makePermutations(
  const vector<string> & subgroups,
  const Samples & samples,
  const string & likelihood,
  const string & analysis,
  const bool & need_qnorm,
  const Covariates & covariates,
  const Grid & iGridL,
  const Grid & iGridS,
  const PriorMatrices & iPriorM,
  const string & error_model,
  const float & prop_cov_errors,
  size_t & nb_permutations,
  const size_t & seed,
  const int & trick,
  const size_t & trick_cutoff,
  const int & perm_sep,
  const string & permbf,
  const bool & use_max_bf,
  map<string, Gene>::iterator & itG_begin,
  map<string, Gene>::iterator & itG_end)
{

  gsl_rng * rngPerm = NULL, * rngTrick = NULL;
  gsl_rng_env_setup();
  rngPerm = gsl_rng_alloc (gsl_rng_default);
  if(rngPerm == NULL){
    cerr << "ERROR: can't allocate memory for the RNG" << endl;
    exit(EXIT_FAILURE);
  }
  if(trick != 0){
    rngTrick = gsl_rng_alloc (gsl_rng_default);
    if(rngTrick == NULL){
      cerr << "ERROR: can't allocate memory for the RNG" << endl;
      exit(EXIT_FAILURE);
    }
  }

  if(analysis == "sep" && perm_sep != 0)
    makePermutationsSep(subgroups, samples, likelihood, need_qnorm, covariates,
			nb_permutations, seed, trick, trick_cutoff,
			perm_sep, rngPerm, rngTrick, itG_begin, itG_end);
  if(analysis == "join" && permbf != "none")
    makePermutationsJoin(subgroups, samples, likelihood, need_qnorm, covariates,
                         iGridL, iGridS, iPriorM, error_model, prop_cov_errors,
                         nb_permutations, seed, trick, trick_cutoff, permbf,
                         use_max_bf, rngPerm, rngTrick, itG_begin, itG_end);

  gsl_rng_free(rngPerm);
  if(trick != 0)
    gsl_rng_free(rngTrick);
}

void writeResSstats(
  const string & out_prefix,
  const vector<string> & subgroups,
  const map<string,Gene>::iterator & itG_begin,
  const map<string,Gene>::iterator & itG_end,
  const map<string,Snp> & snp2object,
  const string & header_opt)
{
  string subgroup, sep = "\t";

  for(size_t s = 0; s < subgroups.size(); ++s){
    subgroup = subgroups[s];
    stringstream ssOutFile, ssTxt;
    ssOutFile << out_prefix << "_sumstats_" << subgroup << ".txt.gz";
    gzFile outStream;
    string file_mode;
    if (header_opt == "none")
      file_mode = "ab";
    else
      file_mode = "wb";
    openFile(ssOutFile.str(), outStream, file_mode.c_str());
    size_t nb_lines = 1;
    if (header_opt != "none") {
      ssTxt << "gene" << sep << "snp" << sep << "maf" << sep << "n" << sep
            << "pve" << sep << "sigmahat" << sep << "betahat.geno" << sep
            << "sebetahat.geno" << sep << "betapval.geno";
      ssTxt << endl;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
    if (header_opt == "only") {
      closeFile(ssOutFile.str(), outStream);
      continue;
    }

    ssTxt.precision(6);
    ssTxt.setf(ios::scientific);
    for(map<string,Gene>::const_iterator it_gene = itG_begin;
	 it_gene != itG_end; ++it_gene){
      if(! it_gene->second.HasAtLeastOneCisSnp(subgroup))
	continue;
      for(vector<GeneSnpPair>::const_iterator it_pair
	     = it_gene->second.BeginPair();
	   it_pair != it_gene->second.EndPair(); ++it_pair){
	if(! it_pair->HasResults(subgroup))
	  continue;
	ssTxt.str("");
	++nb_lines;
	ssTxt << it_gene->second.GetName()
	      << sep << it_pair->GetSnpName()
	      << sep << snp2object.find(it_pair->GetSnpName())->second.GetMinorAlleleFreq(subgroup)
	      << sep << it_pair->GetSampleSize(subgroup)
	      << sep << it_pair->GetPve(subgroup)
	      << sep << it_pair->GetSigmahat(subgroup)
	      << sep << it_pair->GetBetahatGeno(subgroup)
	      << sep << it_pair->GetSebetahatGeno(subgroup)
	      << sep << it_pair->GetBetapvalGeno(subgroup);
	ssTxt << "\n";
	gzwriteLine (outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      } // end of loop over cis snps

    } // end of loop over genes

    closeFile(ssOutFile.str(), outStream);

  } // end of loop over subgroups
}

void writeResSepPermPval(
  const string & out_prefix,
  const map<string, Gene>::iterator & itG_begin,
  const map<string, Gene>::iterator & itG_end,
  const vector<string> & subgroups,
  const size_t & seed,
  const string & header_opt)
{
  string sep = "\t";

  for(vector<string>::const_iterator it_sbgrp  = subgroups.begin();
      it_sbgrp != subgroups.end(); ++it_sbgrp){
    stringstream ssOutFile, ssTxt;
    ssOutFile << out_prefix << "_sepPermPvals_" << *it_sbgrp << ".txt.gz";
    gzFile outStream;
    string file_mode;
    if (header_opt == "none")
      file_mode = "ab";
    else
      file_mode = "wb";
    openFile(ssOutFile.str(), outStream, file_mode.c_str());
    size_t nb_lines = 0;
    if (header_opt != "none") {
      ssTxt << "# seed=" << seed << endl;
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      ssTxt.str("");
      ssTxt << "gene" << sep << "nb.snps" << sep << "sep.perm.pval"
            << sep << "nb.permutations" << sep << "true.min.pval" << endl;
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    } else
      nb_lines = 2;
    if (header_opt == "only") {
      closeFile(ssOutFile.str(), outStream);
      continue;
    }

    for(map<string,Gene>::const_iterator it_gene = itG_begin;
	 it_gene != itG_end; ++it_gene){
      if(it_gene->second.GetNbGeneSnpPairs() > 0){
	ssTxt.str("");
	++nb_lines;
	ssTxt << it_gene->second.GetName()
	      << sep << it_gene->second.GetNbGeneSnpPairs(*it_sbgrp)
	      << sep << it_gene->second.GetPermutationPvalueSep(*it_sbgrp)
	      << sep << it_gene->second.GetNbPermutationsSep(*it_sbgrp);
	ssTxt.setf(ios::scientific);
	ssTxt << sep << setprecision(6) << it_gene->second.GetTrueMinPval(*it_sbgrp);
	ssTxt.unsetf(ios::scientific);
	ssTxt << "\n";
	gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      }
    }

    closeFile(ssOutFile.str(), outStream);
  }
}

void writeResSepPermPval(
  const string & out_prefix,
  const map<string, Gene>::iterator & itG_begin,
  const map<string, Gene>::iterator & itG_end,
  const size_t & seed,
  const string & header_opt)
{
  string sep = "\t";

  stringstream ssOutFile, ssTxt;
  ssOutFile << out_prefix << "_sepPermPvals.txt.gz";
  gzFile outStream;
  string file_mode;
  if (header_opt == "none")
    file_mode = "ab";
  else
    file_mode = "wb";
  openFile(ssOutFile.str(), outStream, file_mode.c_str());
  size_t nb_lines = 0;
  if (header_opt != "none") {
    ssTxt << "# seed=" << seed << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    ssTxt.str("");
    ssTxt << "gene" << sep << "nb.snps" << sep << "sep.perm.pval"
          << sep << "nb.permutations" << sep << "true.min.pval" << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  } else
    nb_lines = 2;
  if (header_opt == "only") {
    closeFile(ssOutFile.str(), outStream);
    return;
  }

  for(map<string,Gene>::const_iterator it_gene = itG_begin;
      it_gene != itG_end; ++it_gene){
    if(it_gene->second.GetNbGeneSnpPairs() > 0){
      ssTxt.str("");
      ++nb_lines;
      ssTxt << it_gene->second.GetName()
	    << sep << it_gene->second.GetNbGeneSnpPairs()
	    << sep << it_gene->second.GetPermutationPvalueSep()
	    << sep << it_gene->second.GetNbPermutationsSep();
      ssTxt.setf(ios::scientific);
      ssTxt << sep << setprecision(6) << it_gene->second.GetTrueMinPval();
      ssTxt.unsetf(ios::scientific);
      ssTxt << "\n";
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }

  closeFile(ssOutFile.str(), outStream);
}

/** \brief
 *  \note the 'const' BFs use the large grid while the other use the small one,
 *  but the header line lists the large grid only
 */
void writeResAbfsRaw(
  const string & out_prefix,
  const map<string,Gene>::iterator & itG_begin,
  const map<string,Gene>::iterator & itG_end,
  const size_t & nb_subgroups,
  const Grid & iGridL,
  const Grid & iGridS,
  const PriorMatrices & iPriorM,
  const string & bfs,
  const string & header_opt)
{
  string sep = "\t";

  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << out_prefix << "_l10abfs_raw.txt.gz";
  gzFile outStream;
  string file_mode;
  if (header_opt == "none")
    file_mode = "ab";
  else
    file_mode = "wb";
  openFile(ssOutFile.str(), outStream, file_mode.c_str());
  size_t nb_lines = 1;
  if (header_opt != "none") {
    // write header line
    ssTxt << "gene" << sep << "snp" << sep << "config";
    for(size_t i = 0; i < iGridL.size(); ++i)
      ssTxt << sep << "l10abf.grid" << (i+1);
    ssTxt << endl;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  }
  if (header_opt == "only") {
    closeFile(ssOutFile.str(), outStream);
    return;
  }

  // write results
  ssTxt.precision(6);
  ssTxt.setf(ios::scientific);
  for(map<string,Gene>::const_iterator it_gene = itG_begin;
      it_gene != itG_end; ++it_gene){
    for(vector<GeneSnpPair>::const_iterator it_pair
	   = it_gene->second.BeginPair(); it_pair != it_gene->second.EndPair();
	 ++it_pair){

      // write gen BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen");
	  it != it_pair->EndUnweightedAbf("gen"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);

      // write gen-fix BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen-fix";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen-fix");
	  it != it_pair->EndUnweightedAbf("gen-fix"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);

      // write gen-maxh BFs (large grid)
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << "gen-maxh";
      for(vector<double>::const_iterator it
	    = it_pair->BeginUnweightedAbf("gen-maxh");
	  it != it_pair->EndUnweightedAbf("gen-maxh"); ++it)
	ssTxt << sep << *it;
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);

      // write the BFs for each config (small grid)
      if(bfs != "gen"){
      	for(size_t k = 1; k <= nb_subgroups; ++k){
      	  comb = gsl_combination_calloc(nb_subgroups, k);
      	  if(comb == NULL){
      	    cerr << "ERROR: can't allocate memory for the combination"
      		 << endl;
      	    exit(EXIT_FAILURE);
      	  }
      	  while(true){
      	    ssTxt.str("");
      	    ssTxt << it_gene->first
      		  << sep << it_pair->GetSnpName();
      	    ssConfig.str("");
      	    ssConfig << gsl_combination_get(comb, 0) + 1;
      	    if(comb->k > 1)
      	      for(size_t i = 1; i < k; ++i)
      		ssConfig << "-" << gsl_combination_get(comb, i) + 1;
      	    ssTxt << sep << ssConfig.str();
      	    for(size_t j = 0; j < iGridL.size(); ++j){
      	      if(j < iGridS.size())
      		ssTxt << sep << *(it_pair->BeginUnweightedAbf(ssConfig.str()) + j);
      	      else
      		ssTxt << sep << NaN;
      	    }
      	    ssTxt << "\n";
      	    ++nb_lines;
      	    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
      	    if(gsl_combination_next(comb) != GSL_SUCCESS)
      	      break;
      	  }
      	  gsl_combination_free(comb);
          if(bfs == "sin" || bfs == "customized")
      	    break;
      	}
      }
      // write the BFs for each customized prior (with customized grid)
      if(bfs == "customized") {
        for (size_t m = 0; m < iPriorM.Wg_names.size(); ++m) {
          ssTxt.str("");
      	  ssTxt << it_gene->first
                << sep << it_pair->GetSnpName();
          ssTxt << sep << iPriorM.Wg_names[m];
          for (size_t j = 0; j < iGridL.size(); ++j) {
            if (j < iPriorM.Wg_scalars.size())
      		ssTxt << sep << *(it_pair->BeginUnweightedAbf(iPriorM.Wg_names[m]) + j);
      	      else
      		ssTxt << sep << NaN;
          }
      	  ssTxt << "\n";
      	  ++nb_lines;
      	  gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
        }
      }
    }
  }

  closeFile(ssOutFile.str(), outStream);
}

void writeResAbfsAvgGrids(
  const string & out_prefix,
  const map<string, Gene>::iterator & itG_begin,
  const map<string, Gene>::iterator & itG_end,
  const size_t & nb_subgroups,
  const vector<string> & customized_priors,
  const string & bfs,
  const string & error_model,
  const string & header_opt)
{
  string sep = "\t";

  gsl_combination * comb;
  stringstream ssOutFile, ssConfig, ssTxt;
  ssOutFile << out_prefix << "_l10abfs_avg-grids.txt.gz";
  gzFile outStream;
  string file_mode;
  if (header_opt == "none")
    file_mode = "ab";
  else
    file_mode = "wb";
  openFile(ssOutFile.str(), outStream, file_mode.c_str());
  size_t nb_lines = 1;
  if (header_opt != "none") {
    // write header line
    ssTxt << "gene" << sep << "snp" << sep << "nb.subgroups";
    ssTxt << sep << "l10abf.gen"
          << sep << "l10abf.gen.fix"
          << sep << "l10abf.gen.maxh";
    if(bfs == "sin" || bfs == "all" || bfs == "customized")
      ssTxt << sep << "l10abf.gen.sin";
    if(bfs == "all")
      ssTxt << sep << "l10abf.all";
    if(bfs == "customized")
      ssTxt << sep << "l10abf.customized";
    if(bfs != "gen"){
      for(size_t k = 1; k <= nb_subgroups; ++k){
        comb = gsl_combination_calloc(nb_subgroups, k);
        if(comb == NULL){
          cerr << "ERROR: can't allocate memory for the combination" << endl;
          exit(EXIT_FAILURE);
        }
        while(true){
          ssTxt << sep << "l10abf." << gsl_combination_get(comb, 0) + 1;
          if(comb->k > 1)
            for(size_t i = 1; i < k; ++i)
              ssTxt << "-" << gsl_combination_get(comb, i) + 1;
          if(gsl_combination_next(comb) != GSL_SUCCESS)
            break;
        }
        gsl_combination_free(comb);
        if(bfs == "sin" || bfs == "customized")
          break;
      }
    }
    if(bfs == "customized") {
      for (size_t m = 0; m < customized_priors.size(); ++m) {
        ssTxt << sep << "l10abf." << customized_priors[m];
      }
    }
    ssTxt << endl;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  }
  if (header_opt == "only") {
    closeFile(ssOutFile.str(), outStream);
    return;
  }

  // write results
  ssTxt.precision(6);
  ssTxt.setf(ios::scientific);
  for(map<string, Gene>::const_iterator it_gene = itG_begin;
       it_gene != itG_end; ++it_gene){
    for(vector<GeneSnpPair>::const_iterator it_pair
	   = it_gene->second.BeginPair(); it_pair != it_gene->second.EndPair();
	 ++it_pair){
      ssTxt.str("");
      ssTxt << it_gene->first
	    << sep << it_pair->GetSnpName()
	    << sep << it_pair->GetNbSubgroups()
	    << sep << it_pair->GetWeightedAbf("gen")
	    << sep << it_pair->GetWeightedAbf("gen-fix")
	    << sep << it_pair->GetWeightedAbf("gen-maxh");
      if(bfs == "sin" || bfs == "all" || bfs == "customized")
	ssTxt << sep << it_pair->GetWeightedAbf("gen-sin");
      if(bfs == "all")
	ssTxt << sep << it_pair->GetWeightedAbf("all");
      if(bfs == "customized")
	ssTxt << sep << it_pair->GetWeightedAbf("customized");
      if(bfs != "gen"){
	for(size_t k = 1; k <= nb_subgroups; ++k){
	  comb = gsl_combination_calloc(nb_subgroups, k);
	  if(comb == NULL){
	    cerr << "ERROR: can't allocate memory for the combination"
		 << endl;
	    exit(EXIT_FAILURE);
	  }
	  while(true){
	    ssConfig.str("");
	    ssConfig << gsl_combination_get(comb, 0) + 1;
	    if(comb->k > 1)
	      for(size_t i = 1; i < k; ++i)
		ssConfig << "-" << gsl_combination_get(comb, i) + 1;
	    ssTxt << sep << it_pair->GetWeightedAbf(ssConfig.str());
	    if(gsl_combination_next(comb) != GSL_SUCCESS)
	      break;
	  }
	  gsl_combination_free(comb);
	  if(bfs == "sin" || bfs == "customized")
	    break;
	}
      }
      if (bfs == "customized") {
        for (size_t m = 0; m < customized_priors.size(); ++m) {
          ssTxt << sep << it_pair->GetWeightedAbf(customized_priors[m]);
        }
      }
      ssTxt << "\n";
      ++nb_lines;
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }

  closeFile(ssOutFile.str(), outStream);
}

void writeResJoinPermPval(
  const string & out_prefix,
  const map<string,Gene>::iterator & itG_begin,
  const map<string,Gene>::iterator & itG_end,
  const size_t & seed,
  const string & permbf,
  const bool & use_max_bf,
  const string & header_opt)
{
  string sep = "\t";

  stringstream ssOutFile, ssTxt;
  ssOutFile << out_prefix << "_joinPermPvals.txt.gz";
  gzFile outStream;
  string file_mode;
  if (header_opt == "none")
    file_mode = "ab";
  else
    file_mode = "wb";
  openFile(ssOutFile.str(), outStream, file_mode.c_str());
  size_t nb_lines = 0;
  if (header_opt != "none") {
    ssTxt << "# perm.bf=" << permbf << " seed=" << seed << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);

    ssTxt.str("");
    ssTxt << "gene" << sep << "nb.snps" << sep << "join.perm.pval"
          << sep << "nb.permutations" << sep << "true.l10abf"
          << sep << "med.perm.l10abf" << endl;
    ++nb_lines;
    gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
  } else
    nb_lines = 2;
  if (header_opt == "only") {
    closeFile(ssOutFile.str(), outStream);
    return;
  }

  for(map<string,Gene>::const_iterator it_gene = itG_begin;
      it_gene != itG_end; ++it_gene){
    if(it_gene->second.GetNbGeneSnpPairs() > 0){
      ssTxt.str("");
      ++nb_lines;
      ssTxt << it_gene->second.GetName()
	    << sep << it_gene->second.GetNbGeneSnpPairs()
	    << sep << it_gene->second.GetPermutationPvalueJoin()
	    << sep << it_gene->second.GetNbPermutationsJoin();
      ssTxt.setf(ios::scientific);
      ssTxt << setprecision(6) << sep << it_gene->second.GetTrueL10Abf(use_max_bf)
	    << sep << it_gene->second.GetMedianPermL10Abf();
      ssTxt.unsetf(ios::scientific);
      ssTxt << "\n";
      gzwriteLine(outStream, ssTxt.str(), ssOutFile.str(), nb_lines);
    }
  }

  closeFile(ssOutFile.str(), outStream);
}

void writeRes(
  const string & out_prefix,
  const bool & save_sstats,
  const bool & save_weighted_abfs,
  const vector<string> & subgroups,
  const map<string, Gene>::iterator & itG_begin,
  const map<string, Gene>::iterator & itG_end,
  const map<string, Snp> & snp2object,
  const string & analysis,
  const Grid & iGridL,
  const Grid & iGridS,
  const PriorMatrices & iPriorM,
  const string & bfs,
  const size_t & nb_permutations,
  const int & perm_sep,
  const string & error_model,
  const size_t & seed,
  const string & permbf,
  const bool & use_max_bf,
  const string & header_opt)
{
  if(analysis == "sep" ||
     (analysis == "join" && save_sstats &&
      error_model != "mvlr"))
    writeResSstats(out_prefix, subgroups, itG_begin, itG_end, snp2object,
                   header_opt);

  if(analysis == "sep" && nb_permutations > 0 &&
     perm_sep != 0){
    if(perm_sep == 1)
      writeResSepPermPval(out_prefix, itG_begin, itG_end, seed, header_opt);
    else if(perm_sep == 2)
      writeResSepPermPval(out_prefix, itG_begin, itG_end, subgroups, seed,
                          header_opt);
  }


  if(analysis == "join"){
    writeResAbfsRaw(out_prefix, itG_begin, itG_end, subgroups.size(),
                    iGridL, iGridS, iPriorM, bfs, header_opt);
    if(save_weighted_abfs)
      writeResAbfsAvgGrids(out_prefix, itG_begin, itG_end, subgroups.size(),
                           iPriorM.Wg_names, bfs, error_model, header_opt);
  }

  if(analysis == "join" && nb_permutations > 0)
    writeResJoinPermPval(out_prefix, itG_begin, itG_end, seed, permbf,
                         use_max_bf, header_opt);
}

void run(const string & file_genopaths,
	 const string & file_snpcoords,
	 const string & file_exppaths,
	 const string & file_genecoords,
	 const string & anchor,
	 const size_t & radius,
	 const string & file_sstats,
	 const string & out_prefix,
	 const bool & save_sstats,
	 const bool & save_weighted_abfs,
	 const string & likelihood,
	 const string & analysis,
	 const bool & need_qnorm,
	 const float & min_maf,
	 const string & file_covarpaths,
	 const string & file_largegrid,
	 const string & file_smallgrid,
   const string & file_prior,
   const string & file_prior_scalar,
	 const string & bfs,
	 const string & error_model,
	 const float & prop_cov_errors,
	 size_t & nb_permutations,
	 const size_t & seed,
	 const int & trick,
	 const size_t & trick_cutoff,
	 const int & perm_sep,
	 const string & permbf,
	 const bool & use_max_bf,
	 const int & nb_threads,
	 const string & file_snpstokeep,
	 const vector<string> & subgroups_tokeep,
     const int & write_group_size,
     const int & verbose)
{
  set<string> sSnpsToKeep;
  if(! file_snpstokeep.empty()){
    loadSnpsToKeep(file_snpstokeep, verbose, sSnpsToKeep);
    if(sSnpsToKeep.empty())
      return;
  }

  vector<string> subgroups;
  Samples samples;
  map<string,Snp> snp2object;
  map<string,vector<Snp*> > mChr2VecPtSnps;
  Covariates covariates;
  map<string,Gene> gene2object;

  if(file_sstats.empty())
    loadRawInputData(file_genopaths, file_snpcoords, file_exppaths,
		     file_genecoords, anchor, radius, min_maf, file_covarpaths,
		     error_model, subgroups_tokeep, sSnpsToKeep, verbose,
		     subgroups, samples, snp2object, mChr2VecPtSnps,
		     covariates, gene2object);
  else
    loadSummaryStats(file_sstats, verbose, subgroups, gene2object, snp2object);
  if(gene2object.empty() || snp2object.empty())
    return;

  // load grid / customized priors
  Grid iGridL(file_largegrid, true, verbose);
  Grid iGridS(file_smallgrid, false, verbose);
  PriorMatrices iPriorM(file_prior, file_prior_scalar, subgroups.size(), verbose);
  if (iPriorM.Wg_scalars.size() > iGridL.size()) {
    cerr << "WARNING: --gridM is larger than --gridL; thus will be truncated to the size of --gridL" << endl;
  }

  // write header
  writeRes(out_prefix, save_sstats, save_weighted_abfs, subgroups,
           gene2object.begin(), gene2object.end(),
           snp2object, analysis, iGridL, iGridS, iPriorM, bfs, nb_permutations,
           perm_sep, error_model, seed, permbf, use_max_bf, "only");

  size_t countGenes = 0;
  size_t totalGenes = gene2object.size();
  bool is_perm = nb_permutations > 0 && (perm_sep != 0 || permbf != "none");
  if(verbose > 0){
    cout << "test for association between each pair gene-SNP ..." << endl
	 << "analysis=" << analysis
	 << " likelihood=" << likelihood
	 << " error_model=" << error_model;
    if(error_model != "uvlr") // i.e. if 'mvlr' or 'hybrid'
      cout << " prop_cov_errors=" << prop_cov_errors;
    if(file_sstats.empty())
      cout << " anchor=" << anchor << " radius=" << radius;
    if(is_perm) {
      cout << endl << "permutation"<< (nb_permutations > 1 ? "s=" : "=")
           << nb_permutations << " seed=" << seed;
      if(trick != 0){
        cout << " trick=" << trick
             << " trick_cutoff=" << trick_cutoff;
      }
      if(analysis == "sep")
        cout << " perm_sep=" << perm_sep;
      else if(analysis == "join")
        cout << " perm_bf=" << permbf;
      cout << " threads=" << omp_get_max_threads();
    }
    cout << endl << flush;
  }
  clock_t startTime = clock();
  size_t nbAnalyzedGenes = 0, nbAnalyzedPairs = 0;
  if(verbose == 1)
    progressBar("", 0, totalGenes);
  for(map<string,Gene>::iterator itG = gene2object.begin();
      itG != gene2object.end();) {
    map<string,Gene>::iterator itG_begin = itG;
    size_t step_size = min((int)distance(itG, gene2object.end()),
                           write_group_size);
    advance(itG, step_size);
    testForAssociations(file_sstats.empty(), mChr2VecPtSnps, anchor, radius,
                        subgroups, samples, likelihood, analysis,
                        need_qnorm, covariates, iGridL, iGridS, iPriorM, bfs,
                        error_model, prop_cov_errors, verbose,
                        itG_begin, itG, nbAnalyzedGenes, nbAnalyzedPairs);
    if(is_perm){
      omp_set_num_threads(nb_threads);
      makePermutations(subgroups, samples, likelihood, analysis, need_qnorm,
                       covariates, iGridL, iGridS, iPriorM, error_model,
                       prop_cov_errors,
                       nb_permutations, seed, trick, trick_cutoff, perm_sep,
                       permbf, use_max_bf, itG_begin, itG);
    }
    // write results
    writeRes(out_prefix, save_sstats, save_weighted_abfs, subgroups,
             itG_begin, itG, snp2object, analysis, iGridL, iGridS, iPriorM,
             bfs, nb_permutations, perm_sep, error_model, seed, permbf,
             use_max_bf, "none");
    // progress tracker
    countGenes += step_size;
    if(verbose == 1)
      progressBar("", countGenes, totalGenes);
    // delete processed object to release RAM
    gene2object.erase(itG_begin, itG);
  }
  if(verbose > 0){
    if(verbose == 1)
      cout << " (" << fixed << setprecision(2) << getElapsedTime(startTime)
	   << " sec)" << endl << flush;
    cout << "nb of analyzed gene-SNP pairs: " << nbAnalyzedPairs
	 << " (" << nbAnalyzedGenes << " genes)" << endl;
  }
}

int main(int argc, char ** argv)
{
  int verbose = 1, trick = 0, perm_sep = 0, nb_threads = 1,
    write_group_size = 10;
  size_t radius = 100000, nb_types = string::npos, nb_permutations = 0,
    seed = string::npos, trick_cutoff = 10;
  float min_maf = 0.0, prop_cov_errors = 0.5;
  bool save_sstats = false, save_weighted_abfs = false, need_qnorm = false,
    use_max_bf = false;
  string file_genopaths, file_snpcoords, file_exppaths, file_genecoords,
    anchor = "TSS", file_sstats, out_prefix, likelihood = "normal",
    analysis, file_covarpaths, file_largegrid, file_smallgrid,
    file_prior, file_prior_scalar, bfs = "gen", error_model = "uvlr",
    permbf = "none", file_snpstokeep;
  vector<string> subgroups_tokeep;

  parseCmdLine(argc, argv, file_genopaths, file_snpcoords, file_exppaths,
	       file_genecoords, anchor, radius, file_sstats, out_prefix,
	       save_sstats, save_weighted_abfs, likelihood, analysis,
	       need_qnorm, min_maf, file_covarpaths, file_largegrid,
         file_smallgrid, file_prior, file_prior_scalar, bfs, error_model,
         prop_cov_errors, nb_types,
	       nb_permutations, seed, trick, trick_cutoff, perm_sep,
         permbf, use_max_bf, nb_threads, file_snpstokeep,
	       subgroups_tokeep, write_group_size, verbose);

  time_t time_start, time_end;
  if(verbose > 0){
    time(&time_start);
    cout << "START " << basename(argv[0])
	 << " " << getDateTime(time_start) << endl
	 << "version " << VERSION << " compiled " << __DATE__
	 << " " << __TIME__ << endl
	 << "cmd-line: " << getCmdLine(argc, argv) << endl
	 << "cwd: " << getCurrentDirectory() << endl << flush;
  }

  run(file_genopaths, file_snpcoords, file_exppaths, file_genecoords, anchor,
      radius, file_sstats, out_prefix, save_sstats, save_weighted_abfs,
      likelihood, analysis, need_qnorm, min_maf, file_covarpaths,
      file_largegrid, file_smallgrid, file_prior, file_prior_scalar,
      bfs, error_model, prop_cov_errors,
      nb_permutations, seed, trick, trick_cutoff, perm_sep, permbf,
      use_max_bf, nb_threads, file_snpstokeep, subgroups_tokeep,
      write_group_size, verbose);

  if(verbose > 0){
    time (&time_end);
    cout << "END " << basename(argv[0])
	 << " " << getDateTime(time_end) << endl
	 << "elapsed -> " << getElapsedTime(time_start, time_end) << endl
	 << "max.mem -> " << getMaxMemUsedByProcess2Str() << endl;
  }

  return EXIT_SUCCESS;
}
