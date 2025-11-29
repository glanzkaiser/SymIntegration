/*
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

/*
 * zlib License
 *
 * Regularized Incomplete Beta Function
 *
 * Copyright (c) 2016, 2017 Lewis Van Winkle
 * http://CodePlea.com
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgement in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

#include <random>
#include <iostream>
#include <vector>

#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DECLARE

int fibonacciseries(int);
double rising_pochhammer(double, int);
double falling_pochhammer(double, int);

Symbolic divisions(Symbolic, Symbolic);
Symbolic factorialsym(Symbolic);
Symbolic factorials(int);
Symbolic combinationss(int, int);
double divisionint(double, double);

double factoriald(int);
double combinationsd(int, int);

vector<string> loadStringVector(const string&);
void printStringVector(vector<string>);
void saveVectordouble(vector<double>, const string&);

double binomialpmf(int, int, double);
double binomialcdf(int, int, double);
double binomialmean(int, int, double);
double binomialvar(int, int, double);
Symbolic binomialmgf(int, int, double);

double negativebinomialpmf(int, int, double);
double negativebinomialcdf(int, int, double);
double negativebinomialmean(int, int, double);
double negativebinomialvar(int, int, double);
Symbolic negativebinomialmgf(const Symbolic &, const Symbolic &, const Symbolic &);

double geometricpmf(int, double);
double geometricmean(int, double);
double geometricvar(int, double);
Symbolic geometricmgf(const Symbolic &, const Symbolic &);

double hypergeometricpmf(int, int, int, int);
double hypergeometricmean(int, int, int, int);
double hypergeometricvar(int, int, int, int);

double poissonpmf(int, int);
double poissoncdf(int, int);
double poissonmean(int, int);
double poissonvar(int, int);
Symbolic poissonmgf(const Symbolic &, const Symbolic &);

double uniformpdf(double, double, double);
double uniformcdf(double, double, double);
Symbolic uniformmgf(double, double, double);
double uniformmean(double, double, double);
double uniformvar(double, double, double);

Symbolic normalpdf(double, double, double);
double normalcdf(double, double, double);
Symbolic normalmgf(double, double, double);
double normalmean(double, double, double);
double normalvar(double, double, double);

double zquantile(double);
double tquantile(double, double, double);
double Fquantile(double, double, double);
double chisquaredquantile(double, double);

double gammapdf(double, double, double);
double gammacdf(double, double, double);
Symbolic gammamgf(double, double, double);
double gammamean(double, double, double);
double gammavar(double, double, double);

Symbolic exponentialpdf(double, double);
Symbolic exponentialcdf(double, double);
Symbolic exponentialmgf(double, double);
double exponentialmean(double, double);
double exponentialvar(double, double);

double betapdf(double, double, double);
double betacdf(double, double, double);
Symbolic betamgf(double, double, double);
double betamean(double, double, double);
double betavar(double, double, double);

double chisquaredpdf(double, double);
double chisquaredcdf(double, double);
Symbolic chisquaredmgf(double, double);
double chisquaredmean(double, double);
double chisquaredvar(double, double);

double cauchypdf(double);

double lowergamma(double, double);
double incbeta(double, double, double);

double Fpdf(double, double, double);
double Fcdf(double, double, double);
double Fmean(double, double, double);
double Fvar(double, double, double);

double tpdf(double, double);
double tcdf(double, double);
double tmean(double, double);
double tvar(double, double);

Symbolic laplacepdf(double, double);
Symbolic laplacecdf(double, double);
Symbolic laplacemgf(double, double);
Symbolic laplacemean(double, double);
Symbolic laplacevar(double, double);

Symbolic logisticpdf(double, double);
Symbolic logisticcdf(double, double);
Symbolic logisticmgf(double, double);
Symbolic logisticmean(double, double);
Symbolic logisticvar(double, double);

double descriptivestatistics(vector<double>);

double rpearson(const SymbolicMatrix&, int);
Symbolic regressionline(const SymbolicMatrix&, int);
vector<vector<double>> multipleregression(vector<vector<double>> &, vector<vector<double>> &);

int randomnumberint(int, int, int);
double randomnumberreal(double, double, int);

double randomnumbernormal(double, double, int);
double randomnumbergamma(double, double, int);
double randomnumberbeta(double, double, int);

double generateBeta(double, double, std::mt19937&);
std::vector<double> vrandn_bernoulli(double, int);
std::vector<double> vrandn_binomial(double, int);
std::vector<double> vrandn_normal(double, double, int);
std::vector<double> vrandn_exponential(double, int);
std::vector<double> vrandn_chisquared(double, int);
std::vector<double> vrandn_gamma(double, double, int);
std::vector<double> vrandn_beta(double, double, int);
std::vector<double> vrandn_fdist(double, double, int);
std::vector<double> vrandn_tdist(double, int);
std::vector<double> vrandn_erlang(double, double, int);

double calculateMean(vector<double>);
double calculateCovariance(vector<double>, vector<double>); 
SymbolicMatrix covariancematrix(vector<vector<double>>);

vector<double> normalquantile(vector<double>, double);

void confidenceinterval_onesampletwosides(vector<double>, double, double);
void confidenceinterval_onesampleonesided(vector<double>, double, double);
void confidenceinterval_onesampletwosides(vector<double>, double);
void confidenceinterval_onesampleonesided(vector<double>, double);
void predictioninterval(vector<double>, double, double);
void predictioninterval(vector<double>, double);

void confidenceinterval_knownsigma(vector<double>, vector<double>, double, double, double);
void confidenceinterval_sameunknownsigma(vector<double>, vector<double>, double);
void confidenceinterval_unequalunknownsigma(vector<double>, vector<double>, double);
void confidenceinterval_paired(vector<double>, vector<double>, double);
void confidenceinterval_proportion(vector<string>, const string&, double);
void confidenceinterval_differencebetweentwoproportions(vector<string>, vector<string>, const string&, double);
void confidenceinterval_variance(vector<double>, double);
void confidenceinterval_ratiotwovariances(vector<double>, vector<double>, double);

void hypothesistest(vector<double>, double, double, double, double);
void hypothesistest_lefttailed(vector<double>, double, double, double, double);
void hypothesistest_righttailed(vector<double>, double, double, double, double);

void hypothesistest_knownvariances_righttailed(vector<double>, vector<double>, double, double, double, double, double);
void hypothesistest_knownvariances_lefttailed(vector<double>, vector<double>, double, double, double, double, double);
void hypothesistest_knownvariances_twotailed(vector<double>, vector<double>, double, double, double, double, double);
void hypothesistest_equalunknownvariances_righttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_equalunknownvariances_lefttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_equalunknownvariances_twotailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_unequalunknownvariances_righttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_unequalunknownvariances_lefttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_unequalunknownvariances_twotailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_paired_righttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_paired_lefttailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_paired_twotailed(vector<double>, vector<double>, double, double, double);
void hypothesistest_proportion_righttailed(vector<string>, const string&, double, double);
void hypothesistest_proportion_lefttailed(vector<string>, const string&, double, double);
void hypothesistest_proportion_twotailed(vector<string>, const string&, double, double);
void hypothesistest_twoproportions_righttailed(vector<string>, vector<string>, const string&, double);
void hypothesistest_twoproportions_lefttailed(vector<string>, vector<string>, const string&, double);
void hypothesistest_twoproportions_twotailed(vector<string>, vector<string>, const string&, double);
void hypothesistest_onesamplevariance_righttailed(vector<double>, double, double);
void hypothesistest_onesamplevariance_lefttailed(vector<double>, double, double);
void hypothesistest_onesamplevariance_twotailed(vector<double>, double, double);
void hypothesistest_twosamplesvariances_righttailed(vector<double>, vector<double>, double);
void hypothesistest_twosamplesvariances_lefttailed(vector<double>, vector<double>, double);
void hypothesistest_twosamplesvariances_twotailed(vector<double>, vector<double>, double);
void hypothesistest_categoricaldata(vector<vector<double>>, double);

void countstring(vector<string>, const string&);
void addstring(vector<string>, const string&);
void regressionline(vector<vector<double>>, double);
void ANOVA(vector<vector<double>>);
#endif
#endif


#endif
