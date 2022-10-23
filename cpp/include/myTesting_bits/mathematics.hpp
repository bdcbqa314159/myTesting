#ifndef MATHEMATICS
#define MATHEMATICS

#include "standard.hpp"
#include "RealFunction.hpp"

double min(const std::vector<double> &v);
double max(const std::vector<double> &v);
double sum(const std::vector<double> &v);
double mean(const std::vector<double> &v);
double prctile(const std::vector<double> &v, double percentage);
double standardDeviation(const std::vector<double> &v, bool population = false);

std::vector<double> sort(const std::vector<double> &v);
std::vector<double> linspace(double from, double to, int numPoints);
std::vector<double> randUniform(int n);
std::vector<double> randN(int n);

void rng(const std::string &setting);

void plot(const std::string &fileName,
          const std::vector<double> &x,
          const std::vector<double> &y);

void hist(const std::string &fileName,
          const std::vector<double> &values,
          int numBuckets = 10);

double normcdf(double x);
double normInv(double x);
double integral(RealFunction &f, double a, double b, int nSteps);

double integralToInfinity(RealFunction &f, double x, int nSteps);

double integralFromInfinity(RealFunction &f, double x, int nSteps);

double integralOverR(RealFunction &f, int nSteps);

void testMathematics();

#endif