#include "../include/myTesting_bits/mathematics.hpp"
#include "../include/myTesting_bits/LineChart.hpp"
#include "../include/myTesting_bits/Histogram.hpp"
#include "../include/myTesting_bits/standard.hpp"
#include "../include/myTesting_bits/globalConstants.hpp"
#include "../include/myTesting_bits/testing.hpp"

using namespace std;

vector<double> linspace(double from, double to, int numPoints)
{
    ASSERT(numPoints >= 2);
    vector<double> ret(numPoints, 0.);
    double step = (to - from) / (numPoints - 1);

    double current = from;

    for (int i = 0; i < numPoints; i++)
    {
        ret[i] = current;
        current += step;
    }
    return ret;
}

double min(const vector<double> &v)
{
    vector<double> copy = v;
    sort(copy.begin(), copy.end());
    return copy.at(0);
}

double max(const vector<double> &v)
{
    vector<double> copy = v;
    sort(copy.begin(), copy.end());
    return copy.back();
}

vector<double> sort(const vector<double> &v)
{
    vector<double> copy = v;
    sort(copy.begin(), copy.end());
    return copy;
}

double sum(const vector<double> &v)
{
    double total = 0.0;
    int n = v.size();
    for (int i = 0; i < n; i++)
    {
        total += v[i];
    }
    return total;
}

double mean(const vector<double> &v)
{
    int n = v.size();
    ASSERT(n > 0);
    return sum(v) / n;
}

double standardDeviation(const vector<double> &v, bool population)
{
    int n = v.size();
    double total = 0.0;
    double totalSq = 0.0;
    for (int i = 0; i < n; i++)
    {
        total += v[i];
        totalSq += v[i] * v[i];
    }
    if (population)
    {
        ASSERT(n > 0);
        return sqrt((totalSq - total * total / n) / n);
    }
    else
    {
        ASSERT(n > 1);
        return sqrt((totalSq - total * total / n) / (n - 1));
    }
}

double prctile(const vector<double> &v, double percentage)
{
    // See the text for a precise specification

    ASSERT(percentage >= 0.0);
    ASSERT(percentage <= 100.0);
    int n = v.size();
    vector<double> sorted = sort(v);

    int indexBelow = static_cast<int>(n * percentage / 100.0 - 0.5);
    int indexAbove = indexBelow + 1;
    if (indexAbove > n - 1)
    {
        return sorted[n - 1];
    }
    if (indexBelow < 0)
    {
        return sorted[0];
    }
    double valueBelow = sorted[indexBelow];
    double valueAbove = sorted[indexAbove];
    double percentageBelow = 100.0 * (indexBelow + 0.5) / n;
    double percentageAbove = 100.0 * (indexAbove + 0.5) / n;
    if (percentage <= percentageBelow)
    {
        return valueBelow;
    }
    if (percentage >= percentageAbove)
    {
        return valueAbove;
    }
    double correction = (percentage - percentageBelow) * (valueAbove - valueBelow) / (percentageAbove - percentageBelow);
    return valueBelow + correction;
}

static mt19937 mersenneTwister;

vector<double> randUniform(int n)
{
    vector<double> result(n);

    for (int i = 0; i < n; i++)
    {
        result[i] = (mersenneTwister() + 0.5) / (mersenneTwister.max() + 1.0);
    }

    return result;
}

void rng(const string &description)
{
    ASSERT(description == "default");
    mersenneTwister.seed(mt19937::default_seed);
}

vector<double> randN(int n)
{

    vector<double> result = randUniform(n);

    for (int i = 0; i < n; i++)
    {
        result[i] = normInv(result[i]);
    }

    return result;
}

void plot(const string &file,
          const vector<double> &x,
          const vector<double> &y)
{
    LineChart lc;
    lc.setSeries(x, y);
    lc.writeAsHTML(file);
}

/**
 *  Convenience method for generating plots
 */
void hist(const string &file,
          const vector<double> &data,
          int numBuckets)
{
    Histogram h;
    h.setData(data);
    h.setNumBuckets(numBuckets);
    h.writeAsHTML(file);
}

static const double ROOT_2_PI = sqrt(2. * PI);

static inline double hornerFunction(double x, double a0)
{
    return a0;
}

static inline double hornerFunction(double x, double a0, double a1)
{
    return a0 + x * hornerFunction(x, a1);
}

static inline double hornerFunction(double x, double a0, double a1, double a2)
{
    return a0 + x * hornerFunction(x, a1, a2);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3)
{
    return a0 + x * hornerFunction(x, a1, a2, a3);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7, a8);
}

double normcdf(double x)
{
    if (x < 0)
    {
        return 1. - normcdf(-x);
    }

    const double a0 = 0.;
    const double a1 = 0.319381530;
    const double a2 = -0.356563782;
    const double a3 = 1.781477937;
    const double a4 = -1.821255978;
    const double a5 = 1.330274429;

    double k = 1. / (1. + 0.2316419 * x);
    double poly = hornerFunction(k, a0, a1, a2, a3, a4, a5);

    double approx = 1. - 1. / ROOT_2_PI * exp(-0.5 * x * x) * poly;

    return approx;
}

double normInv(double x)
{
    const double a0 = 2.50662823884;
    const double a1 = -18.61500062529;
    const double a2 = 41.39119773534;
    const double a3 = -25.44106049637;
    const double b1 = -8.47351093090;
    const double b2 = 23.08336743743;
    const double b3 = -21.06224101826;
    const double b4 = 3.13082909833;
    const double c0 = 0.3374754822726147;
    const double c1 = 0.9761690190917186;
    const double c2 = 0.1607979714918209;
    const double c3 = 0.0276438810333863;
    const double c4 = 0.0038405729373609;
    const double c5 = 0.0003951896511919;
    const double c6 = 0.0000321767881768;
    const double c7 = 0.0000002888167364;
    const double c8 = 0.0000003960315187;

    double y = x - 0.5;
    double r{};
    if (y < 0.42 && y > -0.42)
    {
        r = y * y;
        return y * hornerFunction(r, a0, a1, a2, a3) / hornerFunction(r, 1., b1, b2, b3, b4);
    }
    else
    {
        if (y < 0.)
        {
            r = x;
        }
        else
        {
            r = 1. - x;
        }

        double s = log(-log(r));
        double t = hornerFunction(s, c0, c1, c2, c3, c4, c5, c6, c7, c8);

        if (x > 0.5)
            return t;
        else
            return -t;
    }
}

double integral(RealFunction &f, double a, double b, int nPoints)
{
    double h = (b - a) / nPoints;
    double x = a + 0.5 * h;
    double total = 0.;
    for (int i = 0; i < nPoints; i++)
    {
        double y = f(x);
        total += y;
        x += h;
    }
    return h * total;
}

double integralToInfinity(RealFunction &f, double x, int nPoints)
{

    class Integrand : public RealFunction
    {
    public:
        double x;
        RealFunction &f;

        virtual double operator()(double s) override
        {
            return 1 / (s * s) * f(1 / s + x - 1);
        }

        Integrand(double x, RealFunction &f) : x{x}, f{f}
        {
        }
    };

    Integrand i(x, f);
    return integral(i, 0, 1, nPoints);
}

double integralFromInfinity(RealFunction &f, double x, int nPoints)
{

    class ReversedIntegrand : public RealFunction
    {
    public:
        RealFunction &f;
        virtual double operator()(double x) override
        {
            return f(-x);
        }
        ReversedIntegrand(RealFunction &f) : f{f} {}
    };

    ReversedIntegrand i(f);
    return integralToInfinity(i, -x, nPoints);
}

double integralOverR(RealFunction &f, int nPoints)
{
    return integralToInfinity(f, 0, nPoints) + integralToInfinity(f, 0, nPoints);
}

static vector<double> createTestVector()
{
    vector<double> v = {1, 5, 3, 9, 7};
    return v;
}

static void testLinspace()
{
    vector<double> result = linspace(1.0, 10.0, 4);
    ASSERT_APPROX_EQUAL(result[0], 1.0, 0.001);
    ASSERT_APPROX_EQUAL(result[1], 4.0, 0.001);
    ASSERT_APPROX_EQUAL(result[2], 7.0, 0.001);
    ASSERT_APPROX_EQUAL(result[3], 10.0, 0.001);
}

static void testMin()
{
    ASSERT_APPROX_EQUAL(min(createTestVector()), 1.0, 0.001);
}

static void testMax()
{
    vector<double> v = createTestVector();
    sort(v.begin(), v.end());
    cout << v.back() << endl;
    cout << max(createTestVector()) << endl;
    ASSERT_APPROX_EQUAL(max(createTestVector()), 9.0, 0.001);
}

static void testMean()
{
    ASSERT_APPROX_EQUAL(mean(createTestVector()), 5.0, 0.001);
}

static void testStandardDeviation()
{
    ASSERT_APPROX_EQUAL(standardDeviation(createTestVector()), 3.1623, 0.001);
    ASSERT_APPROX_EQUAL(standardDeviation(createTestVector(), true), 2.8284, 0.001);
}

static void testRandUniform()
{
    rng("default");
    vector<double> v = randUniform(1000);
    ASSERT(((int)v.size()) == 1000);
    ASSERT_APPROX_EQUAL(mean(v), 0.5, 0.1);
    ASSERT(max(v) < 1.0);
    ASSERT(min(v) > 0.0);
}

static void testRandn()
{
    rng("default");
    vector<double> v = randN(10000);
    ASSERT(((int)v.size()) == 10000);
    ASSERT_APPROX_EQUAL(mean(v), 0.0, 0.1);
    ASSERT_APPROX_EQUAL(standardDeviation(v), 1.0, 0.1);
}

static void testPrctile()
{
    const vector<double> v = createTestVector();
    ASSERT_APPROX_EQUAL(prctile(v, 100.0), 9.0, 0.001);
    ASSERT_APPROX_EQUAL(prctile(v, 0.0), 1.0, 0.001);
    ASSERT_APPROX_EQUAL(prctile(v, 50.0), 5.0, 0.001);
    ASSERT_APPROX_EQUAL(prctile(v, 17.0), 1.7, 0.001);
    ASSERT_APPROX_EQUAL(prctile(v, 62.0), 6.2, 0.001);
}

static void testNormCdf()
{
    ASSERT(normcdf(0.3) > 0);
    ASSERT(normcdf(0.3) < 1);
    ASSERT_APPROX_EQUAL(normcdf(-1e10), 0., 0.001);
    ASSERT_APPROX_EQUAL(normcdf(1e10), 1., 0.001);
    ASSERT(normcdf(0.3) < normcdf(0.5));
    ASSERT_APPROX_EQUAL(normcdf(0.3), 1 - normcdf(-0.3), 0.01);
    ASSERT_APPROX_EQUAL(normcdf(0.0), 0.5, 0.0001);

    ASSERT_APPROX_EQUAL(normcdf(normInv(0.3)),
                        0.3, 0.0001);

    ASSERT_APPROX_EQUAL(normcdf(1.96), 0.975, 0.001);
}

static void testNormInv()
{
    ASSERT_APPROX_EQUAL(normInv(0.975), 1.96, 0.01);
}

/**
 *  When you create a small class like this, using
 *  nested classes is easier.
 */
static void testIntegral()
{

    class Sin : public RealFunction
    {
        virtual double operator()(double x) override
        {
            return sin(x);
        }
    };

    Sin integrand;
    double actual = integral(integrand, 1, 3, 1000);
    double expected = -cos(3.0) + cos(1.0);
    ASSERT_APPROX_EQUAL(actual, expected, 0.000001);
}

/*  PDF of the normal distribution */
class NormalPDF : public RealFunction
{
public:
    virtual double operator()(double x) override
    {
        return exp(-0.5 * x * x) / ROOT_2_PI;
    }
};

static void testIntegrateNormal()
{
    NormalPDF pdf;
    double result = integral(pdf, -1.96, 1.96, 10000);
    ASSERT_APPROX_EQUAL(result, 0.95, 0.001);
}

static void testIntegrateToInfinity()
{
    NormalPDF pdf;
    double x = -1.96;
    double result = integralToInfinity(pdf, x, 10000);
    ASSERT_APPROX_EQUAL(result, normcdf(-x), 0.001);
    result = integralFromInfinity(pdf, x, 10000);
    ASSERT_APPROX_EQUAL(result, normcdf(x), 0.001);
    result = integralOverR(pdf, 10000);
    ASSERT_APPROX_EQUAL(result, 1.0, 0.001);
}

void testMathematics()
{
    TEST(testMean);
    TEST(testStandardDeviation);
    TEST(testMin);
    TEST(testMax);
    TEST(testRandUniform);
    TEST(testPrctile);
    TEST(testNormInv);
    TEST(testNormCdf);
    TEST(testIntegral);
    TEST(testIntegrateNormal);
    TEST(testIntegrateToInfinity);
}