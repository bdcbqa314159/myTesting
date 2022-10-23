#ifndef REAL_FUNCTION_2D
#define REAL_FUNCTION_2D

class RealFunction2D
{
public:
    virtual ~RealFunction2D() {}
    virtual double operator()(double x, double y) = 0;
};

#endif