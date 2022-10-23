#ifndef REAL_FUNCTION
#define REAL_FUNCTION

class RealFunction
{
public:
    virtual ~RealFunction() {}
    virtual double operator()(double x) = 0;
};

#endif