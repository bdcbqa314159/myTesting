#ifndef LINE_CHART
#define LINE_CHART

#include "standard.hpp"

class LineChart
{
public:
    LineChart();
    void setTitle(const std::string &title);
    void setSeries(const std::vector<double> &x, const std::vector<double> &y);
    void writeAsHTML(std::ostream &out) const;
    void writeAsHTML(const std::string &file) const;

private:
    std::string title;
    std::vector<double> x;
    std::vector<double> y;
};

void testLineChart();

#endif