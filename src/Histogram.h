#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <string>

class Histogram
{
public:
    Histogram(double binWidth, double maxDiameter);
    void addVolume(double diameter, double volume);
    void saveCSV(const std::string &filename) const;
    void print() const;
private:
    double binWidth;
    double maxDiameter;
    size_t numBins;
    std::vector<double> volumeBins;
};

#endif