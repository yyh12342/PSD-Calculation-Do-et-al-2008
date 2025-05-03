#include "Histogram.h"

#include <cmath>
#include <fstream>
#include <iostream>

Histogram::Histogram(double binWidth_, double maxDiameter_)
    : binWidth(binWidth_), maxDiameter(maxDiameter_)
{
    if (binWidth_ <= 0 || maxDiameter_ <= 0)
    {
        numBins = 0;
    }
    else
    {
        numBins = static_cast<size_t>(std::ceil(maxDiameter_ / binWidth_));
        volumeBins.assign(numBins, 0.0);
    }
}

// 각 지름이 속한 bin에 volume 누적
void Histogram::addVolume(double diameter, double volume)
{
    if (numBins == 0)
    {
        return;
    }

    if (diameter < 0)
    {
        diameter = 0;
    }

    if (diameter > maxDiameter)
    {
        diameter = maxDiameter;
    }

    size_t index = static_cast<size_t>(diameter / binWidth);

    if (index >= numBins)
    {
        index = numBins - 1;
    }

    volumeBins[index] += volume;
}

// csv로 저장
void Histogram::saveCSV(const std::string &filename) const
{
    std::ofstream ofs(filename);

    if (!ofs)
    {
        std::cerr << "Error: Could not open output file " << filename << std::endl;

        return;
    }

    ofs << "Diameter(nm),Volume(nm^3)\n";

    for (size_t i = 0; i < numBins; ++i)
    {
        double diam_center = (i + 0.5) * binWidth;

        if (diam_center > maxDiameter)
        {
            diam_center = maxDiameter;
        }

        ofs << diam_center << "," << volumeBins[i] << "\n";
    }

    ofs.close();
}

// 따로 콘솔에 출력
void Histogram::print() const
{
    std::cout << "Diameter(nm),Volume(nm^3)\n";

    for (size_t i = 0; i < numBins; ++i)
    {
        double diam_center = (i + 0.5) * binWidth;

        if (diam_center > maxDiameter)
        {
            diam_center = maxDiameter;
        }

        std::cout << diam_center << "," << volumeBins[i] << "\n";
    }
}