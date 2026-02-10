#pragma once

#include <cmath>
#include <utility>

#ifndef M_PI
#    define M_PI    3.14159265358979323846
#    define M_PI_2  1.57079632679489661923
#endif
#define TWO_PI (M_PI * 2.0)

// if LUT_INTERPOLATION is 1 we use smaller lookup tables with linear interpolation.
// This increases accuracy by about 100x while only using 1/8 of the memory.
// The smaller tables easily fit into L1 cache, which can make up for the extra
// instructions that are needed for the linear interpolation.
#ifndef LUT_INTERPOLATION
#    define LUT_INTERPOLATION 1
#endif

namespace detail {

#if LUT_INTERPOLATION
constexpr size_t kSineTableSize = 1024;
constexpr size_t kSineTableRealSize = kSineTableSize + kSineTableSize / 4;
constexpr size_t kPolarTableSize = 256;
#else
constexpr size_t kSineTableSize = 8192;
constexpr size_t kSineTableRealSize = kSineTableSize;
constexpr size_t kPolarTableSize = 2048;
#endif
constexpr size_t kHannTableSize = 1024;

constexpr size_t kSineTableMask = kSineTableSize - 1;
constexpr double kSinePhaseToIndex = static_cast<double>(kSineTableSize) / TWO_PI;
constexpr size_t kPolarTableSize2 = kPolarTableSize / 2;

// QUESTION: should we use floats to reduce memory/cache usage?
// NOTE: the sine table has additional pi/2 values so it also works as a cosine table.
inline double gSineTable[kSineTableRealSize + 1];
inline double gMagTable[kPolarTableSize + 1];
inline double gPhaseTable[kPolarTableSize + 1];
inline double gHannTable[kHannTableSize + 1];

// NOTE: 'findex' is supposed to be in the range of the table
inline double readTableLin(const double* tab, double findex) {
    size_t index = static_cast<size_t>(findex);
    double fract = findex - static_cast<double>(index);
    double a = tab[index];
    double b = tab[index + 1];
    return a + (b - a) * fract;
}

#if LUT_INTERPOLATION

inline double fastCos(double rad) {
    // cos() is symmetric around 0
    double findex = rad >= 0.f ? rad * kSinePhaseToIndex : -rad * kSinePhaseToIndex;
    // first find the fractional part
    size_t index = static_cast<size_t>(findex);
    double fract = findex - static_cast<double>(index);
    // then clamp the index
    index &= kSineTableMask;
    // start reading at pi/2 to get a cosine wave
    const auto* tab = gSineTable + kSineTableSize / 4;
    double a = tab[index];
    double b = tab[index + 1];
    return a + (b - a) * fract;
}

inline double fastSin(double rad) {
    return fastCos(rad - M_PI_2);
}

#else

inline double fastSin(double rad) {
    size_t index = static_cast<size_t>(rad * kSinePhaseToIndex) & kSineTableMask;
    return gSineTable[index];
}

inline double fastCos(double rad) {
#if 1
    size_t index = static_cast<size_t>(rad * kSinePhaseToIndex);
    index = (index + kSineTableSize / 4) & kSineTableMask;
#else
    size_t index = static_cast<size_t>(rad * kSinePhaseToIndex) & kSineTableMask;
    index = (index + (kSineTableSize / 4)) & kSineTableMask;
#endif
    return gSineTable[index];
}

#endif

inline bool lutInitialized = []() {
    using namespace detail;

    for (int i = 0; i <= kSineTableRealSize; ++i) {
        double phase = i * TWO_PI / static_cast<double>(kSineTableSize);
        gSineTable[i] = std::sin(phase);
    }

    for (int i = 0; i <= kPolarTableSize; ++i) {
        // compute atan over the range of -1, 1
        double offset = static_cast<double>(kPolarTableSize2);
        double slope = (static_cast<double>(i) - offset) / offset;
        double angle = std::atan(slope);
        gPhaseTable[i] = angle;
        // cos(α) = real / mag -> mag / real = 1 / cos(α)
        // since we have tan(α) (= imag / real) we can compute mag with gMagTable[index] * real.
        gMagTable[i] = 1.0 / std::cos(angle);
    }

    for (int i = 0; i <= kHannTableSize; ++i) {
        double phase = i * TWO_PI / static_cast<double>(kSineTableSize);
        gHannTable[i] = 0.5 - std::cos(phase) * 0.5;
    }

    return true;
}();

} // namespace detail

/*! @brief takes a complex number and returns the magnitude and phase */
inline std::pair<double, double> complexToPolar(double real, double imag) {
    using namespace detail;

    double absreal = std::abs(real);
    double absimag = std::abs(imag);
    double mag, phase;
    // If abs(y) > abs(x) the slope would be larger than 1.0; in this case
    // we can swap x and y using the the following identity:
    // atan(y/x) = -atan(x/y) + sng(y)*pi/2
    // This makes sure that the slope stays in the interval -1.0, 1.0 and
    // we can do a simple table lookup.
    if (absreal > absimag) {
        double slope = imag / real;
#if LUT_INTERPOLATION
        double findex = slope * kPolarTableSize2 + kPolarTableSize2;
        mag = readTableLin(gMagTable, findex) * absreal;
        phase = readTableLin(gPhaseTable, findex);
#else
        size_t index = static_cast<size_t>(slope * kPolarTableSize2) + kPolarTableSize2;
        mag = gMagTable[index] * absreal;
        phase = gPhaseTable[index];
#endif
        if (real < 0) {
            phase += M_PI;
        }
    } else if (absimag > 0) {
        double slope = real / imag;
#if LUT_INTERPOLATION
        double findex = slope * kPolarTableSize2 + kPolarTableSize2;
        mag = readTableLin(gMagTable, findex) * absimag;
        phase = readTableLin(gPhaseTable, findex);
#else
        size_t index = static_cast<size_t>(slope * kPolarTableSize2) + kPolarTableSize2;
        mag = gMagTable[index] * absimag;
        phase = gPhaseTable[index];
#endif
        if (imag > 0) {
            phase = M_PI_2 - phase;
        } else {
            phase = M_PI + M_PI_2 - phase;
        }
    } else {
        mag = 0;
        phase = 0;
    }
    return { mag, phase };
}

/*! @brief takes magnitude and phase and returns a complex number */
inline std::pair<double, double> polarToComplex(double mag, double phase) {
    // fast approximation
    double real = mag * detail::fastCos(phase);
    double imag = mag * detail::fastSin(phase);
    return { real, imag };
}
