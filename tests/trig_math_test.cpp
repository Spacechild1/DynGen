#include "math_utils.h"

#include <iostream>
#include <random>

double radToDeg(double rad) {
    return rad / M_PI * 180.0;
}

int main() {
    std::mt19937 randEngine(0);

#if LUT_INTERPOLATION
    std::cout << "with LUT interpolation" << std::endl;
#else
    std::cout << "without LUT interpolation" << std::endl;
#endif

    const double maxMagValue = 512.0;
#if LUT_INTERPOLATION
    const double maxAllowedPhaseErrorDegrees = 0.0005;
    const double maxAllowedMagErrorNormalized = 0.00001;
#else
    const double maxAllowedPhaseErrorDegrees = 0.1;
    const double maxAllowedMagErrorNormalized = 0.001;
#endif
    const size_t numIterations = 10000;

    // test sin/cos over the range of -4pi, 4pi
    std::uniform_real_distribution<double> randPhaseDist(M_PI * -4.0, M_PI * 4.0);
    std::uniform_real_distribution<double> randMagDist(-512.0, 512.0);

    auto randPhase = [&]() { return randPhaseDist(randEngine); };
    auto randMag = [&]() { return randMagDist(randEngine); };

    // test sine
    {
        double minError = std::numeric_limits<double>::max();
        double maxError = 0.0;
        double avgError = 0.0;
        for (size_t i = 0; i < numIterations; ++i) {
            double phase = randPhase();
            double a = std::sin(phase);
            double b = detail::fastSin(phase);
            double err = std::abs(a - b);
#if 0
            std::cout << "std::sin: " << a << ", fastSin: " << b
                      << ", err: " << err << std::endl;
#endif
            if (err > maxError) {
                maxError = err;
            }
            if (err < minError) {
                minError = err;
            }
            avgError += err;
        }
        avgError /= numIterations;

        std::cout << "sin() min. error: " << radToDeg(minError) << "°"
                  << ", max. error: " << radToDeg(maxError) << "°"
                  << ", avg. error: " << radToDeg(avgError) << "°"
                  << std::endl;

        if (radToDeg(maxError) > maxAllowedPhaseErrorDegrees) {
            return 1;
        }
    }

    // test cosine
    {
        double minError = std::numeric_limits<double>::max();
        double maxError = 0.0;
        double avgError = 0.0;
        for (size_t i = 0; i < numIterations; ++i) {
            double phase = randPhase();
            double a = std::cos(phase);
            double b = detail::fastCos(phase);
            double err = std::abs(a - b);
#if 0
            std::cout << "std::cos: " << a << ", fastCos: " << b
                      << ", err: " << err << std::endl;
#endif
            if (err > maxError) {
                maxError = err;
            }
            if (err < minError) {
                minError = err;
            }
            avgError += err;
        }
        avgError /= numIterations;

        std::cout << "cos() min. error: " << radToDeg(minError) << "°"
                  << ", max. error: " << radToDeg(maxError) << "°"
                  << ", avg. error: " << radToDeg(avgError) << "°"
                  << std::endl;

        if (radToDeg(maxError) > maxAllowedPhaseErrorDegrees) {
            return 1;
        }
    }

    // test complex to polar
    {
        double minMagError = std::numeric_limits<double>::max();
        double maxMagError = 0.0;
        double avgMagError = 0.0;
        double minPhaseError = std::numeric_limits<double>::max();
        double maxPhaseError = 0.0;
        double avgPhaseError = 0.0;
        for (size_t i = 0; i < numIterations; ++i) {
            double real = randMag();
            double imag = randMag();

            double mag1 = std::sqrt(real * real + imag * imag);
            double phase1 = std::atan2(imag, real);

            auto [mag2, phase2] = complexToPolar(real, imag);

            auto magErr = std::abs(mag1 - mag2) / maxMagValue;
            auto phaseErr = std::abs(phase1 - phase2);
            if (phaseErr >= M_PI) {
                phaseErr = std::abs(phaseErr - M_PI * 2.0);
            }

            if (magErr > maxMagError) {
                maxMagError = magErr;
            }
            if (magErr < minMagError) {
                minMagError = magErr;
            }
            avgMagError += magErr;

            if (phaseErr > maxPhaseError) {
                maxPhaseError = phaseErr;
            }
            if (phaseErr < minPhaseError) {
                minPhaseError = phaseErr;
            }
            avgPhaseError += phaseErr;
        }

        avgMagError /= numIterations;
        avgPhaseError /= numIterations;

        std::cout << "complexToPolar(): min. mag error: " << minMagError
                  << ", max. mag error: " << maxMagError << ", avg. mag error: "
                  << avgMagError << ",\nmin. phase error: " << radToDeg(minPhaseError) << "°"
                  << ", max. phase error: " << radToDeg(maxPhaseError) << "°"
                  << ", avg. phase error: " << radToDeg(avgPhaseError) << "°"
                  << std::endl;

        if (maxMagError > maxAllowedMagErrorNormalized) {
            return 1;
        }

        if (radToDeg(maxPhaseError) > maxAllowedPhaseErrorDegrees) {
            return 1;
        }
    }

    // test polar to complex
    {
        double minMagError = std::numeric_limits<double>::max();
        double maxMagError = 0.0;
        double avgMagError = 0.0;
        for (size_t i = 0; i < numIterations; ++i) {
            double mag = randMag();
            double phase = randPhase();

            double real1 = mag * std::cos(phase);
            double imag1 = mag * std::sin(phase);

            auto [real2, imag2] = polarToComplex(mag, phase);

            auto realErr = std::abs(real1 - real2) / maxMagValue;
            auto imagErr = std::abs(imag1 - imag2) / maxMagValue;
            auto magErr = std::max(realErr, imagErr);

            if (magErr > maxMagError) {
                maxMagError = magErr;
            }
            if (magErr < minMagError) {
                minMagError = magErr;
            }
            avgMagError += magErr;
        }

        avgMagError /= numIterations;

        std::cout << "polarToComplex(): min. mag error: " << minMagError
                  << ", max. mag error: " << maxMagError << ", avg. mag error: "
                  << avgMagError << std::endl;

        if (maxMagError > maxAllowedMagErrorNormalized) {
            return 1;
        }
    }

    return 0;
}
