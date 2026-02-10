#include "math_utils.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <vector>

int main() {
    std::array bufferSizes = {
        16, 32, 64, 128, 256, 512, 1024,
        2048, 4096, 8192, 16384, 32768
    };

    std::vector<double> buffer;

    const double maxAllowedErrorNormalized = 0.00001;

    // test Hann window with several window sizes
    for (const auto bufferSize : bufferSizes) {
        double minError = std::numeric_limits<double>::max();
        double maxError = 0.0;
        double avgError = 0.0;

        buffer.resize(bufferSize);
        std::fill(buffer.begin(), buffer.end(), 1.0);

        applyHannWindow<16>(buffer.data(), buffer.size());

        const double toPhase = M_PI * 2.0 / bufferSize;

        for (int i = 0; i < bufferSize; ++i) {
            double a = buffer[i];
            double b = std::cos(i * toPhase) * -0.5 + 0.5;
            double err = std::abs(a - b);
#if 0
            std::cout << "LUT: " << a << ", std::cos(): " << b << ", err: " << err << std::endl;
#endif
            if (err > maxError) {
                maxError = err;
            }
            if (err < minError) {
                minError = err;
            }
            avgError += err;
        }

        avgError /= bufferSize;

        std::cout << "Apply Hann window to " << bufferSize << " samples."
                  << " Min. error: " << minError << ", max. error: " << maxError
                  << ", avg. error: " << avgError << std::endl;

        if (maxError > maxAllowedErrorNormalized) {
            return 1;
        }
    }

    // test sine window with several window sizes
    for (const auto bufferSize : bufferSizes) {
        double minError = std::numeric_limits<double>::max();
        double maxError = 0.0;
        double avgError = 0.0;

        buffer.resize(bufferSize);
        std::fill(buffer.begin(), buffer.end(), 1.0);

        applySineWindow<16>(buffer.data(), buffer.size());

        const double toPhase = M_PI / bufferSize;

        for (int i = 0; i < bufferSize; ++i) {
            double a = buffer[i];
            double b = std::sin(i * toPhase);
            double err = std::abs(a - b);
#if 0
            std::cout << "LUT: " << a << ", std::sin(): " << b << ", err: " << err << std::endl;
#endif
            if (err > maxError) {
                maxError = err;
            }
            if (err < minError) {
                minError = err;
            }
            avgError += err;
        }

        avgError /= bufferSize;

        std::cout << "Apply sine window to " << bufferSize << " samples."
                  << " Min. error: " << minError << ", max. error: " << maxError
                  << ", avg. error: " << avgError << std::endl;

        if (maxError > maxAllowedErrorNormalized) {
            return 1;
        }
    }

    return 0;
}
