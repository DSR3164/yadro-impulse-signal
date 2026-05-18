#include "dsp.hpp"

#include <cstdint>
#include <vector>
#include <ccomplex>
#include <algorithm>
#include <cmath>
#include <string>
#include <numbers>

void generate_sin(std::vector<double> &signal, double frequency, double phase, double sample_rate, double seconds)
{
    size_t samples_count = static_cast<size_t>(seconds * sample_rate);
    double t = 1.0 / sample_rate;
    signal.resize(samples_count);
    constexpr double TWO_PI = 2.0 * std::numbers::pi_v<double>; // гарантированно поддерживается на всех компиляторах в отличие от M_PI
    // Привожу фазу в диапазон [0, 2pi)
    if (std::abs(phase) >= TWO_PI)
        phase = std::fmod(phase, TWO_PI);

    for (size_t n = 0; n < samples_count; ++n)
        signal[n] = std::sin(TWO_PI * frequency * (t * static_cast<double>(n)) + phase);
}

std::vector<double> fir_lowpass_sinc_windowed(size_t L, double f_c)
{
    std::vector<double> h(L);
    double m = L - 1;

    constexpr double pi = std::numbers::pi_v<double>;

    for (size_t i = 0; i < L; ++i)
    {
        // Центрирую индекс относительно нуля
        double n = i - m / 2.0;
        if (std::abs(n) < 1e-9)
            h[i] = 2.0 * f_c;
        else
            h[i] = std::sin(2.0 * pi * f_c * n) / (pi * n);

        // Применяю окно Блэкмана-Харриса для подавления боковых лепестков
        double window = 0.35875 - 0.48829 * std::cos(2 * pi * i / m) + 0.14128 * std::cos(4 * pi * i / m) - 0.01168 * std::cos(6 * pi * i / m);
        h[i] *= window;
    }

    return h;
}

void filter(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &y)
{
    const size_t nb = b.size();
    const size_t na = a.size() + nb - 1;

    y.resize(na);
    std::fill(y.begin(), y.end(), 0.0);

    for (size_t n = 0; n < na; ++n)
    {
        double acc = 0.0;
        size_t limit = std::min(nb, n + 1);
        for (size_t m = 0; m < limit; ++m)
            acc += a[n - m] * b[m];
        y[n] = acc;
    }
}

std::vector<double> upsampling(std::vector<double> &signal, std::vector<double> &upsampled_signal)
{
    size_t upsample_coef = 2;
    int h_len = 21;
    double cutoff = 0.1;
    std::vector<double> temp(signal.size() * upsample_coef);
    size_t delay = h_len / 2 + 3;
    auto h = fir_lowpass_sinc_windowed(h_len, cutoff);

    save_signal("../../files/h.bin", h); // сохраняю импульсную характеристику для анадиза
    upsample(signal, temp, upsample_coef);
    filter(temp, h, upsampled_signal);

    // Убираю задержку фильтра и обрезаю что бы было четное кол-во периодов синуса
    // для того что бы убрать спектральные утечки на FFT для более легкого анализа работы фильтра
    upsampled_signal.erase(upsampled_signal.begin(), upsampled_signal.begin() + delay + 8);
    upsampled_signal.erase(upsampled_signal.end() - delay - 10, upsampled_signal.end());
    return temp;
}

void filter_int(const std::vector<int16_t> &a, const std::vector<double> &b, std::vector<int16_t> &y)
{
    const size_t nb = b.size();
    const size_t na = a.size() + nb - 1;

    y.resize(na);
    std::fill(y.begin(), y.end(), 0);

    for (size_t n = 0; n < na; ++n)
    {
        double acc = 0.0;
        size_t limit = std::min(nb, n + 1);
        for (size_t m = 0; m < limit; ++m)
            acc += static_cast<double>(a[n - m]) * b[m];

        acc = std::round(acc);
        if (acc > 32767.0)
            acc = 32767.0;
        else if (acc < -32768.0)
            acc = -32768.0;

        y[n] = static_cast<int16_t>(acc);
    }
}

std::vector<int16_t> upsampling_int(std::vector<int16_t> &signal, std::vector<int16_t> &upsampled_signal)
{
    size_t upsample_coef = 2;
    int h_len = 21;
    double cutoff = 0.1;
    std::vector<int16_t> temp(signal.size() * upsample_coef);
    size_t delay = h_len / 2 + 3;
    auto h = fir_lowpass_sinc_windowed(h_len, cutoff);
    save_signal("../../files/h.bin", h);
    upsample(signal, temp, upsample_coef);
    filter_int(temp, h, upsampled_signal);

    // Убираю задержку фильтра и обрезаю что бы было четное кол-во периодов синуса
    // для того что бы убрать спектральные утечки на FFT для более легкого анализа работы фильтра
    upsampled_signal.erase(upsampled_signal.begin(), upsampled_signal.begin() + delay + 8);
    upsampled_signal.erase(upsampled_signal.end() - delay - 10, upsampled_signal.end());
    return temp;
}