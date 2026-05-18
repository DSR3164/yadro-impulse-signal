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

static constexpr int SCALE_BITS = 15;
static constexpr int32_t SCALE = (1 << SCALE_BITS);

std::vector<int16_t> fir_lowpass_sinc_windowed_int(int h_len, double cutoff)
{
    // Получаем вещественные коэффициенты из существующей функции
    auto h_double = fir_lowpass_sinc_windowed(h_len, cutoff);

    std::vector<int16_t> h_int(h_len);
    for (int i = 0; i < h_len; ++i)
    {
        double val = h_double[i] * SCALE;
        val = std::round(val);
        if (val > 32767.0)
            val = 32767.0;
        if (val < -32768.0)
            val = -32768.0;
        h_int[i] = static_cast<int16_t>(val);
    }
    return h_int;
}

void filter(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &y)
{
    const size_t a_size = a.size();
    const size_t nb = b.size();
    const size_t na = a.size() + nb - 1;

    y.resize(na);
    std::fill(y.begin(), y.end(), 0.0);

    for (size_t n = 0; n < na; ++n)
    {
        double acc = 0.0;
        size_t limit = std::min(nb, n + 1);
        for (size_t m = 0; m < limit; ++m)
            if (n - m < a_size)
                acc += a[n - m] * b[m];
        y[n] = acc;
    }
}

std::vector<double> upsampling(std::vector<double> &signal, std::vector<double> &upsampled_signal, double target_f, double sample_rate)
{
    size_t upsample_coef = 2;
    int h_len = 21;
    double cutoff = target_f / sample_rate / upsample_coef;
    std::vector<double> temp(signal.size() * upsample_coef);
    size_t delay = (h_len - 1) / 2;
    auto h = fir_lowpass_sinc_windowed(h_len, cutoff);

    save_signal("../../files/h.bin", h); // сохраняю импульсную характеристику для анадиза
    upsample(signal, temp, upsample_coef);
    filter(temp, h, upsampled_signal);

    // Убираю задержку фильтра и обрезаю что бы было четное кол-во периодов синуса
    // для того что бы убрать спектральные утечки на FFT для более легкого анализа работы фильтра
    upsampled_signal.erase(upsampled_signal.begin(), upsampled_signal.begin() + delay);
    upsampled_signal.erase(upsampled_signal.end() - delay, upsampled_signal.end());
    return temp;
}

void filter_int(const std::vector<int16_t> &a, const std::vector<int16_t> &b, std::vector<int16_t> &y)
{
    const size_t nb = b.size();
    const size_t na = a.size() + nb - 1;

    y.resize(na);
    std::fill(y.begin(), y.end(), 0);

    for (size_t n = 0; n < na; ++n)
    {
        int64_t acc = 0;
        size_t limit = std::min(nb, n + 1);

        for (size_t m = 0; m < limit; ++m)
            if (n - m < na - nb + 1)
                acc += static_cast<int64_t>(a[n - m]) * b[m];

        // Сдвиг вместо деления (компилятор всё равно так сделает, но явно)
        acc = (acc + (SCALE >> 1)) >> SCALE_BITS; // округление к ближайшему

        if (acc > 32767)
            acc = 32767;
        if (acc < -32768)
            acc = -32768;

        y[n] = static_cast<int16_t>(acc);
    }
}

std::vector<int16_t> upsampling_int(std::vector<int16_t> &signal, std::vector<int16_t> &upsampled_signal, double target_f, double sample_rate)
{
    size_t upsample_coef = 2;
    int h_len = 21;
    double cutoff = target_f / sample_rate / upsample_coef;
    std::vector<int16_t> temp(signal.size() * upsample_coef);
    size_t delay = (h_len - 1) / 2;
    auto h = fir_lowpass_sinc_windowed_int(h_len, cutoff);
    save_signal("../../files/h.bin", h);

    upsample(signal, temp, upsample_coef);
    filter_int(temp, h, upsampled_signal);

    upsampled_signal.erase(upsampled_signal.begin(), upsampled_signal.begin() + delay);
    upsampled_signal.erase(upsampled_signal.end() - delay, upsampled_signal.end());
    return temp;
}

std::vector<double> perform_float_test()
{
    size_t max_f = 50;
    double power_acc = 0.0;
    std::vector<double> buffer, sin_ref, sin;
    std::vector<double> errors_float(max_f);
    for (size_t i = 1; i <= max_f; ++i)
    {
        power_acc = 0.0;
        generate_sin(buffer, static_cast<double>(i), 0.0, 100, 10);
        generate_sin(sin_ref, static_cast<double>(i), 0.0, 200, 10);
        upsampling(buffer, sin, static_cast<double>(i), 100);

        size_t samples_per_period = static_cast<size_t>(100 / i);

        cutoff(sin, samples_per_period);
        cutoff(sin_ref, samples_per_period);

        auto [min_it1, max_it1] = std::minmax_element(sin.begin(), sin.end());
        double amplitude = std::max(std::abs(*min_it1), std::abs(*max_it1));
        auto [min_it_ref, max_it_ref] = std::minmax_element(sin_ref.begin(), sin_ref.end());
        double amplitude_ref = std::max(std::abs(*min_it_ref), std::abs(*max_it_ref));

        for (auto &x : sin)
            x /= amplitude;
        for (auto &x : sin_ref)
            x /= amplitude_ref;

        size_t min = std::min(sin.size(), sin_ref.size());

        for (size_t n = 0; n < min; ++n)
            power_acc += std::pow(sin_ref[n], 2.0);

        for (size_t n = 0; n < min; ++n)
            errors_float[i] += std::pow((sin_ref[n] - sin[n]), 2.0);

        power_acc /= static_cast<double>(sin.size());
        errors_float[i] /= static_cast<double>(sin.size());

        errors_float[i] = 10 * std::log10(power_acc / errors_float[i]);
    }
    return errors_float;
}

std::vector<double> perform_int_test()
{
    size_t max_f = 50;
    double power_acc = 0.0;
    std::vector<double> buffer, sin_ref, sin;
    std::vector<int16_t> buffer_clamped, sin_ref_clamped, sin_clamped;
    std::vector<double> errors_int(max_f);
    for (size_t i = 1; i <= max_f; ++i)
    {
        power_acc = 0.0;
        generate_sin(buffer, static_cast<double>(i), 0.0, 100, 10);
        generate_sin(sin_ref, static_cast<double>(i), 0.0, 200, 10);

        clamp_to(buffer, buffer_clamped);
        clamp_to(sin_ref, sin_ref_clamped);
        upsampling_int(buffer_clamped, sin_clamped, static_cast<double>(i), 100);

        size_t samples_per_period = static_cast<size_t>(100 / i);

        cutoff(sin_clamped, samples_per_period);
        cutoff(sin_ref_clamped, samples_per_period);

        sin.resize(sin_clamped.size());
        sin_ref.resize(sin_ref_clamped.size());

        for (size_t k = 0; k < sin.size(); ++k)
            sin[k] = static_cast<double>(sin_clamped[k]);

        for (size_t k = 0; k < sin_ref.size(); ++k)
            sin_ref[k] = static_cast<double>(sin_ref_clamped[k]);

        auto [min_it1, max_it1] = std::minmax_element(sin.begin(), sin.end());
        double amplitude = std::max(std::abs(*min_it1), std::abs(*max_it1));
        auto [min_it_ref, max_it_ref] = std::minmax_element(sin_ref.begin(), sin_ref.end());
        double amplitude_ref = std::max(std::abs(*min_it_ref), std::abs(*max_it_ref));

        for (auto &x : sin)
            x /= amplitude;
        for (auto &x : sin_ref)
            x /= amplitude_ref;

        size_t min = std::min(sin.size(), sin_ref.size());

        for (size_t n = 0; n < min; ++n)
            power_acc += std::pow(sin_ref[n], 2.0);

        for (size_t n = 0; n < min; ++n)
            errors_int[i] += std::pow((sin_ref[n] - sin[n]), 2.0);

        power_acc /= static_cast<double>(sin.size());
        errors_int[i] /= static_cast<double>(sin.size());

        errors_int[i] = 10 * std::log10(power_acc / errors_int[i]);
    }
    return errors_int;
}