#pragma once

#include <cstdint>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <type_traits>

template <typename T> // Что бы не писать функцию под каждый тип
void save_signal(const std::string &filename, const std::vector<T> &signal)
{
    if (signal.empty())
        return;
    std::ofstream ofs(filename, std::ios::out | std::ios::binary); // Записываю в бинарном виде
    if (!ofs.is_open())
        return;
    char type_char = std::is_floating_point_v<T> ? 'f' : (std::is_unsigned_v<T> ? 'u' : 'i');
    uint8_t type_size = sizeof(T);
    ofs.write(&type_char, 1);                                 // Записываю тип данных
    ofs.write(reinterpret_cast<const char *>(&type_size), 1); // Записываю размер типа - далее питон сам поймет в каком виде данные
    ofs.write(reinterpret_cast<const char *>(signal.data()), signal.size() * type_size);
}

template <typename T>
void clamp_to(const std::vector<double> &signal, std::vector<T> &clamped)
{
    if (signal.empty())
        return;

    size_t signal_size = signal.size();
    clamped.resize(signal_size);
    auto min = *std::min_element(signal.begin(), signal.end());

    double max_target = static_cast<double>(std::numeric_limits<T>::max());
    double scale = -(max_target / min);
    for (size_t i = 0; i < signal_size; ++i)
        clamped[i] = static_cast<T>(signal[i] * scale);
}

template <typename T>
void upsample(const std::vector<T> &symbols, std::vector<T> &upsampled, size_t up)
{
    if (upsampled.size() < symbols.size() * up)
    {
        std::cerr << "Wrong upsampled vector size!\n";
        return;
    }
    std::fill(upsampled.begin(), upsampled.end(), 0);

    for (size_t i = 0; i < symbols.size(); ++i)
        upsampled[i * up] = symbols[i];
}

void generate_sin(std::vector<double> &signal, double frequency = 50.0f, double phase = 0.0f, double sample_rate = 100.0f, double seconds = 1.0f);

std::vector<double> upsampling(std::vector<double> &signal, std::vector<double> &upsampled_signal);
std::vector<int16_t> upsampling_int(std::vector<int16_t> &signal, std::vector<int16_t> &upsampled_signal);

std::vector<double> fir_lowpass_sinc_windowed(size_t numTaps, double f_c);
void filter(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &y);
void filter_int(const std::vector<int16_t> &a, const std::vector<double> &b, std::vector<int16_t> &y);