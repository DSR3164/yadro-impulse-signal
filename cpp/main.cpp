#include "dsp.hpp"

int main()
{
    std::vector<double> sin;
    std::vector<int16_t> clamped;
    std::vector<int16_t> clamped_upsampled;
    std::vector<double> upsampled;

    generate_sin(sin, 25, 1.6, 100, 10);
    clamp_to(sin, clamped);
    auto temp = upsampling(sin, upsampled);
    upsampling_int(clamped, clamped_upsampled);

    cutoff(sin, 4);
    cutoff(clamped, 4);
    cutoff(upsampled, 4);
    cutoff(clamped_upsampled, 4);

    save_signal("../../files/signal.bin", sin);
    save_signal("../../files/signal_2fc.bin", temp);
    save_signal("../../files/signal_c.bin", clamped);
    save_signal("../../files/signal_u.bin", upsampled);
    save_signal("../../files/signal_cu.bin", clamped_upsampled);

    auto errors_int = perform_int_test();
    auto errors_float = perform_float_test();

    save_signal("../../files/errors_int.bin", errors_int);
    save_signal("../../files/errors_float.bin", errors_float);

    return 0;
}
