#include "dsp.hpp"

int main()
{
    std::vector<double> sin;
    std::vector<int16_t> clamped;
    std::vector<int16_t> clamped_upsampled;
    std::vector<double> upsampled;

    generate_sin(sin, 25, 0.0, 100, 10);
    clamp_to(sin, clamped);
    auto temp = upsampling(sin, upsampled);
    upsampling_int(clamped, clamped_upsampled);

    save_signal("../../files/signal.bin", sin);
    save_signal("../../files/signal_2fc.bin", temp);
    save_signal("../../files/signal_c.bin", clamped);
    save_signal("../../files/signal_u.bin", upsampled);
    save_signal("../../files/signal_cu.bin", clamped_upsampled);

    return 0;
}
