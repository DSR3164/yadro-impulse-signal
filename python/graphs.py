import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'axes.facecolor': '#1e1e1e',
    'figure.facecolor': '#1e1e1e',
    'text.color': 'white',
    'axes.labelcolor': 'white',
    'xtick.color': 'white',
    'ytick.color': 'white',
    'grid.color': "#747474",
    'font.size': 9,
    'savefig.facecolor': '#1e1e1e',
})

def apply_dark_theme(fig):
    fig.canvas.get_tk_widget().configure(bg='#1e1e1e')
    fig.canvas.manager.window.configure(bg='#1e1e1e')

_original_figure = plt.figure

def figure(*args, **kwargs):
    fig = _original_figure(*args, **kwargs)
    apply_dark_theme(fig)
    return fig

plt.figure = figure

def load_signal(filename):
    with open(filename, 'rb') as f:
        dtype_str = f.read(1).decode() + str(ord(f.read(1)))
        data = np.fromfile(f, dtype=dtype_str)
        print(f"{filename} \t{data.dtype}")
        return data
sin = load_signal("files/signal.bin") # Чистый синус

sin_2fc = load_signal("files/signal_2fc.bin") # Синус с нулями
sin_u = load_signal("files/signal_u.bin") # Синус с повышенной частотой дискретизации

sin_c = load_signal("files/signal_c.bin") # Квантованный синус
sin_cu = load_signal("files/signal_cu.bin") # Квантованный синус с повышенной частотой дискретизации

h = load_signal("files/h.bin") # Импульсная характеристика фильтра

errors_float = load_signal("files/errors_float.bin")
errors_int = load_signal("files/errors_int.bin")

fc1 = 1 / 100
fc2 = fc1 / 2

sin_normalized_to_int16 = sin / np.max(np.abs(sin)) * np.max(np.abs(sin_c))
sin_с_normalized_to_int16 = sin_c

len_min = min(len(sin_с_normalized_to_int16), len(sin_normalized_to_int16))
sin_cu_t = sin_с_normalized_to_int16[:len_min]
sin_u_t = sin_normalized_to_int16[:len_min]

error = sin_u_t - sin_cu_t
signal_power = np.mean(np.float64(sin_cu_t) ** 2)
noise_power = np.mean(error ** 2)

snr = 10 * np.log10(signal_power / noise_power)
print(f"SQNR (мощность сигнала к ошибке квантования): {snr:.2f}dB ")

Fsinc = np.fft.fftshift(np.fft.fft(h))
Xsinc = np.fft.fftshift(np.fft.fftfreq(len(h), fc1))

F = np.fft.fftshift(np.fft.fft(sin))
X = np.fft.fftshift(np.fft.fftfreq(len(sin), fc1))

Ff = np.fft.fftshift(np.fft.fft(sin_u))
Xf = np.fft.fftshift(np.fft.fftfreq(len(sin_u), fc2))

Fcf = np.fft.fftshift(np.fft.fft(sin_cu))
Xcf = np.fft.fftshift(np.fft.fftfreq(len(sin_cu), fc2))

Ffc2 = np.fft.fftshift(np.fft.fft(sin_2fc))
Xfc2 = np.fft.fftshift(np.fft.fftfreq(len(sin_2fc), fc2))

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(errors_float, label="Errors")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude dB")
plt.tight_layout()
plt.savefig('docs/Errors float')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(errors_int, label="Errors")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude dB")
plt.tight_layout()
plt.savefig('docs/Errors int')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(errors_float, label="Errors float")
plt.plot(errors_int, label="Errors int")
plt.grid()
plt.legend()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude dB")
plt.tight_layout()
plt.savefig('docs/Errors')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(Xsinc, 20 * np.log10(np.abs(Fsinc) + 1e-12), label="Filter spectrum")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.tight_layout()
plt.savefig('docs/Filter spectrum')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(X, 20 * np.log10(np.abs(F) + 1e-12), label="Original sin spectrum")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.tight_layout()
plt.savefig('docs/Original sin spectrum')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(Xf, 20 * np.log10(np.abs(Ff) + 1e-12), label="Upsampled sin spectrum")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.tight_layout()
plt.savefig('docs/Upsampled sin spectrum')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(h, label="Filter")
plt.grid()
plt.tight_layout()
plt.savefig('docs/Filter')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(sin, "o-", label="Original sin")
plt.grid()
plt.tight_layout()
plt.xlim([0, 1 / fc1])
plt.savefig('docs/Original sin')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(Xfc2, 20 * np.log10(np.abs(Ffc2) + 1e-12), label="Original sin with 2fc")
plt.grid()
plt.tight_layout()
plt.savefig('docs/Original sin with 2fc spectrum')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(sin_c, label="Clamped sin")
plt.grid()
plt.tight_layout()
plt.xlim([0, 1 / fc1])
plt.savefig('docs/Clamped sin')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(sin_u, "o-", label="Upsampled sin")
plt.grid()
plt.tight_layout()
plt.xlim([0, 1 / fc2])
plt.savefig('docs/Upsampled sin')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(sin_cu, "o-", label="Clamped Upsampled sin")
plt.grid()
plt.tight_layout()
plt.xlim([0, 1 / fc2])
plt.savefig('docs/Clamped Upsampled sin')

plt.figure(figsize = [12, 4], dpi=100)
plt.plot(Xcf, 20 * np.log10(np.abs(Fcf) + 1e-12), label="Clamped Upsampled sin spectrum")
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.tight_layout()
plt.savefig('docs/Clamped Upsampled sin spectrum')
