from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import numpy as np
from typing import Tuple

def fluctuation_spectrum(x: np.ndarray, y: np.ndarray, z: np.ndarray, Lbox: float) -> Tuple[np.ndarray, np.ndarray]:
    # Интерполяция на регулярную сетку
    grid_size = 64
    xi = np.linspace(0, Lbox, grid_size)
    yi = np.linspace(0, Lbox, grid_size)
    Xi, Yi = np.meshgrid(xi, yi)

    # Используем метод 'linear' для более стабильной интерполяции
    zi = griddata((x, y), z, (Xi, Yi), method='linear', fill_value=0)

    # Функция высоты h(x, y)
    h = zi

    # Сглаживание данных для уменьшения шума
    h_smoothed = gaussian_filter(h, sigma=1)

    # Двумерное Фурье-преобразование
    H = np.fft.fftshift(np.fft.fft2(h_smoothed))

    # Вычисление модуля волнового вектора
    kx = 2 * np.pi * np.fft.fftfreq(grid_size, d=Lbox/grid_size)
    ky = 2 * np.pi * np.fft.fftfreq(grid_size, d=Lbox/grid_size)
    Kx, Ky = np.meshgrid(kx, ky)
    k = np.sqrt(Kx**2 + Ky**2)

    # Бинирование по модулю волнового вектора
    k_bins = np.logspace(np.log10(k[k > 0].min()), np.log10(k.max()), 50)
    k_values = (k_bins[:-1] + k_bins[1:]) / 2  # Центры бинов
    C_k = np.zeros_like(k_values)

    for i in range(len(k_values)):
        mask = (k >= k_bins[i]) & (k < k_bins[i+1])
        if np.any(mask):
            C_k[i] = np.mean(np.abs(H[mask])**2)

    # Удаление NaN и inf
    mask = np.isfinite(C_k) & (k_values > 0)
    k_values_clean = k_values[mask]
    return 1/k_values_clean, C_k[mask]
   