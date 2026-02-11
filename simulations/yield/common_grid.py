from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def make_k_grid_full(n: int, l: float) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    dx = l / n
    freq = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    kx, ky = np.meshgrid(freq, freq, indexing="xy")
    k = np.sqrt(kx**2 + ky**2)
    return kx, ky, k


def make_k_grid_rfft(n: int, l: float) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    dx = l / n
    kx = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    ky = 2.0 * np.pi * np.fft.rfftfreq(n, d=dx)
    kx2d, ky2d = np.meshgrid(kx, ky, indexing="ij")
    k = np.sqrt(kx2d**2 + ky2d**2)
    return kx2d, ky2d, k


def _random_half_spectrum(
    rng: np.random.Generator,
    power: NDArray[np.float64],
) -> NDArray[np.complex128]:
    n, n_half = power.shape
    coeffs = np.zeros((n, n_half), dtype=np.complex128)
    nyquist = n_half - 1

    for i in range(n):
        ni = (-i) % n
        if i > ni:
            continue
        for j in range(n_half):
            variance = power[i, j]
            if variance <= 0.0:
                value = 0.0
            else:
                if i == ni and (j == 0 or j == nyquist):
                    value = rng.normal(0.0, np.sqrt(variance))
                else:
                    std = np.sqrt(variance / 2.0)
                    value = rng.normal(0.0, std) + 1j * rng.normal(0.0, std)
            coeffs[i, j] = value
            if i < ni:
                coeffs[ni, j] = np.conj(value)
    return coeffs


def generate_grf_field(
    rng: np.random.Generator,
    n: int,
    l: float,
    power_spectrum,
) -> NDArray[np.float64]:
    kx, ky, k = make_k_grid_rfft(n, l)
    power = power_spectrum(k)
    power = np.where(np.isfinite(power), power, 0.0)
    power[k == 0.0] = 0.0
    coeffs = _random_half_spectrum(rng, power)
    field = np.fft.irfftn(coeffs, s=(n, n))
    return field.astype(np.float64, copy=False)


def spectral_gradient(
    field: NDArray[np.float64],
    l: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    n = field.shape[0]
    kx, ky, _ = make_k_grid_rfft(n, l)
    field_hat = np.fft.rfftn(field)
    grad_x = np.fft.irfftn(1j * kx * field_hat, s=(n, n))
    grad_y = np.fft.irfftn(1j * ky * field_hat, s=(n, n))
    return grad_x.real.astype(np.float64, copy=False), grad_y.real.astype(np.float64, copy=False)
