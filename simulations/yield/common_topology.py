from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def bilinear_interpolate(
    field: NDArray[np.float64],
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    x_min: float,
    y_min: float,
    dx: float,
) -> NDArray[np.float64]:
    n = field.shape[0]
    j = (x - x_min) / dx
    i = (y - y_min) / dx

    i0 = np.floor(i).astype(int)
    j0 = np.floor(j).astype(int)
    i0 = np.clip(i0, 0, n - 2)
    j0 = np.clip(j0, 0, n - 2)
    i1 = i0 + 1
    j1 = j0 + 1

    di = i - i0
    dj = j - j0

    f00 = field[i0, j0]
    f10 = field[i1, j0]
    f01 = field[i0, j1]
    f11 = field[i1, j1]

    f0 = f00 * (1.0 - di) + f10 * di
    f1 = f01 * (1.0 - di) + f11 * di
    return f0 * (1.0 - dj) + f1 * dj


def sample_complex_on_circle(
    real_field: NDArray[np.float64],
    imag_field: NDArray[np.float64],
    x0: float,
    y0: float,
    radius: float,
    num_points: int,
    x_min: float,
    y_min: float,
    dx: float,
) -> NDArray[np.complex128]:
    angles = np.linspace(0.0, 2.0 * np.pi, num_points, endpoint=False)
    xs = x0 + radius * np.cos(angles)
    ys = y0 + radius * np.sin(angles)
    real_samples = bilinear_interpolate(real_field, xs, ys, x_min, y_min, dx)
    imag_samples = bilinear_interpolate(imag_field, xs, ys, x_min, y_min, dx)
    return real_samples + 1j * imag_samples


def sample_complex_from_polar_on_circle(
    amplitude: NDArray[np.float64],
    phase: NDArray[np.float64],
    x0: float,
    y0: float,
    radius: float,
    num_points: int,
    x_min: float,
    y_min: float,
    dx: float,
) -> NDArray[np.complex128]:
    angles = np.linspace(0.0, 2.0 * np.pi, num_points, endpoint=False)
    xs = x0 + radius * np.cos(angles)
    ys = y0 + radius * np.sin(angles)
    w_samples = bilinear_interpolate(amplitude, xs, ys, x_min, y_min, dx)
    phi_samples = bilinear_interpolate(phase, xs, ys, x_min, y_min, dx)
    return w_samples * np.exp(1j * phi_samples)


def winding_number_from_samples(samples: NDArray[np.complex128]) -> int:
    phases = np.angle(samples)
    deltas = np.diff(phases, append=phases[0])
    deltas = (deltas + np.pi) % (2.0 * np.pi) - np.pi
    total = float(np.sum(deltas))
    return int(np.rint(total / (2.0 * np.pi)))


def min_abs_in_disk(
    field: NDArray[np.float64],
    r_grid: NDArray[np.float64],
    radius: float,
) -> float:
    mask = r_grid <= radius
    if not np.any(mask):
        return float("inf")
    return float(np.min(np.abs(field[mask])))
