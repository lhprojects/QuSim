
#pragma once

#include <Windows.h>
#include <eigen/Eigen/Dense>
#include <QuSim.h>

inline void DrawPsi2D(Gdiplus::Graphics& graphics,
	Eigen::MatrixXcd const& psi,
	long left, long top, long w, long h)
{
	auto it = std::max_element(psi.data(), psi.data() + psi.size(), [](Complex l, Complex r) {
		return abs2(l) < abs2(r);
		});
	Real vmax = 1.1 * abs(*it);

	Bitmap bitmap(w, h);
	Rect lockrect(0, 0, w, h);
	BitmapData bitData;
	bitmap.LockBits(&lockrect, ImageLockModeWrite, PixelFormat32bppARGB, &bitData);

	uint8_t* da = (uint8_t*)bitData.Scan0;
	int nx = psi.rows();
	int ny = psi.cols();

	for (int i = 0; i < (int)bitData.Width; ++i) {
		for (int j = 0; j < (int)bitData.Height; ++j) {

            int ix = (int)(1.0 * i / bitData.Width * nx);
            int iy = (int)(1.0 * j / bitData.Height * ny);

			Complex psi_ = psi(ix, iy);
			int b = (int)(abs(psi_) / vmax * 255);
			b = (b >= 255) ? 255 : b;
			double a = 0.5 * real(psi_) / abs(psi_) + 0.5;

			da[j * bitData.Stride + 4 * i] = a * 255; // blue
			da[j * bitData.Stride + 4 * i + 1] = (1 - a) * 255; // green
			da[j * bitData.Stride + 4 * i + 2] = 0; // red
			da[j * bitData.Stride + 4 * i + 3] = b;
		}
	}

	bitmap.UnlockBits(&bitData);
	CachedBitmap cb(&bitmap, &graphics);
	graphics.DrawCachedBitmap(&cb, 0, 0);



}

inline void DrawPotential2D(Gdiplus::Graphics& graphics,
	Eigen::MatrixXcd const& V,
	long left, long top, long w, long h)
{

	auto it = std::max_element(V.data(), V.data() + V.size(), [](Complex a, Complex b) {
		return abs(a) < abs(b);
		});

	Real vmax = 1.1*abs(*it);

	Bitmap bitmap(w, h);
	Rect lockrect(0, 0, w, h);
	BitmapData bitData;
	bitmap.LockBits(&lockrect, ImageLockModeWrite, PixelFormat32bppARGB, &bitData);

	uint8_t* da = (uint8_t*)bitData.Scan0;

	for (int i = 0; i < (int)bitData.Width; ++i) {
		for (int j = 0; j < (int)bitData.Height; ++j) {

			int nx = V.rows();
			int ny = V.cols();

			int ix = (int)(1.0 * i / bitData.Width * nx);
			int iy = (int)(1.0 * j / bitData.Height * ny);

			int b = (int)(abs(V(ix, iy).real()) / vmax * 255);
			int r = (int)(abs(V(ix, iy).imag()) / vmax * 255);

            da[j * bitData.Stride + 4 * i] = 0;
            da[j * bitData.Stride + 4 * i + 1] = r;
            da[j * bitData.Stride + 4 * i + 2] = b;
            da[j * bitData.Stride + 4 * i + 3] = 255;
		}
	}

	bitmap.UnlockBits(&bitData);
	CachedBitmap cb(&bitmap, &graphics);
	graphics.DrawCachedBitmap(&cb, 0, 0);

}
