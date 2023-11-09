#include "../kernel.hpp"

string kernel_color_utils() {
	return R(
	) + "#ifdef GRAPHICS" + R(
	) + R(int iron_color(float x) { // coloring scheme (float 0-255 -> int color)
		x = clamp(360.0f - x * 360.0f / 255.0f, 0.0f, 360.0f);
		float r = 255.0f, g = 0.0f, b = 0.0f;
		if (x < 60.0f) { // white - yellow
			g = 255.0f;
			b = 255.0f - 255.0f * x / 60.0f;
		}
		else if (x < 180.0f) { // yellow - red
			g = 255.0f - 255.0f * (x - 60.0f) / 120.0f;
		}
		else if (x < 270.0f) { // red - violet
			r = 255.0f - 255.0f * (x - 180.0f) / 180.0f;
			b = 255.0f * (x - 180.0f) / 90.0f;
		}
		else { // violet - black
			r = 255.0f - 255.0f * (x - 180.0f) / 180.0f;
			b = 255.0f - 255.0f * (x - 270.0f) / 90.0f;
		}
		return (((int)r) << 16) | (((int)g) << 8) | ((int)b);
	}
	) + R(int rainbow_color(float x) { // coloring scheme (float 0-255 -> int color)
		x = clamp(360.0f - x * 360.0f / 255.0f, 0.0f, 360.0f);
		float r = 0.0f, g = 0.0f, b = 0.0f; // black
		if (x < 60.0f) { // red - yellow
			r = 255.0f;
			g = 255.0f * x / 60.0f;
		}
		else if (x >= 60.0f && x < 120.0f) { // yellow - green
			r = 255.0f - 255.0f * (x - 60.0f) / 60.0f;
			g = 255.0f;
		}
		else if (x >= 120.0f && x < 180.0f) { // green - cyan
			g = 255.0f;
			b = 255.0f * (x - 120.0f) / 60.0f;
		}
		else if (x >= 180.0f && x < 240.0f) { // cyan - blue
			g = 255.0f - 255.0f * (x - 180.0f) / 60.0f;
			b = 255.0f;
		}
		else if (x >= 240.0f && x < 300.0f) { // blue - violet
			r = (255.0f * (x - 240.0f) / 60.0f) / 2.0f;
			b = 255.0f;
		}
		else { // violet - black
			r = (255.0f - 255.0f * (x - 300.0f) / 60.0f) / 2.0f;
			b = 255.0f - 255.0f * (x - 300.0f) / 60.0f;
		}
		return (((int)r) << 16) | (((int)g) << 8) | ((int)b);
	}
	) + R(int color_dim(const int c, const float x) {
		const int r = clamp((int)fma((float)((c >> 16) & 255), x, 0.5f), 0, 255);
		const int g = clamp((int)fma((float)((c >> 8) & 255), x, 0.5f), 0, 255);
		const int b = clamp((int)fma((float)(c & 255), x, 0.5f), 0, 255);
		return (r & 255) << 16 | (g & 255) << 8 | (b & 255);
	}
	) + R(int color_mix(const int c1, const int c2, const float w) {
		const uchar4 cc1 = as_uchar4(c1), cc2 = as_uchar4(c2);
		const float3 fc1 = (float3)((float)cc1.x, (float)cc1.y, (float)cc1.z), fc2 = (float3)((float)cc2.x, (float)cc2.y, (float)cc2.z);
		const float3 fcm = fma(w, fc1, fma(1.0f - w, fc2, (float3)(0.5f, 0.5f, 0.5f)));
		return as_int((uchar4)((uchar)fcm.x, (uchar)fcm.y, (uchar)fcm.z, (uchar)0u));
	}
	) + R(int color_mix_3(const int c0, const int c1, const int c2, const float w0, const float w1, const float w2) { // w0+w1+w2 = 1
		const uchar4 cc0 = as_uchar4(c0), cc1 = as_uchar4(c1), cc2 = as_uchar4(c2);
		const float3 fc0 = (float3)((float)cc0.x, (float)cc0.y, (float)cc0.z), fc1 = (float3)((float)cc1.x, (float)cc1.y, (float)cc1.z), fc2 = (float3)((float)cc2.x, (float)cc2.y, (float)cc2.z);
		const float3 fcm = fma(w0, fc0, fma(w1, fc1, fma(w2, fc2, (float3)(0.5f, 0.5f, 0.5f))));
		return as_int((uchar4)((uchar)fcm.x, (uchar)fcm.y, (uchar)fcm.z, (uchar)0u));
	}
	) + R(int hsv_to_rgb(const float h, const float s, const float v) {
		const float c = v * s;
		const float x = c * (1.0f - fabs(fmod(h / 60.0f, 2.0f) - 1.0f));
		const float m = v - c;
		float r = 0.0f, g = 0.0f, b = 0.0f;
		if (0.0f <= h && h < 60.0f) { r = c; g = x; }
		else if (h < 120.0f) { r = x; g = c; }
		else if (h < 180.0f) { g = c; b = x; }
		else if (h < 240.0f) { g = x; b = c; }
		else if (h < 300.0f) { r = x; b = c; }
		else if (h < 360.0f) { r = c; b = x; }
		return (int)((r + m) * 255.0f) << 16 | (int)((g + m) * 255.0f) << 8 | (int)((b + m) * 255.0f);
	}
	) + "#endif" + R( // GRAPHICS
	);
}