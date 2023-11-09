#include "../kernel.hpp"

string kernel_graphics_utils() {
	return R( 
	) + "#ifdef GRAPHICS" + R(
	) + R(int lighting(const int c, const float3 p, const float3 normal, const float* camera_cache) { // calculate lighting of triangle
		const float dis = camera_cache[1]; // fetch camera parameters (rotation matrix, camera position, etc.)
		const float posx = camera_cache[2] - def_domain_offset_x;
		const float posy = camera_cache[3] - def_domain_offset_y;
		const float posz = camera_cache[4] - def_domain_offset_z;
		const float Rzx = camera_cache[11];
		const float Rzy = camera_cache[12];
		const float Rzz = camera_cache[13];
		const uchar cr = c >> 16 & 255, cg = c >> 8 & 255, cb = c & 255;
		const float nl2 = sq(normal.x) + sq(normal.y) + sq(normal.z); // only one native_sqrt instead of two
		const float dx = p.x - fma(Rzx, dis, posx); // direction of light source is viewing direction
		const float dy = p.y - fma(Rzy, dis, posy);
		const float dz = p.z - fma(Rzz, dis, posz);
		const float dl2 = sq(dx) + sq(dy) + sq(dz);
		const float br = max(1.5f * fabs(normal.x * dx + normal.y * dy + normal.z * dz) * rsqrt(nl2 * dl2), 0.3f);
		return min((int)(br * cr), 255) << 16 | min((int)(br * cg), 255) << 8 | min((int)(br * cb), 255);
	}
	) + R(bool is_off_screen(const int x, const int y, const int stereo) {
		switch (stereo) {
		default: return x < 0 || x >= def_screen_width || y < 0 || y >= def_screen_height; // entire screen
		case -1: return x < 0 || x >= def_screen_width / 2 || y < 0 || y >= def_screen_height; // left half
		case +1: return x < def_screen_width / 2 || x >= def_screen_width || y < 0 || y >= def_screen_height; // right half
		}
	}
	) + R(void draw(const int x, const int y, const float z, const int color, global int* bitmap, volatile global int* zbuffer, const int stereo) {
		const int index = x + y * def_screen_width, iz = (int)(z * (2147483647.0f / 10000.0f)); // use int z-buffer and atomic_max to minimize noise in image
		if (!is_off_screen(x, y, stereo) && iz > atomic_max(&zbuffer[index], iz)) bitmap[index] = color; // only draw if point is on screen and first in zbuffer
	}
	) + R(bool convert(int* rx, int* ry, float* rz, const float3 p, const float* camera_cache, const int stereo) { // 3D -> 2D
		const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
		const float dis = camera_cache[1];
		const float posx = camera_cache[2];
		const float posy = camera_cache[3];
		const float posz = camera_cache[4];
		const float Rxx = camera_cache[5];
		const float Rxy = camera_cache[6];
		const float Rxz = camera_cache[7];
		const float Ryx = camera_cache[8];
		const float Ryy = camera_cache[9];
		const float Ryz = camera_cache[10];
		const float Rzx = camera_cache[11];
		const float Rzy = camera_cache[12];
		const float Rzz = camera_cache[13];
		const float eye_distance = vload_half(28, (half*)camera_cache);
		float3 t, r;
		t.x = p.x + def_domain_offset_x - posx - (float)stereo * eye_distance / zoom * Rxx; // transformation
		t.y = p.y + def_domain_offset_y - posy - (float)stereo * eye_distance / zoom * Rxy;
		t.z = p.z + def_domain_offset_z - posz;
		r.z = Rzx * t.x + Rzy * t.y + Rzz * t.z; // z-position for z-buffer
		const float rs = zoom * dis / (dis - r.z * zoom); // perspective (reciprocal is more efficient)
		if (rs <= 0.0f) return false; // point is behins camera
		const float tv = ((as_int(camera_cache[14]) >> 30) & 0x1) && stereo != 0 ? 0.5f : 1.0f;
		r.x = ((Rxx * t.x + Rxy * t.y + Rxz * t.z) * rs + (float)stereo * eye_distance) * tv + (0.5f + (float)stereo * 0.25f) * (float)def_screen_width; // x position on screen
		r.y = (Ryx * t.x + Ryy * t.y + Ryz * t.z) * rs + 0.5f * (float)def_screen_height; // y position on screen
		*rx = (int)(r.x + 0.5f);
		*ry = (int)(r.y + 0.5f);
		*rz = r.z;
		return true;
	}
	) + R(void convert_circle(float3 p, const float r, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
		int rx, ry; float rz;
		if (convert(&rx, &ry, &rz, p, camera_cache, stereo)) {
			const float zoom = camera_cache[0];
			const float dis = camera_cache[1];
			const float rs = zoom * dis / (dis - rz * zoom);
			const int radius = (int)(rs * r + 0.5f);
			switch (stereo) {
			default: if (rx < -radius || rx >= (int)def_screen_width + radius || ry < -radius || ry >= (int)def_screen_height + radius) return; break; // cancel drawing if circle is off screen
			case -1: if (rx < -radius || rx >= (int)def_screen_width / 2 + radius || ry < -radius || ry >= (int)def_screen_height + radius) return; break;
			case +1: if (rx < (int)def_screen_width / 2 - radius || rx >= (int)def_screen_width + radius || ry < -radius || ry >= (int)def_screen_height + radius) return; break;
			}
			int d = -radius, x = radius, y = 0; // Bresenham algorithm for circle
			while (x >= y) {
				draw(rx + x, ry + y, rz, color, bitmap, zbuffer, stereo);
				draw(rx - x, ry + y, rz, color, bitmap, zbuffer, stereo);
				draw(rx + x, ry - y, rz, color, bitmap, zbuffer, stereo);
				draw(rx - x, ry - y, rz, color, bitmap, zbuffer, stereo);
				draw(rx + y, ry + x, rz, color, bitmap, zbuffer, stereo);
				draw(rx - y, ry + x, rz, color, bitmap, zbuffer, stereo);
				draw(rx + y, ry - x, rz, color, bitmap, zbuffer, stereo);
				draw(rx - y, ry - x, rz, color, bitmap, zbuffer, stereo);
				d += 2 * y + 1;
				y++;
				if (d > 0) d -= 2 * (--x);
			}
		}
	}
	) + R(void convert_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
		int r0x, r0y, r1x, r1y; float r0z, r1z;
		if (convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo)
			&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo))) { // cancel drawing if both points are off screen
			int x = r0x, y = r0y; // Bresenham algorithm
			const float z = 0.5f * (r0z + r1z); // approximate line z position for each pixel to be equal
			const int dx = abs(r1x - r0x), sx = 2 * (r0x < r1x) - 1;
			const int dy = -abs(r1y - r0y), sy = 2 * (r0y < r1y) - 1;
			int err = dx + dy;
			while (x != r1x || y != r1y) {
				draw(x, y, z, color, bitmap, zbuffer, stereo);
				const int e2 = 2 * err;
				if (e2 > dy) { err += dy; x += sx; }
				if (e2 < dx) { err += dx; y += sy; }
			}
		}
	}
	) + R(void convert_triangle(float3 p0, float3 p1, float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
		int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
		if (convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
			&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
			if (r0x * (r1y - r2y) + r1x * (r2y - r0y) + r2x * (r0y - r1y) > 100000 || (r0y == r1y && r0y == r2y)) return; // return for large triangle area or degenerate triangles
			if (r0y > r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; } // sort vertices ascending by y
			if (r0y > r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; }
			if (r1y > r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; }
			const float z = (r0z + r1z + r2z) / 3.0f; // approximate triangle z position for each pixel to be equal
			for (int y = r0y; y < r1y; y++) { // Bresenham algorithm (lower triangle half)
				const int xA = r0x + (r2x - r0x) * (y - r0y) / (r2y - r0y);
				const int xB = r0x + (r1x - r0x) * (y - r0y) / (r1y - r0y);
				for (int x = min(xA, xB); x < max(xA, xB); x++) {
					draw(x, y, z, color, bitmap, zbuffer, stereo);
				}
			}
			for (int y = r1y; y < r2y; y++) { // Bresenham algorithm (upper triangle half)
				const int xA = r0x + (r2x - r0x) * (y - r0y) / (r2y - r0y);
				const int xB = r1x + (r2x - r1x) * (y - r1y) / (r2y - r1y);
				for (int x = min(xA, xB); x < max(xA, xB); x++) {
					draw(x, y, z, color, bitmap, zbuffer, stereo);
				}
			}
		}
	}
	) + R(void convert_triangle_interpolated(float3 p0, float3 p1, float3 p2, int c0, int c1, int c2, const float* camera_cache, global int* bitmap, global int* zbuffer, const int stereo) { // 3D -> 2D
		int r0x, r0y, r1x, r1y, r2x, r2y; float r0z, r1z, r2z;
		if (convert(&r0x, &r0y, &r0z, p0, camera_cache, stereo) && convert(&r1x, &r1y, &r1z, p1, camera_cache, stereo) && convert(&r2x, &r2y, &r2z, p2, camera_cache, stereo)
			&& !(is_off_screen(r0x, r0y, stereo) && is_off_screen(r1x, r1y, stereo) && is_off_screen(r2x, r2y, stereo))) { // cancel drawing if all points are off screen
			if (r0x * (r1y - r2y) + r1x * (r2y - r0y) + r2x * (r0y - r1y) > 100000 || (r0y == r1y && r0y == r2y)) return; // return for large triangle area or degenerate triangles
			if (r0y > r1y) { const int xt = r0x; const int yt = r0y; r0x = r1x; r0y = r1y; r1x = xt; r1y = yt; const int ct = c0; c0 = c1; c1 = ct; } // sort vertices ascending by y
			if (r0y > r2y) { const int xt = r0x; const int yt = r0y; r0x = r2x; r0y = r2y; r2x = xt; r2y = yt; const int ct = c0; c0 = c2; c2 = ct; }
			if (r1y > r2y) { const int xt = r1x; const int yt = r1y; r1x = r2x; r1y = r2y; r2x = xt; r2y = yt; const int ct = c1; c1 = c2; c2 = ct; }
			const float z = (r0z + r1z + r2z) / 3.0f; // approximate triangle z position for each pixel to be equal
			const float d = (float)((r1y - r2y) * (r0x - r2x) + (r2x - r1x) * (r0y - r2y));
			for (int y = r0y; y < r1y; y++) { // Bresenham algorithm (lower triangle half)
				const int xA = r0x + (r2x - r0x) * (y - r0y) / (r2y - r0y);
				const int xB = r0x + (r1x - r0x) * (y - r0y) / (r1y - r0y);
				for (int x = min(xA, xB); x < max(xA, xB); x++) {
					const float w0 = (float)((r1y - r2y) * (x - r2x) + (r2x - r1x) * (y - r2y)) / d; // barycentric coordinates
					const float w1 = (float)((r2y - r0y) * (x - r2x) + (r0x - r2x) * (y - r2y)) / d;
					const float w2 = 1.0f - w0 - w1;
					const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
					draw(x, y, z, color, bitmap, zbuffer, stereo);
				}
			}
			for (int y = r1y; y < r2y; y++) { // Bresenham algorithm (upper triangle half)
				const int xA = r0x + (r2x - r0x) * (y - r0y) / (r2y - r0y);
				const int xB = r1x + (r2x - r1x) * (y - r1y) / (r2y - r1y);
				for (int x = min(xA, xB); x < max(xA, xB); x++) {
					const float w0 = (float)((r1y - r2y) * (x - r2x) + (r2x - r1x) * (y - r2y)) / d; // barycentric coordinates
					const float w1 = (float)((r2y - r0y) * (x - r2x) + (r0x - r2x) * (y - r2y)) / d;
					const float w2 = 1.0f - w0 - w1;
					const int color = color_mix_3(c0, c1, c2, w0, w1, w2); // interpolate color
					draw(x, y, z, color, bitmap, zbuffer, stereo);
				}
			}
		}
	}
	) + R(void draw_point(const float3 p, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
		const bool vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		int rx, ry; float rz;
		if (!vr) {
			if (convert(&rx, &ry, &rz, p, camera_cache, 0)) draw(rx, ry, rz, color, bitmap, zbuffer, 0);
		}
		else {
			if (convert(&rx, &ry, &rz, p, camera_cache, -1)) draw(rx, ry, rz, color, bitmap, zbuffer, -1); // left eye
			if (convert(&rx, &ry, &rz, p, camera_cache, +1)) draw(rx, ry, rz, color, bitmap, zbuffer, +1); // right eye
		}
	}
	) + R(void draw_circle(const float3 p, const float r, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
		const bool vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		if (!vr) {
			convert_circle(p, r, color, camera_cache, bitmap, zbuffer, 0);
		}
		else {
			convert_circle(p, r, color, camera_cache, bitmap, zbuffer, -1); // left eye
			convert_circle(p, r, color, camera_cache, bitmap, zbuffer, +1); // right eye
		}
	}
	) + R(void draw_line(const float3 p0, const float3 p1, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
		const bool vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		if (!vr) {
			convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, 0);
		}
		else {
			convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, -1); // left eye
			convert_line(p0, p1, color, camera_cache, bitmap, zbuffer, +1); // right eye
		}
	}
	) + R(void draw_triangle(const float3 p0, const float3 p1, const float3 p2, const int color, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
		const bool vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		if (!vr) {
			convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, 0);
		}
		else {
			convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, -1); // left eye
			convert_triangle(p0, p1, p2, color, camera_cache, bitmap, zbuffer, +1); // right eye
		}
	}
	) + R(void draw_triangle_interpolated(const float3 p0, const float3 p1, const float3 p2, const int c0, const int c1, const int c2, const float* camera_cache, global int* bitmap, global int* zbuffer) { // 3D -> 2D
		const bool vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		if (!vr) {
			convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, 0);
		}
		else {
			convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, -1); // left eye
			convert_triangle_interpolated(p0, p1, p2, c0, c1, c2, camera_cache, bitmap, zbuffer, +1); // right eye
		}
	}
	) + R(kernel void graphics_clear(global int* bitmap, global int* zbuffer) {
		const uint n = get_global_id(0);
		bitmap[n] = def_background_color; // black background = 0x000000, use 0xFFFFFF for white background
		zbuffer[n] = -2147483648;
	}
	) + R(constant ushort edge_table_data[128] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, mirror symmetry applied, makes marching-cubes 31% faster
	0x000, 0x109, 0x203, 0x30A, 0x406, 0x50F, 0x605, 0x70C, 0x80C, 0x905, 0xA0F, 0xB06, 0xC0A, 0xD03, 0xE09, 0xF00,
	0x190, 0x099, 0x393, 0x29A, 0x596, 0x49F, 0x795, 0x69C, 0x99C, 0x895, 0xB9F, 0xA96, 0xD9A, 0xC93, 0xF99, 0xE90,
	0x230, 0x339, 0x033, 0x13A, 0x636, 0x73F, 0x435, 0x53C, 0xA3C, 0xB35, 0x83F, 0x936, 0xE3A, 0xF33, 0xC39, 0xD30,
	0x3A0, 0x2A9, 0x1A3, 0x0AA, 0x7A6, 0x6AF, 0x5A5, 0x4AC, 0xBAC, 0xAA5, 0x9AF, 0x8A6, 0xFAA, 0xEA3, 0xDA9, 0xCA0,
	0x460, 0x569, 0x663, 0x76A, 0x066, 0x16F, 0x265, 0x36C, 0xC6C, 0xD65, 0xE6F, 0xF66, 0x86A, 0x963, 0xA69, 0xB60,
	0x5F0, 0x4F9, 0x7F3, 0x6FA, 0x1F6, 0x0FF, 0x3F5, 0x2FC, 0xDFC, 0xCF5, 0xFFF, 0xEF6, 0x9FA, 0x8F3, 0xBF9, 0xAF0,
	0x650, 0x759, 0x453, 0x55A, 0x256, 0x35F, 0x055, 0x15C, 0xE5C, 0xF55, 0xC5F, 0xD56, 0xA5A, 0xB53, 0x859, 0x950,
	0x7C0, 0x6C9, 0x5C3, 0x4CA, 0x3C6, 0x2CF, 0x1C5, 0x0CC, 0xFCC, 0xEC5, 0xDCF, 0xCC6, 0xBCA, 0xAC3, 0x9C9, 0x8C0
		};
	) + R(constant uchar triangle_table_data[1920] = { // source: Paul Bourke, http://paulbourke.net/geometry/polygonise/, termination value 15, bit packed
		255,255,255,255,255,255,255, 15, 56,255,255,255,255,255,255, 16,249,255,255,255,255,255, 31, 56,137,241,255,255,255,255, 33,250,255,255,255,255,255, 15, 56, 33,250,255,255,255,255, 41, 10,146,
		255,255,255,255, 47, 56,162,168,137,255,255,255,179,242,255,255,255,255,255, 15, 43,184,240,255,255,255,255,145, 32,179,255,255,255,255, 31, 43,145,155,184,255,255,255,163,177, 58,255,255,255,
		255, 15, 26,128,138,171,255,255,255,147, 48,155,171,249,255,255,159,168,138,251,255,255,255,255,116,248,255,255,255,255,255, 79,  3, 55,244,255,255,255,255, 16,137,116,255,255,255,255, 79,145,
		116,113, 19,255,255,255, 33,138,116,255,255,255,255, 63,116,  3, 20,162,255,255,255, 41,154, 32, 72,247,255,255, 47,154,146, 39, 55,151,244,255, 72, 55, 43,255,255,255,255,191,116, 43, 36, 64,
		255,255,255,  9,129,116, 50,251,255,255, 79,183, 73,155, 43, 41,241,255,163, 49,171,135,244,255,255, 31,171, 65, 27, 64,183,244,255,116,152,176,185,186, 48,255, 79,183,180,153,171,255,255,255,
		 89,244,255,255,255,255,255,159, 69,128,243,255,255,255,255, 80, 20,  5,255,255,255,255,143, 69, 56, 53, 81,255,255,255, 33,154, 69,255,255,255,255, 63,128, 33, 74, 89,255,255,255, 37, 90, 36,
		  4,242,255,255, 47, 90, 35, 53, 69, 67,248,255, 89, 36,179,255,255,255,255, 15, 43,128, 75, 89,255,255,255, 80,  4, 81, 50,251,255,255, 47, 81, 82, 40,184,132,245,255, 58,171, 49, 89,244,255,
		255, 79, 89,128,129, 26,184,250,255, 69, 80,176,181,186, 48,255, 95,132,133,170,184,255,255,255,121, 88,151,255,255,255,255,159,  3, 89, 83, 55,255,255,255,112,  8,113, 81,247,255,255, 31, 53,
		 83,247,255,255,255,255,121,152,117, 26,242,255,255,175, 33, 89, 80,  3,117,243,255,  8,130, 82, 88,167, 37,255, 47, 90, 82, 51,117,255,255,255,151,117,152,179,242,255,255,159,117,121,146,  2,
		114,251,255, 50, 11,129,113, 24,117,255,191, 18, 27,119, 81,255,255,255, 89,136,117, 26,163,179,255, 95,  7,  5,121, 11,  1,186, 10,171,176, 48, 90,128,112,117,176, 90,183,245,255,255,255,255,
		106,245,255,255,255,255,255, 15, 56,165,246,255,255,255,255,  9, 81,106,255,255,255,255, 31, 56,145, 88,106,255,255,255, 97, 37, 22,255,255,255,255, 31, 86, 33, 54,128,255,255,255,105,149, 96,
		 32,246,255,255, 95,137,133, 82, 98, 35,248,255, 50,171, 86,255,255,255,255,191,128, 43,160, 86,255,255,255, 16, 41,179,165,246,255,255, 95,106,145,146, 43,137,251,255, 54,107, 53, 21,243,255,
		255, 15,184,176,  5, 21,181,246,255,179,  6, 99, 96,  5,149,255,111,149,150,187,137,255,255,255,165, 70,135,255,255,255,255, 79,  3,116, 99,165,255,255,255,145, 80,106, 72,247,255,255,175, 86,
		145, 23, 55,151,244,255, 22, 98, 21,116,248,255,255, 31, 82, 37, 54, 64, 67,247,255, 72,151, 80, 96,  5, 98,255,127,147,151, 52,146,149, 38,150,179,114, 72,106,245,255,255, 95,106,116, 66,  2,
		114,251,255, 16, 73,135, 50, 91,106,255,159, 18,185,146,180,183, 84,106, 72, 55, 91, 83, 81,107,255, 95,177,181, 22,176,183,  4,180, 80,  9, 86, 48,182, 54, 72,103,149,150, 75,151,183,249,255,
		 74,105,164,255,255,255,255, 79,106,148, 10, 56,255,255,255, 10,161,  6, 70,240,255,255,143, 19, 24,134, 70, 22,250,255, 65, 25, 66, 98,244,255,255, 63,128, 33, 41,148, 98,244,255, 32, 68, 98,
		255,255,255,255,143, 35, 40, 68, 98,255,255,255, 74,169, 70, 43,243,255,255, 15, 40,130, 75,169,164,246,255,179,  2, 97, 96,100,161,255,111, 20, 22, 74, 24, 18,139, 27,105,148, 99, 25,179, 54,
		255,143, 27, 24,176, 22, 25,100, 20,179, 54,  6, 96,244,255,255,111,132,107,248,255,255,255,255,167,118,168,152,250,255,255, 15, 55,160,  7,169,118,250,255,106, 23,122,113, 24,  8,255,175,118,
		122, 17, 55,255,255,255, 33, 22,134,129,137,118,255, 47,150,146, 97,151,144,115,147,135,112, 96,  6,242,255,255,127, 35,118,242,255,255,255,255, 50,171,134,138,137,118,255, 47,112,114, 11,121,
		118,154,122,129, 16,135,161,103,167, 50,187, 18, 27,167, 22,118,241,255,152,134,118, 25,182, 54, 49,  6, 25,107,247,255,255,255,255,135,112, 96,179,176,  6,255,127,107,255,255,255,255,255,255,
		103,251,255,255,255,255,255, 63,128,123,246,255,255,255,255, 16,185,103,255,255,255,255,143,145, 56,177,103,255,255,255, 26, 98,123,255,255,255,255, 31,162,  3,104,123,255,255,255,146, 32,154,
		182,247,255,255,111,123,162,163, 56,154,248,255, 39, 99,114,255,255,255,255,127,128,103, 96,  2,255,255,255,114, 38,115, 16,249,255,255, 31, 38,129, 22,137,120,246,255,122,166,113, 49,247,255,
		255,175,103,113, 26,120,  1,248,255, 48,  7,167,160,105,122,255,127,166,167,136,154,255,255,255,134,180,104,255,255,255,255, 63,182,  3,  6,100,255,255,255,104,139,100,  9,241,255,255,159,100,
		105,147, 19, 59,246,255,134,100,139,162,241,255,255, 31,162,  3, 11,182, 64,246,255,180, 72,182, 32, 41,154,255,175, 57, 58,146, 52, 59, 70, 54, 40,131, 36,100,242,255,255, 15, 36,100,242,255,
		255,255,255,145, 32, 67, 66, 70,131,255, 31, 73, 65, 34,100,255,255,255, 24,131, 22, 72,102, 26,255,175,  1, 10,102, 64,255,255,255,100, 67,131,166,  3,147,154,163, 73,166,244,255,255,255,255,
		148,117,182,255,255,255,255, 15, 56,148,181,103,255,255,255,  5, 81,  4,103,251,255,255,191,103, 56, 52, 69, 19,245,255, 89,164, 33,103,251,255,255,111,123, 33, 10, 56,148,245,255,103, 91,164,
		 36, 74, 32,255, 63,132, 83, 52, 82, 90,178,103, 39,115, 38, 69,249,255,255,159, 69,128,  6, 38,134,247,255, 99, 50,103, 81, 80,  4,255,111,130,134, 39,129,132, 21,133, 89,164, 97,113, 22,115,
		255, 31,166,113, 22,112,120,144, 69,  4, 74, 90, 48,106,122,115,122,166,167, 88,164,132,250,255,150,101,155,139,249,255,255, 63,182, 96,  3,101,144,245,255,176,  8,181, 16, 85,182,255,111, 59,
		 54, 85, 19,255,255,255, 33,154,181,185,184,101,255, 15, 59, 96, 11,105,101, 25,162,139,181,101,  8,165, 37, 32,101, 59, 54, 37, 58, 90,243,255,133, 89,130,101, 50, 40,255,159,101,105,  0, 38,
		255,255,255, 81, 24,  8,101, 56, 40, 38, 24,101, 18,246,255,255,255,255, 49, 22,166,131, 86,150,152,166,  1, 10,150,  5,101,240,255, 48, 88,166,255,255,255,255,175,101,255,255,255,255,255,255,
		 91,122,181,255,255,255,255,191,165,123,133,  3,255,255,255,181, 87,186,145,240,255,255,175, 87,186,151, 24, 56,241,255, 27,178, 23, 87,241,255,255, 15, 56, 33, 23, 87, 39,251,255,121,149,114,
		  9, 34,123,255,127, 37, 39, 91, 41, 35,152, 40, 82, 42, 83,115,245,255,255,143,  2, 88,130, 87, 42,245,255,  9, 81, 58, 53, 55, 42,255,159, 40, 41,129, 39, 42,117, 37, 49, 53, 87,255,255,255,
		255, 15,120,112, 17, 87,255,255,255,  9,147, 83, 53,247,255,255,159,120,149,247,255,255,255,255,133, 84,138,186,248,255,255, 95, 64,181, 80,186, 59,240,255, 16,137,164,168,171, 84,255,175, 75,
		 74,181, 67, 73, 49, 65, 82, 33, 88,178, 72,133,255, 15,180,176, 67,181,178, 81,177, 32,  5,149,178, 69,133,139,149, 84,178,243,255,255,255,255, 82, 58, 37, 67, 53, 72,255, 95, 42, 37, 68,  2,
		255,255,255,163, 50,165,131, 69,133, 16, 89, 42, 37, 20, 41, 73,242,255, 72,133, 53, 83,241,255,255, 15, 84,  1,245,255,255,255,255, 72,133, 53,  9,  5, 83,255,159, 84,255,255,255,255,255,255,
		180, 71,185,169,251,255,255, 15, 56,148,151,123,169,251,255,161, 27, 75, 65,112,180,255, 63, 65, 67, 24, 74, 71,171, 75,180,151, 75, 41,155, 33,255,159, 71,185,151,177,178,  1, 56,123,180, 36,
		 66,240,255,255,191, 71, 75,130, 67, 35,244,255,146, 42,151, 50,119,148,255,159,122,121,164,114,120, 32,112,115, 58, 42, 71, 26, 10,  4, 26, 42,120,244,255,255,255,255,148, 65,113, 23,243,255,
		255, 79, 25, 20,  7, 24,120,241,255,  4,115, 52,255,255,255,255, 79,120,255,255,255,255,255,255,169,168,139,255,255,255,255, 63,144,147,187,169,255,255,255, 16, 10,138,168,251,255,255, 63,161,
		 59,250,255,255,255,255, 33, 27,155,185,248,255,255, 63,144,147, 27,146,178,249,255, 32,139,176,255,255,255,255, 63,178,255,255,255,255,255,255, 50, 40,168,138,249,255,255,159, 42,144,242,255,
		255,255,255, 50, 40,168, 16, 24,138,255, 31, 42,255,255,255,255,255,255, 49,152,129,255,255,255,255, 15, 25,255,255,255,255,255,255, 48,248,255,255,255,255,255,255,255,255,255,255,255,255,255
		};
	) + R(ushort edge_table(const uint i) {
		return edge_table_data[i < 128u ? i : 255u - i];
	}
	) + R(uchar triangle_table(const uint i) {
		return (triangle_table_data[i / 2u] >> (4u * (i % 2u))) & 0xF;
	}
	) + R(float3 interpolate_vertex(const float3 p1, const float3 p2, const float v1, const float v2, const float iso) { // linearly interpolate position where isosurface cuts an edge between 2 vertices
		const float w = (iso - v1) / (v2 - v1);
		return (1.0f - w) * p1 + w * p2;
	}
	) + R(uint marching_cubes(const float* v, const float iso, float3 * triangles) { // input: 8 values v, isovalue; output: returns number of triangles, 15 triangle vertices t
		uint cube = 0u; // determine index of which vertices are inside of the isosurface
		for (uint i = 0u; i < 8u; i++) cube |= (v[i] < iso) << i;
		if (cube == 0u || cube == 255u) return 0u; // cube is entirely inside/outside of the isosurface
		float3 p[8]; // definition of unit cube corners
		p[0] = (float3)(0.0f, 0.0f, 0.0f);
		p[1] = (float3)(1.0f, 0.0f, 0.0f);
		p[2] = (float3)(1.0f, 0.0f, 1.0f);
		p[3] = (float3)(0.0f, 0.0f, 1.0f);
		p[4] = (float3)(0.0f, 1.0f, 0.0f);
		p[5] = (float3)(1.0f, 1.0f, 0.0f);
		p[6] = (float3)(1.0f, 1.0f, 1.0f);
		p[7] = (float3)(0.0f, 1.0f, 1.0f);
		const uint edges = edge_table(cube);
		float3 vertex[12]; // find the vertices where the surface intersects the cube
		if (edges & 1u) vertex[0] = interpolate_vertex(p[0], p[1], v[0], v[1], iso); // calculate vertices on all 12 edges
		if (edges & 2u) vertex[1] = interpolate_vertex(p[1], p[2], v[1], v[2], iso);
		if (edges & 4u) vertex[2] = interpolate_vertex(p[2], p[3], v[2], v[3], iso);
		if (edges & 8u) vertex[3] = interpolate_vertex(p[3], p[0], v[3], v[0], iso);
		if (edges & 16u) vertex[4] = interpolate_vertex(p[4], p[5], v[4], v[5], iso);
		if (edges & 32u) vertex[5] = interpolate_vertex(p[5], p[6], v[5], v[6], iso);
		if (edges & 64u) vertex[6] = interpolate_vertex(p[6], p[7], v[6], v[7], iso);
		if (edges & 128u) vertex[7] = interpolate_vertex(p[7], p[4], v[7], v[4], iso);
		if (edges & 256u) vertex[8] = interpolate_vertex(p[0], p[4], v[0], v[4], iso);
		if (edges & 512u) vertex[9] = interpolate_vertex(p[1], p[5], v[1], v[5], iso);
		if (edges & 1024u) vertex[10] = interpolate_vertex(p[2], p[6], v[2], v[6], iso);
		if (edges & 2048u) vertex[11] = interpolate_vertex(p[3], p[7], v[3], v[7], iso);
		cube *= 15u;
		uint i; // number of triangle vertices
		for (i = 0u; i < 15u && triangle_table(cube + i) != 15u; i += 3u) { // create the triangles
			triangles[i] = vertex[triangle_table(cube + i)];
			triangles[i + 1u] = vertex[triangle_table(cube + i + 1u)];
			triangles[i + 2u] = vertex[triangle_table(cube + i + 2u)];
		}
		return i / 3u; // return number of triangles
	}
	) + R(typedef struct __attribute__((packed)) struct_ray {
		float3 origin;
		float3 direction;
	} ray;
	) + R(float intersect_sphere(const ray r, const float3 center, const float radius) {
		const float3 oc = center - r.origin;
		const float b = dot(oc, r.direction), c = sq(b) - dot(oc, oc) + sq(radius);
		return c < 0.0f ? -1.0f : b - sqrt(c);
	}
	) + R(float intersect_sphere_inside(const ray r, const float3 center, const float radius) {
		const float3 oc = center - r.origin;
		const float b = dot(oc, r.direction), c = sq(b) - dot(oc, oc) + sq(radius);
		return c < 0.0f ? -1.0f : b + sqrt(c);
	}
	) + R(float intersect_triangle(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
		const float3 u = p1 - p0, v = p2 - p0, w = r.origin - p0, h = cross(r.direction, v), q = cross(w, u);
		const float f = 1.0f / dot(u, h), s = f * dot(w, h), t = f * dot(r.direction, q);
		return (f < 0.0f || s < 0.0f || s>1.0f || t < 0.0f || s + t>1.0f) ? -1.0f : f * dot(v, q);
	}
	) + R(float intersect_triangle_bidirectional(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
		const float3 u = p1 - p0, v = p2 - p0, w = r.origin - p0, h = cross(r.direction, v), q = cross(w, u);
		const float f = 1.0f / dot(u, h), s = f * dot(w, h), t = f * dot(r.direction, q);
		return (s < 0.0f || s>1.0f || t < 0.0f || s + t>1.0f) ? -1.0f : f * dot(v, q);
	}
	) + R(float intersect_rhombus(const ray r, const float3 p0, const float3 p1, const float3 p2) { // Moeller-Trumbore algorithm
		const float3 u = p1 - p0, v = p2 - p0, w = r.origin - p0, h = cross(r.direction, v), q = cross(w, u);
		const float f = 1.0f / dot(u, h), s = f * dot(w, h), t = f * dot(r.direction, q);
		return (f < 0.0f || s < 0.0f || s>1.0f || t < 0.0f || t>1.0f) ? -1.0f : f * dot(v, q);
	}
	) + R(float intersect_plane(const ray r, const float3 p0, const float3 p1, const float3 p2) { // ray-triangle intersection, but skip barycentric coordinates
		const float3 u = p1 - p0, v = p2 - p0, w = r.origin - p0, h = cross(r.direction, v);
		const float f = 1.0f / dot(u, h);
		return f < 0.0f ? -1.0f : f * dot(v, cross(w, u));
	}
	) + R(float intersect_plane_always(const ray r, const float3 p0, const float3 p1, const float3 p2) { // ray-triangle intersection, but skip barycentric coordinates and visibility check
		const float3 u = p1 - p0, v = p2 - p0, w = r.origin - p0, h = cross(r.direction, v);
		return dot(v, cross(w, u)) / dot(u, h);
	}
	) + R(bool intersect_cuboid_bool(const ray r, const float3 center, const float Lx, const float Ly, const float Lz) {
		const float3 bmin = center - 0.5f * (float3)(Lx, Ly, Lz);
		const float3 bmax = center + 0.5f * (float3)(Lx, Ly, Lz);
		const float txa = (bmin.x - r.origin.x) / r.direction.x;
		const float txb = (bmax.x - r.origin.x) / r.direction.x;
		const float txmin = fmin(txa, txb);
		const float txmax = fmax(txa, txb);
		const float tya = (bmin.y - r.origin.y) / r.direction.y;
		const float tyb = (bmax.y - r.origin.y) / r.direction.y;
		const float tymin = fmin(tya, tyb);
		const float tymax = fmax(tya, tyb);
		if (txmin > tymax || tymin > txmax) return false;
		const float tza = (bmin.z - r.origin.z) / r.direction.z;
		const float tzb = (bmax.z - r.origin.z) / r.direction.z;
		const float tzmin = fmin(tza, tzb);
		const float tzmax = fmax(tza, tzb);
		return fmax(txmin, tymin) <= tzmax && tzmin <= fmin(txmax, tymax);
	}
	) + R(float intersect_cuboid(const ray r, const float3 center, const float Lx, const float Ly, const float Lz) {
		const float3 bmin = center - 0.5f * (float3)(Lx, Ly, Lz);
		const float3 bmax = center + 0.5f * (float3)(Lx, Ly, Lz);
		if (r.origin.x >= bmin.x && r.origin.y >= bmin.y && r.origin.z >= bmin.z && r.origin.x <= bmax.x && r.origin.y <= bmax.y && r.origin.z <= bmax.z) return 0.0f; // ray origin is within cuboid
		float3 p[8]; // 8 cuboid vertices
		p[0] = (float3)(bmin.x, bmin.y, bmin.z);
		p[1] = (float3)(bmax.x, bmin.y, bmin.z);
		p[2] = (float3)(bmax.x, bmin.y, bmax.z);
		p[3] = (float3)(bmin.x, bmin.y, bmax.z);
		p[4] = (float3)(bmin.x, bmax.y, bmin.z);
		p[5] = (float3)(bmax.x, bmax.y, bmin.z);
		p[6] = (float3)(bmax.x, bmax.y, bmax.z);
		p[7] = (float3)(bmin.x, bmax.y, bmax.z);
		float intersect = -1.0f;
		intersect = fmax(intersect, intersect_rhombus(r, p[0], p[3], p[4])); // test for intersections with the 6 cuboid faces
		intersect = fmax(intersect, intersect_rhombus(r, p[3], p[2], p[7])); // ray will intersect with either 0 or 1 rhombuses
		intersect = fmax(intersect, intersect_rhombus(r, p[2], p[1], p[6]));
		intersect = fmax(intersect, intersect_rhombus(r, p[1], p[0], p[5]));
		intersect = fmax(intersect, intersect_rhombus(r, p[7], p[6], p[4]));
		intersect = fmax(intersect, intersect_rhombus(r, p[1], p[2], p[0]));
		return intersect;
	}
	) + R(float3 reflect(const float3 direction, const float3 normal) {
		return direction - 2.0f * dot(direction, normal) * normal;
	}
	) + R(float3 refract(const float3 direction, const float3 normal, const float n) {
		const float direction_normal = dot(direction, normal);
		const float sqrt_part = sq(n) - 1.0f + sq(direction_normal);
		return sqrt_part >= 0.0f ? (direction - (direction_normal + sqrt(sqrt_part)) * normal) / n : direction - 2.0f * direction_normal * normal; // refraction : total internal reflection
	}
	) + R(ray get_camray(const int x, const int y, const float* camera_cache) {
		const float zoom = camera_cache[0]; // fetch camera parameters (rotation matrix, camera position, etc.)
		const float dis = camera_cache[1];
		const float posx = camera_cache[2];
		const float posy = camera_cache[3];
		const float posz = camera_cache[4];
		const float Rxx = camera_cache[5];
		const float Rxy = camera_cache[6];
		const float Rxz = camera_cache[7];
		const float Ryx = camera_cache[8];
		const float Ryy = camera_cache[9];
		const float Ryz = camera_cache[10];
		const float Rzx = camera_cache[11];
		const float Rzy = camera_cache[12];
		const float Rzz = camera_cache[13];
		const bool  vr = (as_int(camera_cache[14]) >> 31) & 0x1;
		const float rtv = (as_int(camera_cache[14]) >> 30) & 0x1 ? 2.0f : 1.0f;
		const float eye_distance = vload_half(28, (half*)camera_cache);
		const float stereo = (x < (int)def_screen_width / 2 ? -1.0f : 1.0f);
		float3 p0 = (float3)(!vr ? 0.0f : stereo * eye_distance / zoom, 0.0f, dis / zoom);
		float3 p1 = p0 + normalize((float3)(!vr ? (float)(x - (int)def_screen_width / 2) : ((float)(x - (int)def_screen_width / 2) - stereo * (float)(def_screen_width / 4u)) * rtv - stereo * eye_distance, (float)(y - (int)def_screen_height / 2), -dis));
		const float x0 = Rxx * p0.x + Ryx * p0.y + Rzx * p0.z; // reverse rotate p0
		const float y0 = Rxy * p0.x + Ryy * p0.y + Rzy * p0.z;
		const float z0 = Rxz * p0.x + Ryz * p0.y + Rzz * p0.z;
		const float x1 = Rxx * p1.x + Ryx * p1.y + Rzx * p1.z; // reverse rotate p1
		const float y1 = Rxy * p1.x + Ryy * p1.y + Rzy * p1.z;
		const float z1 = Rxz * p1.x + Ryz * p1.y + Rzz * p1.z;
		p0 = (float3)(x0, y0, z0) - (float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
		p1 = (float3)(x1, y1, z1) - (float3)(def_domain_offset_x, def_domain_offset_y, def_domain_offset_z);
		p0.x = p0.x + posx; p0.y = p0.y + posy; p0.z = p0.z + posz; // reverse transformation of p0
		p1.x = p1.x + posx; p1.y = p1.y + posy; p1.z = p1.z + posz; // reverse transformation of p1
		ray camray;
		camray.origin = p0;
		camray.direction = p1 - p0;
		return camray;
	}
	) + R(uint skybox_bottom(const ray r, const uint skybox_color) {
		const float3 p0 = (float3)(0.0f, 0.0f, -0.5f * (float)def_Nz), p1 = (float3)(1.0f, 0.0f, -0.5f * (float)def_Nz), p2 = (float3)(0.0f, 1.0f, -0.5f * (float)def_Nz);
		const float distance = intersect_plane(r, p0, p1, p2);
		if (distance > 0.0f) { // ray intersects with bottom
			const float3 normal = normalize(cross(p1 - p0, p2 - p0));
			float3 intersection = r.origin + distance * r.direction;
			const float scale = 2.0f / fmin((float)def_Nx, (float)def_Ny);
			int a = abs((int)floor(scale * intersection.x));
			int b = abs((int)floor(scale * intersection.y));
			const float r = scale * sqrt(sq(intersection.x) + sq(intersection.y));
			return color_mix((a % 2 == b % 2) * 0xFFFFFF, skybox_color, clamp(2.0f / r, 0.0f, 1.0f));
		}
		else {
			return skybox_color;
		}
	}
	) + R(uint skybox_color_bw(const float x, const float y) {
		return color_dim(0xFFFFFF, 1.0f - y);
	}
	) + R(uint skybox_color_hsv(const float x, const float y) {
		const float h = fmod(x * 360.0f + 120.0f, 360.0f);
		const float s = y > 0.5f ? 1.0f : 2.0f * y;
		const float v = y > 0.5f ? 2.0f - 2.0f * y : 1.0f;
		return hsv_to_rgb(h, s, v);
	}
	) + R(uint skybox_color_sunset(const float x, const float y) {
		return color_mix(255 << 16 | 175 << 8 | 55, y < 0.5f ? 55 << 16 | 111 << 8 | 255 : 0, 2.0f * (0.5f - fabs(y - 0.5f)));
	}
	) + R(uint skybox_color_grid(const float x, const float y) {
		int a = (int)(36.0f * x);
		int b = (int)(18.0f * y);
		return 0xFFFFFF * (a % 2 == b % 2);
	}
	) + R(uint skybox_color(const ray r, const global int* skybox) {
		//const float x = fma(atan2(r.direction.x, r.direction.y),  0.5f/3.1415927f, 0.5f);
		//const float y = fma(asin (r.direction.z               ), -1.0f/3.1415927f, 0.5f);
		//return color_mix(skybox_color_hsv(x, y), skybox_color_grid(x, y), 0.95f-0.33f*(2.0f*(0.5f-fabs(y-0.5f))));
		//return skybox_color_sunset(x, y);
		const float fu = (float)def_skybox_width * fma(atan2(r.direction.x, r.direction.y), 0.5f / 3.1415927f, 0.5f);
		const float fv = (float)def_skybox_height * fma(asin(r.direction.z), -1.0f / 3.1415927f, 0.5f);
		const int ua = clamp((int)fu, 0, (int)def_skybox_width - 1), va = clamp((int)fv, 0, (int)def_skybox_height - 1), ub = (ua + 1) % def_skybox_width, vb = min(va + 1, (int)def_skybox_height - 1); // bilinear interpolation positions
		const uint s00 = skybox[ua + va * def_skybox_width], s01 = skybox[ua + vb * def_skybox_width], s10 = skybox[ub + va * def_skybox_width], s11 = skybox[ub + vb * def_skybox_width];
		const float u1 = fu - (float)ua, v1 = fv - (float)va, u0 = 1.0f - u1, v0 = 1.0f - v1; // interpolation factors
		return color_mix(color_mix(s00, s01, v0), color_mix(s10, s11, v0), u0); // perform bilinear interpolation
	}
	) + R(uint last_ray_reflectivity(const ray reflection, const ray transmission, const uint last_color, const float reflectivity, const global int* skybox) {
		return color_mix(skybox_color(reflection, skybox), skybox_color(transmission, skybox), reflectivity);
	}
	) + R(float ray_grid_traverse(const ray r, const global float* phi, const global uchar * flags, float3 * normal, const uint Nx, const uint Ny, const uint Nz) {
		const float3 pa = r.origin;
		const float3 pb = r.origin + r.direction;
		const float xa = pa.x - 0.5f + 0.5f * (float)Nx, ya = pa.y - 0.5f + 0.5f * (float)Ny, za = pa.z - 0.5f + 0.5f * (float)Nz; // start point
		const float xb = pb.x - 0.5f + 0.5f * (float)Nx, yb = pb.y - 0.5f + 0.5f * (float)Ny, zb = pb.z - 0.5f + 0.5f * (float)Nz; // end point
		const int dx = (int)sign(xb - xa), dy = (int)sign(yb - ya), dz = (int)sign(zb - za); // fast ray-grid-traversal
		const float fxa = xa - floor(xa), fya = ya - floor(ya), fza = za - floor(za);
		int3 xyz = (int3)(floor(xa), floor(ya), floor(za));
		const float tdx = fmin((float)dx / (xb - xa), 1E7f);
		const float tdy = fmin((float)dy / (yb - ya), 1E7f);
		const float tdz = fmin((float)dz / (zb - za), 1E7f);
		float tmx = tdx * (dx > 0 ? 1.0f - fxa : fxa);
		float tmy = tdy * (dy > 0 ? 1.0f - fya : fya);
		float tmz = tdz * (dz > 0 ? 1.0f - fza : fza);
		while (true) {
			if (tmx < tmy) {
				if (tmx < tmz) { xyz.x += dx; tmx += tdx; }
				else { xyz.z += dz; tmz += tdz; }
			}
			else {
				if (tmy < tmz) { xyz.y += dy; tmy += tdy; }
				else { xyz.z += dz; tmz += tdz; }
			}
			if (xyz.x < -1 || xyz.y < -1 || xyz.z < -1 || xyz.x >= (int)Nx || xyz.y >= (int)Ny || xyz.z >= (int)Nz) break;
			else if (xyz.x < 0 || xyz.y < 0 || xyz.z < 0 || xyz.x >= (int)Nx - 1 || xyz.y >= (int)Ny - 1 || xyz.z >= (int)Nz - 1) continue;
			const uint x0 = (uint)xyz.x; // cube stencil
			const uint xp = (uint)xyz.x + 1u;
			const uint y0 = (uint)xyz.y * Nx;
			const uint yp = ((uint)xyz.y + 1u) * Nx;
			const uint z0 = (uint)xyz.z * Ny * Nx;
			const uint zp = ((uint)xyz.z + 1u) * Ny * Nx;
			uint j[8];
			j[0] = x0 + y0 + z0; // 000
			j[1] = xp + y0 + z0; // +00
			j[2] = xp + y0 + zp; // +0+
			j[3] = x0 + y0 + zp; // 00+
			j[4] = x0 + yp + z0; // 0+0
			j[5] = xp + yp + z0; // ++0
			j[6] = xp + yp + zp; // +++
			j[7] = x0 + yp + zp; // 0++
			uchar flags_cell = 0u; // check with cheap flags if the isosurface goes through the current marching-cubes cell (~15% performance boost)
			for (uint i = 0u; i < 8u; i++) flags_cell |= flags[j[i]];
			if (!(flags_cell & (TYPE_S | TYPE_E | TYPE_I))) continue; // cell is entirely inside/outside of the isosurface
			float v[8];
			for (uint i = 0u; i < 8u; i++) v[i] = phi[j[i]];
			float3 triangles[15]; // maximum of 5 triangles with 3 vertices each
			const uint tn = marching_cubes(v, 0.5f, triangles); // run marching cubes algorithm
			if (tn == 0u) continue; // if returned tn value is non-zero, iterate through triangles
			const float3 offset = (float3)((float)xyz.x + 0.5f - 0.5f * (float)Nx, (float)xyz.y + 0.5f - 0.5f * (float)Ny, (float)xyz.z + 0.5f - 0.5f * (float)Nz);
			for (uint i = 0u; i < tn; i++) {
				const float3 p0 = triangles[3u * i] + offset;
				const float3 p1 = triangles[3u * i + 1u] + offset;
				const float3 p2 = triangles[3u * i + 2u] + offset;
				const float intersect = intersect_triangle_bidirectional(r, p0, p1, p2); // for each triangle, check ray-triangle intersection
				if (intersect > 0.0f) { // intersection found (there can only be exactly 1 intersection)
					const uint xq = ((uint)xyz.x + 2u) % Nx; // central difference stencil on each cube corner point
					const uint xm = ((uint)xyz.x + Nx - 1u) % Nx;
					const uint yq = (((uint)xyz.y + 2u) % Ny) * Nx;
					const uint ym = (((uint)xyz.y + Ny - 1u) % Ny) * Nx;
					const uint zq = (((uint)xyz.z + 2u) % Nz) * Ny * Nx;
					const uint zm = (((uint)xyz.z + Nz - 1u) % Nz) * Ny * Nx;
					float3 n[8];
					n[0] = (float3)(phi[xm + y0 + z0] - v[1], phi[x0 + ym + z0] - v[4], phi[x0 + y0 + zm] - v[3]); // central difference stencil on each cube corner point
					n[1] = (float3)(v[0] - phi[xq + y0 + z0], phi[xp + ym + z0] - v[5], phi[xp + y0 + zm] - v[2]); // compute normal vectors from gradient
					n[2] = (float3)(v[3] - phi[xq + y0 + zp], phi[xp + ym + zp] - v[6], v[1] - phi[xp + y0 + zq]); // normalize later during trilinear interpolation more efficiently
					n[3] = (float3)(phi[xm + y0 + zp] - v[2], phi[x0 + ym + zp] - v[7], v[0] - phi[x0 + y0 + zq]);
					n[4] = (float3)(phi[xm + yp + z0] - v[5], v[0] - phi[x0 + yq + z0], phi[x0 + yp + zm] - v[7]);
					n[5] = (float3)(v[4] - phi[xq + yp + z0], v[1] - phi[xp + yq + z0], phi[xp + yp + zm] - v[6]);
					n[6] = (float3)(v[7] - phi[xq + yp + zp], v[2] - phi[xp + yq + zp], v[5] - phi[xp + yp + zq]);
					n[7] = (float3)(phi[xm + yp + zp] - v[6], v[3] - phi[x0 + yq + zp], v[4] - phi[x0 + yp + zq]);
					const float3 p = r.origin + intersect * r.direction - offset; // intersection point minus offset
					const float x1 = p.x - floor(p.x), y1 = p.y - floor(p.y), z1 = p.z - floor(p.z), x0 = 1.0f - x1, y0 = 1.0f - y1, z0 = 1.0f - z1; // calculate interpolation factors
					*normal = normalize(
						(x0 * y0 * z0 * rsqrt(fma(n[0].x, n[0].x, fma(n[0].y, n[0].y, fma(n[0].z, n[0].z, 1E-9f))))) * n[0] +
						(x1 * y0 * z0 * rsqrt(fma(n[1].x, n[1].x, fma(n[1].y, n[1].y, fma(n[1].z, n[1].z, 1E-9f))))) * n[1] +
						(x1 * y0 * z1 * rsqrt(fma(n[2].x, n[2].x, fma(n[2].y, n[2].y, fma(n[2].z, n[2].z, 1E-9f))))) * n[2] +
						(x0 * y0 * z1 * rsqrt(fma(n[3].x, n[3].x, fma(n[3].y, n[3].y, fma(n[3].z, n[3].z, 1E-9f))))) * n[3] +
						(x0 * y1 * z0 * rsqrt(fma(n[4].x, n[4].x, fma(n[4].y, n[4].y, fma(n[4].z, n[4].z, 1E-9f))))) * n[4] +
						(x1 * y1 * z0 * rsqrt(fma(n[5].x, n[5].x, fma(n[5].y, n[5].y, fma(n[5].z, n[5].z, 1E-9f))))) * n[5] +
						(x1 * y1 * z1 * rsqrt(fma(n[6].x, n[6].x, fma(n[6].y, n[6].y, fma(n[6].z, n[6].z, 1E-9f))))) * n[6] +
						(x0 * y1 * z1 * rsqrt(fma(n[7].x, n[7].x, fma(n[7].y, n[7].y, fma(n[7].z, n[7].z, 1E-9f))))) * n[7]
					); // perform normalization and trilinear interpolation
					return intersect; // intersection found, exit loop, process transmission ray
				}
			}
		}
		return -1.0f; // no intersection found
	}
	) + R(bool raytrace_phi_mirror(const ray ray_in, ray * ray_reflect, const global float* phi, const global uchar * flags, const global int* skybox, const uint Nx, const uint Ny, const uint Nz) { // only reflection
		float3 normal;
		float d = ray_grid_traverse(ray_in, phi, flags, &normal, Nx, Ny, Nz); // move ray through lattice, at each cell call marching_cubes
		if (d == -1.0f) return false; // no intersection found
		ray_reflect->origin = ray_in.origin + (d - 0.0003163f) * ray_in.direction; // start intersection points a bit in front triangle to avoid self-reflection
		ray_reflect->direction = reflect(ray_in.direction, normal);
		return true;
	}
	) + R(bool raytrace_phi(const ray ray_in, ray * ray_reflect, ray * ray_transmit, float* reflectivity, const global float* phi, const global uchar * flags, const global int* skybox, const uint Nx, const uint Ny, const uint Nz) {
		float3 normal;
		float d = ray_grid_traverse(ray_in, phi, flags, &normal, Nx, Ny, Nz); // move ray through lattice, at each cell call marching_cubes
		if (d == -1.0f) return false; // no intersection found
		ray_reflect->origin = ray_in.origin + (d - 0.0003163f) * ray_in.direction; // start intersection points a bit in front triangle to avoid self-reflection
		ray_reflect->direction = reflect(ray_in.direction, normal); // compute reflection ray
		ray ray_internal; // compute internal ray and transmission ray
		ray_internal.origin = ray_in.origin + (d + 0.0003163f) * ray_in.direction; // start intersection points a bit behind triangle to avoid self-transmission
		ray_internal.direction = refract(ray_in.direction, normal, def_n);
		const bool is_inside = dot(ray_in.direction, normal) > 0.0f; // camera is in fluid
		if (is_inside) { // swap ray_reflect and ray_internal
			const float3 ray_internal_origin = ray_internal.origin;
			ray_internal.origin = ray_reflect->origin; // start intersection points a bit in front triangle to avoid self-reflection
			ray_internal.direction = ray_reflect->direction;
			ray_reflect->origin = ray_internal_origin; // start intersection points a bit behind triangle to avoid self-transmission
			ray_reflect->direction = refract(ray_in.direction, -normal, 1.0f / def_n);
		}
		const float wr = sq(cb(2.0f * acospi(fabs(dot(ray_in.direction, normal))))); // increase reflectivity if ray intersects surface at shallow angle
		*reflectivity = clamp(is_inside ? 1.0f - wr : wr, 0.0f, 1.0f); // ray_reflect and ray_transmit are switched if camera is in fluid
		d = ray_grid_traverse(ray_internal, phi, flags, &normal, Nx, Ny, Nz);
		if (d != -1.0f) { // internal ray intersects isosurface
			const float3 intersection_point = ray_internal.origin + (d + 0.0003163f) * ray_internal.direction; // start intersection points a bit behind triangle to avoid self-transmission
			ray_transmit->origin = intersection_point;
			ray_transmit->direction = refract(ray_internal.direction, -normal, 1.0f / def_n);
		}
		else { // internal ray does not intersect again
			ray_transmit->origin = ray_internal.origin;
			ray_transmit->direction = ray_internal.direction;
		}
		return true;
	}
	) + R(bool is_above_plane(const float3 point, const float3 plane_p, const float3 plane_n) {
		return dot(point - plane_p, plane_n) >= 0.0f;
	}
	) + R(bool is_below_plane(const float3 point, const float3 plane_p, const float3 plane_n) {
		return dot(point - plane_p, plane_n) <= 0.0f;
	}
	) + R(bool is_in_camera_frustrum(const float3 p, const float* camera_cache) { // returns true if point is located in camera frustrum
		const float vr = (float)((as_int(camera_cache[14]) >> 31) & 0x1);
		const ray r00 = get_camray(0, 0, camera_cache); // get 4 edge vectors of frustrum
		const ray r01 = get_camray((int)def_screen_width - 1, 0, camera_cache);
		const ray r10 = get_camray(0, (int)def_screen_height - 1, camera_cache);
		const ray r11 = get_camray((int)def_screen_width - 1, (int)def_screen_height - 1, camera_cache);
		const float3 plane_n_top = cross(r00.direction, r01.direction); // get 4 frustrum planes
		const float3 plane_n_bottom = cross(r11.direction, r10.direction);
		const float3 plane_n_left = cross(r10.direction, r00.direction);
		const float3 plane_n_right = cross(r01.direction, r11.direction);
		const float3 plane_p_top = r00.origin - 2.0f * plane_n_top; // move frustrum planes outward by 2 units
		const float3 plane_p_bottom = r11.origin - 2.0f * plane_n_bottom;
		const float3 plane_p_left = r00.origin - (2.0f + 16.0f * vr) * plane_n_left; // move frustrum planes outward by 2 units, for stereoscopic rendering a bit more
		const float3 plane_p_right = r11.origin - (2.0f + 16.0f * vr) * plane_n_right;
		return is_above_plane(p, plane_p_top, plane_n_top) && is_above_plane(p, plane_p_bottom, plane_n_bottom) && is_above_plane(p, plane_p_left, plane_n_left) && is_above_plane(p, plane_p_right, plane_n_right);
	}
	) + "#endif" + R( // GRAPHICS
	);
}