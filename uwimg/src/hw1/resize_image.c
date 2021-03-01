#include <math.h>
#include "image.h"

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO Fill in
	return get_pixel(im, round(x), round(y),c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix that first line)
	
	image im_nn = make_image(w,h,im.c);

	int i, j, k;
	float w_mapped = 0;
	float h_mapped = 0; 
	float v = 0;

	float w_scale = (float) im.w/w;
	float h_scale = (float) im.h/h;
	for(i = 0 ; i < im.c ; i++){

		for(j = 0 ; j < h ; j++){

			for(k = 0 ; k < w; k++){

				w_mapped = w_scale*(k+0.5) - 0.5;
				h_mapped = h_scale*(j+0.5) - 0.5; 
				v = nn_interpolate(im, w_mapped, h_mapped, i);
				set_pixel(im_nn, k, j, i, v);

			}
		}
	}

    return im_nn;
}

int cap_index(int index, int size) {
  return index < size - 1 ? index >= 0 ? index : 0 : size - 1;
}

float to_original_scale(int new_index, int original_size, int new_size) {
  return ((new_index + 0.5) * (original_size) / (new_size)) - 0.5;
}

void find_bound(int *bound, float current_pos, int size) {
  bound[0] = cap_index((int)current_pos, size);
  bound[1] = cap_index((int)(current_pos + 1), size);
}

void find_distances(float *dis, float current_pos, int *bound) {
  if (bound[0] == bound[1]) {
    dis[0] = 0.5;
    dis[1] = 0.5;
  } else {
    dis[0] = current_pos - bound[0];
    dis[1] = bound[1] - current_pos;
  }
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    
	int rows[2], cols[2], row, col;
  float d_rows[2], d_cols[2], color_value = 0;

  find_bound(rows, y, im.h);
  find_bound(cols, x, im.w);

  find_distances(d_rows, y, rows);
  find_distances(d_cols, x, cols);

  for (row = 0; row < 2; row++) {
    for (col = 0; col < 2; col++) {
      color_value += get_pixel(im, cols[col], rows[row], c) * d_rows[1 - row] *
                     d_cols[1 - col];
    }
  }
  return color_value;

}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    image im_bilinear = make_image(w,h,im.c);

	int i, j, k;
	float w_mapped = 0;
	float h_mapped = 0; 
	float v = 0;
	float w_scale = (float) im.w/w;
	float h_scale = (float) im.h/h;

	for(i = 0 ; i < im.c ; i++){

		for(j = 0 ; j < h ; j++){

			for(k = 0 ; k < w; k++){

				w_mapped = w_scale*(k+0.5) - 0.5;
				h_mapped = h_scale*(j+0.5) - 0.5; 
				v = bilinear_interpolate(im, w_mapped, h_mapped, i);
				set_pixel(im_bilinear, k, j, i, v);

			}
		}
	}

    return im_bilinear;

}

