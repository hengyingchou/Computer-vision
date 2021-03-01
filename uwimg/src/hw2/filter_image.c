#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    // TODO
	int size = im.h*im.h*im.c;
	float sum = 0.0;

	for(int i = 0 ; i < size ; i++){

		sum += im.data[i];

	}

	for(int i = 0 ; i < size ; i++){

		im.data[i] /= sum;

	}

}

image make_box_filter(int w)
{
    // TODO

	image filter = make_image(w,w,1);
	float value = 1.0/(w*w);

	float size = w*w;
	for(int i = 0 ; i < size ; i++){
		filter.data[i] = value;
	}

	return filter;
}

image image_expand(image im, int W, int H){


	
	int i, j, k;
	float value = 0.0;
	image expand = make_image(im.w + 2*W, im.h + 2*H, im.c);

	for(i = 0; i < expand.c; i++){
		for(j = 0; j < expand.h; j++){
			for(k = 0; k < expand.w ; k++){
				set_pixel(expand, k , j, i, 0);
			}
		}
	}

	for(i = 0 ; i < im.c ; i++){
		for(j = 0; j < im.h; j++){
			for(k = 0 ; k < im.w ; k++){
				value = get_pixel(im, k, j, i);
				set_pixel(expand, k+W, j+H , i, value);
			}
		}
	}


	int end_w = expand.w - W - 1;
	


	for(i = 0 ; i < im.c ; i++){
		for(j = 0 ; j < expand.h ; j++){
			for( k = 0 ; k < W ; k++){
				value = get_pixel(expand, W, j,i);
				set_pixel(expand, k, j, i, value);
				value = get_pixel(expand, end_w, j, i);
				set_pixel(expand, k + end_w + 1, j , i , value);
			}
		}
	}

	int end_h = expand.h - H - 1;
	for(i = 0 ; i < im.c ; i++){
		for(j = 0 ; j < expand.w ; j++){
			for(k = 0 ; k < H ; k++){

				value = get_pixel(expand, j, H,i);
				set_pixel(expand, j, k, i, value);
				value = get_pixel(expand, j, end_h , i);
				set_pixel(expand, j, k + end_h + 1 , i , value);

			}
		}
	} 

	return expand;

}

float convolv_process(image im, int col, int row, int c, image filter, int cf){


	int i , j ;
	float value = 0.0;

		for(i = 0; i < filter.h ; i++){
			for(j = 0 ; j < filter.w ; j++){

				value += get_pixel(im, col + j, row + i, c)*get_pixel(filter, j, i ,cf);

			}
		}
	
	return value;

}


image convolve_image(image im, image filter, int preserve)
{
    // TODO

	image object = make_image(im.w , im.h , preserve? im.c : 1); 
	int W = (filter.w - 1)/2;
	int H = (filter.h - 1)/2;

	im = image_expand(im,W,H); 

	int i, j, k;
	int cf;
	float value = 0.0;
	float org = 0.0;
	for(i = 0 ; i < im.c ; i++){
		cf = filter.c == 1? 0: i;
		for(j = 0 ; j < object.h ; j++){
			for(k = 0 ; k < object.w ; k++){
				value = convolv_process(im, k, j, i, filter,cf);
				if(preserve){
					set_pixel(object, k , j , i , value);

				}else{
					org = get_pixel(object, k, j, 0);
					//printf("good\n");
					set_pixel(object, k, j, 0, org+value);
					
				}

			}
		}
	}

	free_image(im);
    return object;
}

image make_highpass_filter()
{
    // TODO
     
	int Hight[9] = {0, -1, 0, -1, 4, -1, 0, -1, 0}; 
	image filter = make_image(3,3,1);
	for(int i = 0 ; i < 9 ; i++){
		filter.data[i] = Hight[i];
	}
	return filter;

    //return make_image(1,1,1);
}

image make_sharpen_filter()
{
    // TODO
    int Sharp[9] = {0, -1, 0, -1, 5, -1, 0, -1, 0}; 
	image filter = make_image(3,3,1);
	for(int i = 0 ; i < 9 ; i++){
		filter.data[i] = Sharp[i];
	}
	return filter;
    
}

image make_emboss_filter()
{
    // TODO
    int Emboss[9] = {-2, -1, 0, -1, 1, 1, 0, 1, 2}; 
	image filter = make_image(3,3,1);
	for(int i = 0 ; i < 9 ; i++){
		filter.data[i] = Emboss[i];
	}
	return filter;

}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO





// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO


image make_gaussian_filter(float sigma){

	int width = (int) ceilf(6*sigma);
	width = width - (width%2) + 1;

	int mid = (width - 1)/2;
	image filter = make_image(width, width,1);

	int i , j ;
	float value;
	float sigma2;

	for(i = 0 ; i < width ; i++){
		for(j = 0 ; j < width ; j++){

			sigma2 = 2*sigma*sigma;
			value = (exp(-((j-mid)*(j-mid) + (i-mid)*(i-mid))/sigma2))/(M_PI*sigma2);

			set_pixel(filter, j, i, 0, value);
		}
	}

	l1_normalize(filter);
	return filter;



}



image add_image(image a, image b)
{
    // TODO
	
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
  image result = make_image(a.w, a.h, a.c);
  int result_size = result.w * result.h * result.c;
  for (int i = 0; i < result_size; i++) {
    result.data[i] = a.data[i] + b.data[i];
  }
  return result;
    //return make_image(1,1,1);
}

image sub_image(image a, image b)
{
    // TODO

	assert(a.w == b.w && a.h == b.h && a.c == b.c);
	image result = make_image(a.w,a.h,a.c);
	int size = a.w*a.h*a.c;
	for(int i = 0 ; i < size; i++){

		result.data[i] = a.data[i] - b.data[i];

	}
	return result;
    //return make_image(1,1,1);
}

image make_gx_filter()
{
    // TODO
	int gx[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1}; 
	image filter = make_image(3,3,1);
	for(int i = 0 ; i < 9 ; i++){
		filter.data[i] = gx[i];
	}
	return filter;

}

image make_gy_filter()
{
    // TODO
	int gy[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1}; 
	image filter = make_image(3,3,1);
	for(int i = 0 ; i < 9 ; i++){
		filter.data[i] = gy[i];
	}
	return filter;
}

void feature_normalize(image im)
{
    // TODO
	float min = im.data[0];
	float max = im.data[0];
	int size = im.w*im.h*im.c;
	int i ;

	for( i = 0 ; i < size ; i++){

		if(im.data[i] > max) max = im.data[i];
		if(im.data[i] < min) min = im.data[i];

	}
	float range = max - min;
	for(i = 0 ; i < size ; i++){

		im.data[i] = range==0 ? 0 : (im.data[i] - min)/range;
	}

}

image *sobel_image(image im)
{
    // TODO
	image gx_filter = make_gx_filter();
	image gy_filter = make_gy_filter();
	image gx = convolve_image(im, gx_filter, 0);
	image gy = convolve_image(im, gy_filter, 0);
	image *result = calloc(2, sizeof(image));
	image g = make_image(im.w, im.h, 1);
	image theta = make_image(im.w, im.h, 1);

	int size = im.w*im.h*1;
	int i ;
	float x,y;

	for(i = 0 ; i < size ; i++){
		x = gx.data[i];
		y = gy.data[i]; 
		g.data[i] = sqrt(x*x + y*y);
		theta.data[i] = atan2f(y,x);
	}

	result[0] = g;
	result[1] =theta; 

	free_image(gx_filter);
	free_image(gy_filter);
	free_image(gx);
	free_image(gy);


    return result;
}

image colorize_sobel(image im)
{
    // TODO
	int channel = 3;
	image *sobel = sobel_image(im);
	image g = sobel[0];
	image theta = sobel[1];
	image result = make_image(im.w, im.h, channel);

	feature_normalize(g);
	feature_normalize(theta);

	float g_data;
	float theta_data;
	int i, j;

	for(i = 0 ; i < im.h; i++){
		for(j = 0 ; j < im.w ; j++){
			g_data = get_pixel(g, j, i, 0);
			theta_data = get_pixel(theta, j, i, 0);

			set_pixel(result, j , i, 0, theta_data);
			set_pixel(result, j , i, 1, g_data);
			set_pixel(result, j , i, 2, g_data);
		}
	}
	hsv_to_rgb(result);

    return result;
}
