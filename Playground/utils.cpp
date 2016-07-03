#include "utils.h"
#include <math.h>

int ImageLoad(char *filename, Image *image) {
    FILE *file;
    unsigned long size;                 // size of the image in bytes.
    unsigned long i;                    // standard counter.
    unsigned short int planes;          // number of planes in image (must be 1)
    unsigned short int bpp;             // number of bits per pixel (must be 24)
    char temp;                          // temporary color storage for bgr-rgb conversion.
    if ((file = fopen(filename, "rb"))==NULL)    {
    	printf("File Not Found : %s\n",filename);
	    return 0;
    }
    fseek(file, 18, SEEK_CUR);
    if ((i = fread(&image->sizeX, 4, 1, file)) != 1) {
     	printf("Error reading width from %s.\n", filename);
	    return 0;
    }
    printf("Width of %s: %lu\n", filename, image->sizeX);
    if ((i = fread(&image->sizeY, 4, 1, file)) != 1) {
	    printf("Error reading height from %s.\n", filename);
	    return 0;
    }
    printf("Height of %s: %lu\n", filename, image->sizeY);
    size = image->sizeX * image->sizeY * 3;
    if ((fread(&planes, 2, 1, file)) != 1) {
	    printf("Error reading planes from %s.\n", filename);
	    return 0;
    }
    if (planes != 1) {
	     printf("Planes from %s is not 1: %u\n", filename, planes);
	     return 0;
    }
    if ((i = fread(&bpp, 2, 1, file)) != 1) {
	     printf("Error reading bpp from %s.\n", filename);
	     return 0;
    }
    if (bpp != 24) {
	      printf("Bpp from %s is not 24: %u\n", filename, bpp);
	      return 0;
    }
    fseek(file, 24, SEEK_CUR);
    image->data = (char *) malloc(size);
    if (image->data == NULL) {
	      printf("Error allocating memory for color-corrected image data");
      	return 0;
    }
    if ((i = fread(image->data, size, 1, file)) != 1) {
	     printf("Error reading image data from %s.\n", filename);
	     return 0;
    }
    for (i=0;i<size;i+=3) { // reverse all of the colors. (bgr -> rgb)
      	temp = image->data[i];
	      image->data[i] = image->data[i+2];
	      image->data[i+2] = temp;
    }
    return 1;
}

float translate(float value,
                float leftMin, float leftMax,
                float rightMin, float rightMax)
{
    float leftSpan = leftMax - leftMin;
    float rightSpan = rightMax - rightMin;
    float valueScaled = float(value - leftMin) / float(leftSpan);
    return rightMin + (valueScaled * rightSpan);
}

