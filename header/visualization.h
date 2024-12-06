#pragma once
#include <png.h>

#include "../header/struct.h"

void ensure_directory_exists(const char* path);
void write_png(const char* filename, int width, int height, png_bytep image);
void create_png(DimI3 ran, double *STV, int cut, int t);
// void create_png(Range ran, int cut, int t);