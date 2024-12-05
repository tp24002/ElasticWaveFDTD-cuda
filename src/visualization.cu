#include <stdio.h>
#include <stdlib.h>
#include <png.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "../header/struct.h"
#include "../header/visualization.h"

void ensure_directory_exists(const char* path) {
    struct stat st = {0};
    if (stat(path, &st) == -1) {
        mkdir(path, 0700);
    }
}

void write_png(const char* filename, int width, int height, png_bytep image) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "Failed to open file %s for writing\n", filename);
        exit(EXIT_FAILURE);
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) {
        fprintf(stderr, "Failed to create PNG write struct\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    png_infop info = png_create_info_struct(png);
    if (!info) {
        fprintf(stderr, "Failed to create PNG info struct\n");
        png_destroy_write_struct(&png, NULL);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    if (setjmp(png_jmpbuf(png))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png, &info);
        fclose(fp);
        exit(EXIT_FAILURE);
    }

    png_init_io(png, fp);

    png_set_IHDR(
        png, info, width, height,
        8, PNG_COLOR_TYPE_RGBA,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE
    );

    png_write_info(png, info);

    png_bytep *row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = image + y * width * 4;
    }
    png_write_image(png, row_pointers);
    png_write_end(png, NULL);

    free(row_pointers);
    fclose(fp);
    png_destroy_write_struct(&png, &info);
}

void create_png(DimI3 ran, double *STV, int cut, int t) {
    ensure_directory_exists("./createpng");
    int width, height;
    width  = ran.x;
    height = ran.y;  
    
    png_bytep image = (png_bytep)malloc(width * height * 4); // RGBA

    // 振幅の最大・最小値
    double max =  1.e-0;
    double min = -1.e-0;

    // 色の定義（カラーマップ: 青→水色→緑→黄→赤）
    int cmap_size = 5;
    int cmap[5][3] = {
        {0, 0, 255},    // 青
        {0, 255, 255},  // 水色
        {0, 255, 0},    // 緑
        {255, 255, 0},  // 黄
        {255, 0, 0}     // 赤
    };

    // 背景を白で初期化
    for (int i = 0; i < width * height * 4; i++) {
        image[i] = 255;
    }

    // ピクセルごとに計算
    
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            int id = cut * width * height + j * width + i;
            int idimage = (j * width + i) * 4;

            // 振幅の正規化
            double normalized = (STV[id] - min) / (max - min);
            if (normalized < 0) normalized = 0;
            if (normalized > 1) normalized = 1;

            // カラーマップの選択
            int idx = (int)(normalized * (cmap_size - 1));
            double frac = (normalized * (cmap_size - 1)) - idx;

            // 色を線形補間
            int r = cmap[idx][0] + frac * (cmap[idx + 1][0] - cmap[idx][0]);
            int g = cmap[idx][1] + frac * (cmap[idx + 1][1] - cmap[idx][1]);
            int b = cmap[idx][2] + frac * (cmap[idx + 1][2] - cmap[idx][2]);

            // ピクセルに色を設定
            image[idimage]     = r;
            image[idimage + 1] = g;
            image[idimage + 2] = b;
            image[idimage + 3] = 255; // アルファ値（透明度）
        }
    }

    // ファイル名と保存
    char filename[256];
    sprintf(filename, "./createpng/frame_%05d.png", t);
    write_png(filename, width, height, image);

    free(image);
}

