#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "vector.h"

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

Vector random_direction() {
    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));
    double x = cos(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double y = sin(2 * M_PI * r1) * sqrt(r2 * (1 - r2));
    double z = 1 - 2 * r2;
    return Vector(x, y, z);
}

void sort_projection(std::vector<std::pair<int, int>> proj) {
    std::sort(proj.begin(), proj.end(),
        [](const std::pair<int, int> &a, std::pair<int, int> const &b) -> bool {
            return a.first > b.first;
        }
    );
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("Error: first argument should be input image, second should be model image, and third should be iterations.\n");
        exit(1);
    }
    const char *input_image_filename = argv[1];
    const char *model_image_filename = argv[2];
    size_t iterations = atoi(argv[3]);

    int input_W, input_H, input_channels;
    unsigned char *input_image = stbi_load(input_image_filename, &input_W, &input_H, &input_channels, 0);
    if (input_image == NULL) {
        printf("Error: failed to load the input image.\n");
        exit(1);
    }

    int model_W, model_H, model_channels;
    unsigned char *model_image = stbi_load(model_image_filename, &model_W, &model_H, &model_channels, 0);
    if (model_image == NULL) {
        printf("Error: failed to load the model image.\n");
        exit(1);
    }

    if (model_W != input_W || model_H != input_H || model_channels != input_channels) {
        printf("Error: input and model must have same dimensions and channels.\n");
        exit(1);
    }

    size_t total_pixels = input_W * input_H;
    std::vector<std::pair<int, int>> projI(total_pixels);
    std::vector<std::pair<int, int>> projM(total_pixels);
    Vector pixel, model_pixel, v;
    
    for (size_t iter = 0; iter < iterations; iter++) {
        v = random_direction();

        // Projections
        for (size_t i = 0; i < total_pixels; i++) {
            unsigned char *I = input_image + input_channels * i;
            unsigned char *M = model_image + model_channels * i;
            pixel = Vector(*I, *(I + 1), *(I + 2));
            model_pixel = Vector(*M, *(M + 1), *(M + 2));
            projI[i] = std::pair<int, int>(dot(pixel, v), i);
            projM[i] = std::pair<int, int>(dot(model_pixel, v), i);
        }

        // Sorting
        sort_projection(projI);
        sort_projection(projM);

        // Advection
        for (size_t i = 0; i < total_pixels; i++) {
            int permutation_index = projI[i].second;
            unsigned char *I = input_image + input_channels * permutation_index;
            pixel = Vector(*I, *(I + 1), *(I + 2)) + (projM[i].first - projI[i].first)*v;
            *I = pixel[0];
            *(I + 1) = pixel[1];
            *(I + 2) = pixel[2];
        }
    }

    stbi_write_png("output.png", input_W, input_H, input_channels, &input_image[0], 0);
    
    return 0;
}