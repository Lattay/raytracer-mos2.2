#include "src/main.hpp" 
 
 
int main() {
    int W = 512;
    int H = 512;
   
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
 
            image[(i*W + j) * 3 + 0] = 255;
            image[(i*W + j) * 3 + 1] = 0;
            image[(i*W + j) * 3 + 2] = 0;
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}
