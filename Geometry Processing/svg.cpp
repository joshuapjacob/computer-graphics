#include "svg.h"
#include <stdio.h>
#include <sstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

// POLYGON ---------------------------------------------------------------------

double Polygon::area() {
    size_t n = this->vertices.size();
    if (n < 3) return 0.0;
    double tmp = 0.0;
    for (uint i = 0; i < n - 1; i++) {
        tmp += vertices[i][0]*vertices[i + 1][1] - vertices[i][1]*vertices[i + 1][0];
    }
    tmp += vertices[n - 1][0]*vertices[0][1] - vertices[n - 1][1]*vertices[0][0];
    return std::abs(tmp/2);
};

Vector Polygon::centroid() {
    double Cx = 0.0;
    double Cy = 0.0;
    size_t n = vertices.size();
    for (uint i = 0; i < n - 1; i++) {
        Cx += (vertices[i][0] + vertices[i+1][0]) *
              (vertices[i][0]*vertices[i+1][1] -
               vertices[i+1][0]*vertices[i][1]);
        Cy += (vertices[i][1] + vertices[i+1][1]) *
              (vertices[i][0]*vertices[i+1][1] -
               vertices[i+1][0]*vertices[i][1]);
    }
    Cx += (vertices[n-1][0] + vertices[0][0]) *
            (vertices[n-1][0]*vertices[0][1] -
            vertices[0][0]*vertices[n-1][1]);
    Cy += (vertices[n-1][1] + vertices[0][1]) *
            (vertices[n-1][0]*vertices[0][1] -
            vertices[0][0]*vertices[n-1][1]);
    Vector centroid(-Cx, -Cy, 0.0);
    return centroid/(this->area()*6.);
}

// SVG -------------------------------------------------------------------------

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol) {
    FILE* f = fopen(filename.c_str(), "w+");
    fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
    for (uint i=0; i<polygons.size(); i++) {
        fprintf(f, "<g>\n");
        fprintf(f, "<polygon points = \"");
        for (uint j = 0; j < polygons[i].vertices.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
        fprintf(f, "</g>\n");
    }
    fprintf(f, "</svg>\n");
    fclose(f);
}

void save_frame(
        const std::vector<Polygon> &cells,
        std::string filename, int frameid,
        uint N
    ) {
    int W = 1000, H = 1000;
    std::vector<unsigned char> image(W*H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {

        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].vertices.size(); j++) {
            bminx = std::min(bminx, cells[i].vertices[j][0]);
            bminy = std::min(bminy, cells[i].vertices[j][1]);
            bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
            bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
        }
        bminx = std::min(W-1., std::max(0., W * bminx));
        bminy = std::min(H-1., std::max(0., H * bminy));
        bmaxx = std::max(W-1., std::max(0., W * bmaxx));
        bmaxy = std::max(H-1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].vertices.size(); j++) {
                    double x0 = cells[i].vertices[j][0] * W;
                    double y0 = cells[i].vertices[j][1] * H;
                    double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
                    double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
                    double det = (x - x0)*(y1-y0) - (y - y0)*(x1-x0);
                    int sign = (det > 0) - (det < 0);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                        if (sign != prevSign) {
                            isInside = false;
                            break;
                        }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
                    double distEdge = std::abs(det)/ edgeLen;
                    double dotp = (x - x0)*(x1 - x0) + (y - y0)*(y1 - y0);
                    if (dotp<0 || dotp>edgeLen*edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    if (i < N) { // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2) {
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 0;
                        image[((H - y - 1)*W + x) * 3 + 2] = 0;
                    }
                }
            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}