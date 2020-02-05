#ifndef RANDOM_TOOLS_H
#define RANDOM_TOOLS_H
#include "Vec.hpp"

const double pi = 3.141592654;

typedef struct { Vec dir; double proba; } Sample;

double roll();
typedef struct { double x; double y; } Vec2;
Vec2 box_muller(double sigma);
Base ortho_prod(Vec const& v);
Sample random_cos(Vec const& n);

#endif
