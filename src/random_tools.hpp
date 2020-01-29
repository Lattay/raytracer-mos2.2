#ifndef RANDOM_TOOLS_H
#define RANDOM_TOOLS_H

double roll();
typedef struct { double x; double y; } Vec2;
Vec2 box_muller(double sigma);

#endif
