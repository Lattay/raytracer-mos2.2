#include <random>
#include <cmath>
#include <iostream>

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0, 1);

const double pi = 3.141592654;

static double roll(){
  return distribution(generator);
}

static double box_muller(double sigma){
  double u1 = roll(), u2 = roll();
  return sigma * sqrt(-2 * log(u1)) * cos(2 * pi * u2);
}

static double gauss(double sigma, double x){

  return exp(-0.5 * x*x / (sigma * sigma)) / (sigma * sqrt(2 * pi));
}

static double cos30(int n){
  double s = 0;
  double sigma = 0.25;
  for(int i = 0; i < n; i++){
    double x = box_muller(sigma);
    s += pow(cos(x), 30) / gauss(sigma, x);
  }
  return s/(double)n;
}

static double cos30xyz(int n){
  double s = 0;
  double sigma = 0.23;
  for(int i = 0; i < n; i++){
    double u1 = roll(), u2 = roll();
    double x = sigma * sqrt(-2 * log(u1)) * cos(2 * pi * u2);
    double y = sigma * sqrt(-2 * log(u1)) * sin(2 * pi * u2);
    u1 = roll(), u2 = roll();
    double z = sigma * sqrt(-2 * log(u1)) * cos(2 * pi * u2);
    s += pow(cos(x*y*z), 30) / gauss(sigma, x) / gauss(sigma, y) / gauss(sigma, z);
  }
  return s/(double)n;
}

int main(){
  std::cout << cos30(1000000) << std::endl;
  std::cout << cos30xyz(1000000) << std::endl;
  return 0;
}
