#include <iostream>
#include "SolarSystem.h"

int main() {
    std::cout << "Hello, World!" << std::endl;

    double Msun = 1;
    double Mplanet = 1e-5;
    double Msat = 1e-10;
    SolarSystem threebody(Msun, Mplanet, Msat);
    threebody.runsystem();
    return 0;
}