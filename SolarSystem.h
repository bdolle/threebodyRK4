//
// Created by Brian Dolle on 11/1/18.
//

#ifndef THREEBODY_SOLARSYSTEM_H
#define THREEBODY_SOLARSYSTEM_H

#include "Coords.h"

class SolarSystem {

public:
    SolarSystem(double M1, double M2, double M3); // constructor
    ~SolarSystem();

    double SunVxDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);
    double SunVyDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);
    double SunVzDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);

    double PlanetVxDot(double xs, double ys, double zs, double xp, double yp, double zp,
                       double xsat, double ysat, double zsat,
                       double Ms, double Mp, double Msat);
    double PlanetVyDot(double xs, double ys, double zs, double xp, double yp, double zp,
                       double xsat, double ysat, double zsat,
                       double Ms, double Mp, double Msat);
    double PlanetVzDot(double xs, double ys, double zs, double xp, double yp, double zp,
                       double xsat, double ysat, double zsat,
                       double Ms, double Mp, double Msat);

    double SatVxDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);
    double SatVyDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);
    double SatVzDot(double xs, double ys, double zs, double xp, double yp, double zp,
                    double xsat, double ysat, double zsat,
                    double Ms, double Mp, double Msat);




    void runsystem();

    Coords sun_;
    Coords planet_;
    Coords sat_;

    double t_, h_;


};



#endif //THREEBODY_SOLARSYSTEM_H
