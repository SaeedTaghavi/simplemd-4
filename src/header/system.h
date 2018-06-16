#ifndef SYSTEM
#define SYSTEM

#include "atom.h"

// Constants to determine seeding type

class System {
private:
    double size; // Box side length
    int N; // Number of particles
    double dt; // Length of the timestep
    double Ek; // Kinetic energy
    double Ep; // Potential energy
    double rc; // Cut-off limit
    double dx;
    double dy;
    int div;
    Atom *first; // Pointer to first Atom
    Atom *current; // Pointer to current Atom
protected:
public:

    // To handle the tree structure

    System() : size(0), N(0), dt(0), Ek(0), Ep(0), rc(0), first(0), current(0), div(0) {
    };
    ~System();
    Atom *AddAtom(double x, double y, double vx, double vy, double fx, double fy, int tag);
    void KillAll();

    // Misc routines for handling the MD system

    void ForceLoop();
    void Initialize(int shape, double kinetic, double cutoff, double timestep, int number, double side, int seed);
    void InitializeFromFile(char *fname);
    void Distance(int k, int l);
    void Kinetic();
    void RemoveLinearMomentum();
    void ScaleKineticEnergy(double dEk);
    void IntegrateVerlet();
    void IntegrateEuler();
    void IntegrateGearPredictor();
    void IntegratePartitionRungeKutta();

    // File handling routines

    void ToNewFile(char *fname);
    void AppendToFile(char *fname);
    void VitalsToFile(char *fname);
    void AppendToTinker(char *fname);
    void NewTinker(char *fname);
    void NewVitalsToFile(char *fname);
    void NewSingleToFile(char *fname);
    void SingleToFile(char *fname);

    // Get data from structure

    int getN();
    double getEk();
    double getEp();
    double getLm();
    double T(); // Return temperature

};

#endif
