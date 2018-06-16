#include "header/system.h"
#include "header/atom.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>

void System::NewSingleToFile(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::trunc);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(7);

    ofile << temp->getx() << "\t" << temp->gety() << std::endl;

    ofile.close();

}

void System::SingleToFile(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::app);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(7);

    ofile << temp->getx() << "\t" << temp->gety() << std::endl;

    ofile.close();

}

double System::getLm() {
    double lmx = 0;
    double lmy = 0;
    Atom *temp = current;
    while (temp != 0) {
        lmx = lmx + temp->getvx();
        lmy = lmy + temp->getvy();
        temp = temp->getnext();
    }
    return (sqrt(lmx * lmx + lmy * lmy));
}

void System::IntegrateEuler() {
    ForceLoop();
    Atom *temp = current;
    while (temp != 0) {
        temp->setx(temp->getx() + temp->getvx() * dt);
        temp->sety(temp->gety() + temp->getvy() * dt);
        temp->setvx(temp->getvx() + temp->getfx() * dt);
        temp->setvy(temp->getvy() + temp->getfy() * dt);
        temp->setx(temp->getx() - size * floor(temp->getx() / size));
        temp->sety(temp->gety() - size * floor(temp->gety() / size));
        temp = temp->getnext();
    }
    Kinetic();
}

void System::IntegrateGearPredictor() {

}

void System::IntegratePartitionRungeKutta() {

    Atom *temp = current;


    // Coefficients for the Runge-Kutta
    int s = 7;
    double B[7];
    double C[7];
    B[0] = 0.0792036964311957;
    B[1] = 0.353172906049774;
    B[2] = -0.0420650803577195;
    B[3] = 1 - 2 * (B[0] + B[1] + B[2]);
    B[4] = B[2];
    B[5] = B[1];
    B[6] = B[0];
    C[0] = 0.209515106613362;
    C[1] = -0.143851773179818;
    C[2] = 0.5 - (C[0] + C[1]);
    C[3] = C[2];
    C[4] = C[1];
    C[5] = C[0];
    C[6] = 0;

    // Calculate Forces

    for (int i = 0; i < s; i++) {
        temp = current;
        while (temp != 0) {
            temp->setx(temp->getx() + dt * B[i] * temp->getvx());
            temp->sety(temp->gety() + dt * B[i] * temp->getvy());
            temp->setx(temp->getx() - size * floor(temp->getx() / size));
            temp->sety(temp->gety() - size * floor(temp->gety() / size));
            temp = temp->getnext();
        }
        ForceLoop();
        temp = current;
        while (temp != 0) {
            temp->setvx(temp->getvx() + dt * C[i] * temp->getfx());
            temp->setvy(temp->getvy() + dt * C[i] * temp->getfy());
            temp = temp->getnext();
        }
    }
    Kinetic();
}

double System::T() {
    double g = 2 * N;
    double kb = 1.38066 * pow(10, -23);
    return (2 * Ek / g / kb);
}

void System::KillAll() {
    Atom *temp;
    while (current != 0) {
        temp = current->getnext();
        delete current;
        current = temp;
    }
}

System::~System() {
    Atom *temp;
    while (current != 0) {
        temp = current->getnext();
        delete current;
        current = temp;
    }
}

void System::Kinetic() {

    Atom *temp;

    Ek = 0;
    temp = current;
    double tx, ty;
    while (temp != 0) {
        tx = temp->getvx(); // tx now in Ã/s
        ty = temp->getvy();
        /*tx=tx/10000000000;
          ty=ty/10000000000; // scale to m/s */
        Ek = Ek + 0.5 * (tx * tx + ty * ty);
        temp = temp->getnext();
    }
}

double System::getEk() {
    return Ek;
}

double System::getEp() {
    return Ep;
}

void System::ToNewFile(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::trunc);
    temp = current;
    ofile << "Initial Kinetic energy: " << getEk() << "\t" << "(temperature " << T() << " K)" << std::endl << std::endl;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(7);
    while (temp != 0) {
        ofile << temp->getx() << "\t" << temp->gety() << "\t\t" << temp->getvx() << "\t" << temp->getvy() << "\t\t" << temp->getfx() << "\t" << temp->getfy() << "\t" << temp->gettag() << std::endl;
        temp = temp->getnext();
    }
    ofile.close();

}

void System::RemoveLinearMomentum() {

    double lx, ly;
    Atom *temp;

    lx = 0;
    ly = 0;
    temp = current;

    while (temp != 0) {
        lx = lx + temp->getvx();
        ly = ly + temp->getvy();
        temp = temp->getnext();
    }

    lx = lx / N;
    ly = ly / N;

    temp = current;

    while (temp != 0) {
        temp->setvx(temp->getvx() - lx);
        temp->setvy(temp->getvy() - ly);
        temp = temp->getnext();
    }
}

void System::ScaleKineticEnergy(double dEk) {

    Kinetic();
    double s = sqrt(dEk / Ek);
    Atom *temp = current;
    while (temp != 0) {
        temp->setvx(temp->getvx() * s);
        temp->setvy(temp->getvy() * s);
        temp = temp->getnext();
    }
}

void System::NewVitalsToFile(char *fname) {
    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::trunc);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(5);

    ofile << Ek << "\t" << Ep << "\t" << Ep + Ek << std::endl;

    ofile.close();
}

void System::VitalsToFile(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::app);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(5);

    ofile << Ek << "\t" << Ep << "\t" << Ep + Ek << std::endl;

    ofile.close();
}

int System::getN() {
    return N;
}

void System::IntegrateVerlet() {

    Atom *temp;
    temp = current;
    while (temp != 0) {
        temp->setx(temp->getx() + temp->getvx() * dt + temp->getfx()*(0.5 * dt * dt));
        temp->setvx(temp->getvx() + temp->getfx() * dt * 0.5);
        temp->setx(temp->getx() - size * floor(temp->getx() / size));
        temp->sety(temp->gety() + temp->getvy() * dt + temp->getfy()*(0.5 * dt * dt));
        temp->setvy(temp->getvy() + temp->getfy() * dt * 0.5);
        temp->sety(temp->gety() - size * floor(temp->gety() / size));
        temp = temp->getnext();
    }
    ForceLoop();
    temp = current;
    while (temp != 0) {
        temp->setvx(temp->getvx() + temp->getfx() * dt * 0.5);
        temp->setvy(temp->getvy() + temp->getfy() * dt * 0.5);
        temp = temp->getnext();
    }
    Kinetic();
}

void System::AppendToFile(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::app);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << std::setprecision(7) << std::endl;
    while (temp != 0) {
        ofile << temp->getx() << "\t" << temp->gety() << "\t\t" << temp->getvx() << "\t" << temp->getvy() << "\t\t" << temp->getfx() << "\t" << temp->getfy() << "\t" << temp->gettag() << std::endl;
        temp = temp->getnext();
    }
    ofile.close();

}

void System::NewTinker(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::trunc);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << "    " << N << std::endl;
    ofile << std::setprecision(5);

    int k = 1;
    while (temp != 0) {
        ofile << "    " << k << "   Xe   " << temp->getx() << "   " << temp->gety() << "   0   5" << std::endl;
        temp = temp->getnext();
        k++;
    }
    ofile.close();
}

void System::AppendToTinker(char *fname) {

    std::ofstream ofile;
    Atom *temp;

    ofile.open(fname, std::ios::app);
    temp = current;
    ofile.setf(std::ios::fixed);
    ofile << "    " << N << std::endl;
    ofile << std::setprecision(5);

    int k = 1;
    while (temp != 0) {
        ofile << "    " << k << "   Xe   " << temp->getx() << "   " << temp->gety() << "   0   5" << std::endl;
        temp = temp->getnext();
        k++;
    }
    ofile.close();

}

Atom* System::AddAtom(double x, double y, double vx, double vy, double fx, double fy, int tag) {

    Atom *temp;

    if (first == 0) { // If the list is empty...
        current = new Atom(x, y, vx, vy, fx, fy, tag, 0); // Initialize first element
        first = current; // First element will forever be the one just created!
    } else {
        temp = current; // Save pointer to previous element
        current = new Atom(x, y, vx, vy, fx, fy, tag, temp); // Initialize new element
    }
    return current;
}

void System::Initialize(int shape, double kinetic, double cutoff, double timestep, int number, double side, int seed) {

    dt = timestep; // Set up system parameters
    N = number;
    size = side;
    rc = cutoff;
    Ek = kinetic;

    if (seed == 0)
        srand((unsigned) time(NULL)); // Seed random number generator with clock
    else
        srand(seed); // Seed random number generator from input number

    double vx, vy; // Some temporary variables
    double r1, r2;

    if (shape == 1) { // Seed in a box shape, requires N to be a perfect square

        for (int i = 1; i <= sqrt(N); i++) {
            for (int j = 1; j <= sqrt(N); j++) {

                // Compute velocities, sort of bounded to be <= Ek.
                // A MB distribution function would be better

                r1 = (double) rand() / RAND_MAX*Ek; // Generate velocity component 1 <= Ek
                r2 = (double) rand() / RAND_MAX * (Ek - r1); // Generate second velocity component

                // Add the atom to the list
                AddAtom((j - 1) * size / sqrt(N), (i - 1) * size / sqrt(N), r1, r2, 0, 0, (int) ceil((i - 1) * sqrt(N) + j));
            }
        }
    } else if (shape == 2) {
        // Seed in a random manner, interdistance checking should be performed to
        // avoid too much overlap
    } else return;
}

void System::InitializeFromFile(char *fname) {
    return;
}

void System::Distance(int k, int l) {

    Atom *temp;
    temp = current;
    double dx1, dy1, dx2, dy2;

    for (int i = N; i >= ((N - k) + 1); i--) {
        dx1 = temp->getx();
        dy1 = temp->gety();
        temp = temp->getnext();
    }

    temp = current;
    for (int i = N; i >= ((N - l) + 1); i--) {
        dx2 = temp->getx();
        dy2 = temp->gety();
        temp = temp->getnext();
    }
    dx = dx1 - dx2;
    dy = dy1 - dy2;
}

void System::ForceLoop() { // Perform an iterative loop over the list elements

    double rc2 = rc*rc; // Cutoff squared to avoid using sqrts in the loop
    Ep = 0;

    Atom *temp = current;

    while (temp != 0) { // Loop over list to zero forces on atoms
        temp->setfx(0);
        temp->setfy(0);
        temp = temp->getnext();
    }

    temp = current;
    Atom *temp2;

    double d2;
    double fdivr;

    int k = 1;
    int l = 2;

    while (temp != 0) {

        temp2 = temp->getnext();

        while (temp2 != 0) {
            Distance(k, l);
            dx = dx - size * floor(0.5 + dx / size);
            dy = dy - size * floor(0.5 + dy / size);
            d2 = dx * dx + dy*dy;
            if (d2 < rc2) {
                Ep = Ep + 4 * (1 / (d2 * d2 * d2 * d2 * d2 * d2) - 1 / (d2 * d2 * d2)) + 1;
                fdivr = 4 * (12 / (d2 * d2 * d2 * d2 * d2 * d2 * d2) - 6 / (d2 * d2 * d2 * d2));
                temp->setfx(temp->getfx() + dx * fdivr);
                temp2->setfx(temp2->getfx() - dx * fdivr);
                temp->setfy(temp->getfy() + dy * fdivr);
                temp2->setfy(temp2->getfy() - dy * fdivr);
            }
            l++;
            temp2 = temp2->getnext();
        }
        k++;
        l = k + 1;
        temp = temp->getnext();

    }
    return;
}
