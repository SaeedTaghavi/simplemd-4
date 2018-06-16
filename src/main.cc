#include "header/system.h"
#include "header/atom.h"
#include <cstdio>
#include <iostream>
#include <math.h>

using namespace std;

int main(void) {

    System Box; // Create the box

    double E = 0;
    double simultime1 = 0.1;
    double simultime2 = 10;
    double verletstep = 0.001;
    double eulerstep = 0.001;
    double rungestep = 0.001;
    int seed = 123456; // Random seed number; 0 to generate from clock

    cout << "Starting simulation..." << endl << endl;
    cout << "Verlet Timestep is " << verletstep << " seconds." << endl;
    cout << "SRK Timestep is " << rungestep << " seconds." << endl;
    cout << "Euler Timestep is " << eulerstep << " seconds." << endl;
    cout << "Random seed is " << seed << " (0 generates seed from clock)" << endl << endl;
    cout << "Now running Velocity Verlet over " << round(simultime1 / verletstep) << " timesteps..." << endl;

    // Quick comparison over t=0.1

    Box.Initialize(1, 10, 1.122462, verletstep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultVerlet100.dat");
    Box.NewVitalsToFile("vitalsVerlet100.dat");
    Box.NewTinker("tinkerVerlet100.dat");
    Box.NewSingleToFile("singleVerlet100.dat");

    for (int k = 1; k <= round(simultime1 / verletstep); k++) {

        Box.IntegrateVerlet();

        Box.AppendToFile("resultVerlet100.dat");
        Box.VitalsToFile("vitalsVerlet100.dat");
        if (k % 1 == 0) {
            Box.AppendToTinker("tinkerVerlet100.dat");
            Box.SingleToFile("singleVerlet100.dat");
        }
        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl;
    cout << "Now running Euler Forward over " << round(simultime1 / eulerstep) << " timesteps..." << endl;

    Box.Initialize(1, 10, 1.122462, eulerstep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultEuler100.dat");
    Box.NewVitalsToFile("vitalsEuler100.dat");
    Box.NewTinker("tinkerEuler100.dat");
    Box.NewSingleToFile("singleEuler100.dat");

    for (int k = 1; k <= round(simultime1 / eulerstep); k++) {

        Box.IntegrateEuler();

        Box.AppendToFile("resultEuler100.dat");
        Box.VitalsToFile("vitalsEuler100.dat");
        if (k % 1 == 0) {
            Box.AppendToTinker("tinkerEuler100.dat");
            Box.SingleToFile("singleEuler100.dat");
        }

        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl;
    cout << "Now running Partition Runge-Kutta over " << round(simultime1 / rungestep) << " timesteps..." << endl;

    Box.Initialize(1, 10, 1.122462, rungestep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultRungeKutta100.dat");
    Box.NewVitalsToFile("vitalsRungeKutta100.dat");
    Box.NewTinker("tinkerRungeKutta100.dat");
    Box.NewSingleToFile("singleRungeKutta100.dat");

    for (int k = 1; k <= round(simultime1 / rungestep); k++) {

        Box.IntegratePartitionRungeKutta();
        Box.AppendToFile("resultRungeKutta100.dat");
        Box.VitalsToFile("vitalsRungeKutta100.dat");
        if (k % 1 == 0) {
            Box.AppendToTinker("tinkerRungeKutta100.dat");
            Box.SingleToFile("singleRungeKutta100.dat");
        }
        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl;
    cout << "Now running Partition Runge-Kutta over " << round(simultime2 / rungestep) << " timesteps..." << endl;

    Box.Initialize(1, 10, 1.122462, rungestep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultRungeKutta10000.dat");
    Box.NewVitalsToFile("vitalsRungeKutta10000.dat");
    Box.NewTinker("tinkerRungeKutta10000.dat");
    Box.NewSingleToFile("singleRungeKutta10000.dat");

    for (int k = 1; k <= round(simultime2 / rungestep); k++) {
        E = Box.getEp() + Box.getEk();
        Box.IntegratePartitionRungeKutta();
        Box.AppendToFile("resultRungeKutta10000.dat");
        Box.VitalsToFile("vitalsRungeKutta10000.dat");
        if (k % 10 == 0)
            Box.AppendToTinker("tinkerRungeKutta10000.dat");
        if (k % 10 == 0)
            Box.SingleToFile("singleRungeKutta10000.dat");

        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
        E = E - Box.getEk() - Box.getEp();
        E = E*E;
        if (k > 200 && E > 10000) {
            cout << "System has diverged (dE > 100) at time " << rungestep * k << "(step " << k << "), terminating run!" << endl;
            break;
        }
    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl;
    cout << "Now running Euler Forward over " << round(simultime2 / eulerstep) << " timesteps..." << endl;

    Box.Initialize(1, 10, 1.122462, eulerstep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultEuler10000.dat");
    Box.NewVitalsToFile("vitalsEuler10000.dat");
    Box.NewTinker("tinkerEuler10000.dat");
    Box.NewSingleToFile("singleEuler10000.dat");

    for (int k = 1; k <= round(simultime2 / eulerstep); k++) {
        E = Box.getEk() + Box.getEp();
        Box.IntegrateEuler();
        Box.AppendToFile("resultEuler10000.dat");
        Box.VitalsToFile("vitalsEuler10000.dat");
        if (k % 10 == 0)
            Box.AppendToTinker("tinkerEuler10000.dat");

        if (k % 10 == 0)
            Box.SingleToFile("singleEuler10000.dat");

        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
        E = E - Box.getEk() - Box.getEp();
        E = E*E;
        if (k > 200 && E > 10000) {
            cout << "System has diverged (dE > 100) at time " << eulerstep * k << "(step " << k << "), terminating run!" << endl;
            break;
        }

    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl;
    cout << "Now running Velocity Verlet over " << round(simultime2 / verletstep) << " timesteps..." << endl;

    Box.Initialize(1, 10, 1.122462, verletstep, 25, 5, seed); // Initialize in a square pattern
    Box.ForceLoop();
    Box.Kinetic();

    Box.ToNewFile("resultVerlet10000.dat");
    Box.NewVitalsToFile("vitalsVerlet10000.dat");
    Box.NewTinker("tinkerVerlet10000.dat");
    Box.NewSingleToFile("singleVerlet10000.dat");

    for (int k = 1; k <= round(simultime2 / verletstep); k++) {
        E = Box.getEp() + Box.getEk();
        Box.IntegrateVerlet();
        Box.AppendToFile("resultVerlet10000.dat");
        Box.VitalsToFile("vitalsVerlet10000.dat");
        if (k % 10 == 0)
            Box.AppendToTinker("tinkerVerlet10000.dat");

        if (k % 10 == 0)
            Box.SingleToFile("singleVerlet10000.dat");


        if (k % 10 == 0 && k <= 200) {
            Box.ScaleKineticEnergy(10 * Box.getN());
            Box.RemoveLinearMomentum();
        }
        E = E - Box.getEk() - Box.getEp();
        E = E*E;
        if (k > 200 && E > 10000) {
            cout << "System has diverged (dE > 100) at time " << verletstep * k << "(step " << k << "), terminating run!" << endl;
            break;
        }
    }

    Box.KillAll(); // Remove old data tree

    cout << "Done." << endl << endl << "Simulation finished!" << endl;

    return 1;
}
