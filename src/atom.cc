#include "header/atom.h"

Atom::Atom(double x1, double y1, double vx1, double vy1, double fx1, double fy1, int tag1, Atom *next1) {
    next = next1;
    x = x1;
    y = y1;
    vx = vx1;
    vy = vy1;
    fx = fx1;
    fy = fy1;
    tag = tag1;
}

void Atom::setx(double x1) {
    x = x1;
}

void Atom::sety(double y1) {
    y = y1;
}

void Atom::setxy(double x1, double y1) {
    x = x1;
    y = y1;
}

void Atom::setvxvy(double vx1, double vy1) {
    vx = vx1;
    vy = vy1;
}

void Atom::setvx(double vx1) {
    vx = vx1;
}

void Atom::setvy(double vy1) {
    vy = vy1;
}

void Atom::setfx(double fx1) {
    fx = fx1;
}

void Atom::setfy(double fy1) {
    fy = fy1;
}

void Atom::setfxfy(double fx1, double fy1) {
    fx = fx1;
    fy = fy1;
}

void Atom::settag(int tag1) {
    tag = tag1;
}

void Atom::setnext(Atom *next1) {
    next = next1;
}

double Atom::getx() {
    return x;
}

double Atom::gety() {
    return y;
}

double Atom::getvx() {
    return vx;
}

double Atom::getvy() {
    return vy;
}

double Atom::getfx() {
    return fx;
}

double Atom::getfy() {
    return fy;
}

int Atom::gettag() {
    return tag;
}

Atom* Atom::getnext() {
    return next;
}



