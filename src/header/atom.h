#ifndef ATOM
#define ATOM

class Atom {
private:
    int tag; // ID number in list
    double x, y; // Position
    double vx, vy; // Velocity
    double fx, fy; // Force
    Atom *next;
protected:
public:

    Atom() : x(0), y(0), vx(0), vy(0), fx(0), fy(0), tag(0), next(0) {
    }
    Atom(double x1, double y1, double vx1, double vy1, double fx1, double fy1, int tag1, Atom *next1);
    void setxy(double x1, double y1); // Various methods to assign stuff
    void setvxvy(double vx1, double vy1);
    void setfxfy(double fx1, double fy1);
    void setx(double x1);
    void sety(double x2);
    void setvx(double vx1);
    void setvy(double vy1);
    void setfx(double fx1);
    void setfy(double fy1);
    void settag(int tag1);
    void setnext(Atom *next);
    double getx(); // Various methods to read the class
    double gety();
    double getvx();
    double getvy();
    double getfx();
    double getfy();
    int gettag();
    Atom* getnext();
};

#endif
