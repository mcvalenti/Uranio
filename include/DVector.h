#ifndef FVECTOR_H
#define FVECTOR_H
#include <iostream>

using namespace std;

class DVector
{
    private:

        unsigned asize;
    public:
        double* vec;
        DVector();
        DVector(unsigned asize);
        DVector(double* vec, unsigned asize);
        DVector(const DVector &other);
        double& operator[](unsigned index);
        const double& operator[](unsigned index)const;


        const DVector operator+(const DVector &other)const;
        DVector& operator+=(const DVector &other);
        DVector& operator=(const DVector &other);
        const DVector operator*(double val)const;
        virtual ~DVector();
        double getValue(unsigned index);

        friend ostream& operator<<(ostream& os, const DVector& vec);
};

#endif // FVECTOR_H
