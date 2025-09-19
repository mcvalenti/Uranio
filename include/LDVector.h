/*
 * LLDVector.h
 *
 *  Created on: 19 may. 2023
 *      Author: ceci
 */

#ifndef LLDVector_H_
#define LLDVector_H_

#include <iostream>

using namespace std;

class LDVector
{
    private:

        unsigned asize;
    public:
        long double* vec;

        // Constructors
        LDVector();
        LDVector(unsigned asize);
        LDVector(long double* vec, unsigned asize);
        LDVector(const LDVector &other);

        // Operators
        long double& operator[](unsigned index);
        const long double& operator[](unsigned index)const;
        const LDVector operator+(const LDVector &other)const;
        const LDVector operator-(const LDVector &other)const;
        LDVector& operator+=(const LDVector &other);
        LDVector& operator=(const LDVector &other);
        const LDVector operator*(long double val)const;
        friend ostream& operator<<(ostream& os, const LDVector& vec);

        // Methods
        //bool areEqual(const LDVector &other, double tolerance);
        //double get_max_absolute(const LDVector &other);
        long double getValue(unsigned index);

        virtual ~LDVector();

};





#endif /* LLDVector_H_ */
