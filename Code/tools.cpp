#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include <NTL/RR.h>
#include <random>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "globals.h"

#ifndef TOOLS
#define TOOLS

namespace tools
{
    string logPoly(ZZX poly)
    {
        string output = "(";
        int deg = NTL::deg(poly);
        ZZ c;
        for (int i = 0; i <= deg; i++)
        {
            GetCoeff(c, poly, (long)deg - i);
            output += to_string(conv<int>(c));
            if (i != deg)
            {
                output += " x^" + to_string(deg - i) + " + ";
            }
        }
        output += ")";
        return output;
    }

    string logPoly(ZZ_pX poly)
    {
        return logPoly(conv<ZZX>(poly));
    }

    ZZ SymMod(const ZZ a, const ZZ m)
    {
        ZZ r = a % m;

        if (r <= (m - 1) / 2)
        {
            return r;
        }
        else
        {
            return r - m;
        }
    }

    void resizePolynomialsCoeff(ZZX *polynomial, ZZ mod)
    {
        for (long i = 0; i <= deg(*polynomial); i++)
        {
            (*polynomial)[i] = SymMod((*polynomial)[i], mod);
        }
    }
}

#endif