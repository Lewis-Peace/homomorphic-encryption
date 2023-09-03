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

using namespace NTL;
using namespace std;

namespace randomness
{
    gsl_rng *r;

    void setSeed()
    {
        r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, time(0));
    }

    int getNumber()
    {
        int sample;
        do
        {
            sample = round(gsl_ran_gaussian(r, 3.19));
        } while (sample < 0 || sample > B);
        return sample;
    }

    ZZ getNumber(ZZ from, ZZ to)
    {
        ZZ_pPush push(abs(from) + abs(to) + 1);
        ZZ_p s;
        NTL::random(s);
        return conv<ZZ>(s) + from;
    }

    /**
     * @brief Returns a poylnomial with all the coefficients equal to
     * a random number obtained from the discrete gaussian distribution
     *
     * @param length
     * @return ZZX
     */
    ZZX NormalPolynomial(int length, long mean = 0, long stddev = sigma)
    {
        ZZX r;
        r.SetLength(length);

        // std::default_random_engine generator;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(mean, stddev);

        int coeff = getNumber();
        for (long i = 0; i < length; i++)
        {
            SetCoeff(r, i, NTL::ZZ((int)distribution(gen)));
        }
        return r;
    }

    /**
     * @brief Returns a polynomial with random coeficients bound to the limit
     * of the arguments from and to
     *
     * @param from
     * @param to
     * @param length
     * @return ZZX
     */
    ZZX RandomPolynomial(ZZ from, ZZ to, int length)
    {
        ZZX r;
        r.SetLength(length);
        for (long i = 0; i < length; i++)
        {
            r[i] = getNumber(from, to);
        }
        return r;
    }

    ZZX RandomPolynomial(int from, int to, int length)
    {
        return randomness::RandomPolynomial(conv<ZZ>(from), conv<ZZ>(to), length);
    }

    void clear()
    {
        gsl_rng_free(r);
    }
}