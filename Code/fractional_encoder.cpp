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
#include "tools.cpp"

using namespace tools;

namespace fractionalEncoder
{
    class Rational
    {
    public:
        NTL::ZZ numerator;
        NTL::ZZ denominator;

        Rational(const long num, const long den)
        {
            NTL::ZZ numerator_ = NTL::ZZ(num);
            NTL::ZZ denominator_ = NTL::ZZ(den);
            NTL::ZZ g = NTL::GCD(numerator_, denominator_);
            this->numerator = numerator_ / g;
            this->denominator = denominator_ / g;
        };

        Rational(NTL::ZZ numerator_, NTL::ZZ denominator_)
        {
            this->numerator = std::move(numerator_);
            this->denominator = std::move(denominator_);
            GCD();
        }

        void GCD()
        {
            NTL::ZZ g = NTL::GCD(this->numerator, this->denominator);
            this->numerator = numerator / g;
            this->denominator = denominator / g;
        }
    };

    class HenselCode
    {
    public:
        NTL::ZZ prime;
        int r;
        NTL::ZZ modulus;
        NTL::ZZ code;
        HenselCode(const NTL::ZZ &prime_, const long &r_exponent, const NTL::ZZ &code)
        {
            this->prime = prime_;
            this->r = r_exponent;
            this->modulus = NTL::power(this->prime, this->r);
            this->code = code;
        }

        static HenselCode Encode(const NTL::ZZ &prime_, const long &r_exponent, Rational &m_prime)
        {
            NTL::ZZ den_inv = InvMod(m_prime.denominator, power(prime_, r_exponent));
            NTL::ZZ hencoded =
                MulMod(m_prime.numerator, den_inv, power(prime_, r_exponent));

            HenselCode h = HenselCode(prime_, r_exponent, hencoded);
            return h;
        }

        static Rational Decode(const NTL::ZZ &prime_, const long &r_exponent, const HenselCode &hensel_code)
        {
            int index = 0;
            NTL::ZZ big_n = SqrRoot((prime_ - 1) / 2);
            NTL::Vec<NTL::ZZ> x;
            x.SetLength(2);
            x[index] = power(prime_, r_exponent);
            x[index + 1] = hensel_code.code;

            NTL::Vec<NTL::ZZ> y;
            y.SetLength(2);
            y[index] = NTL::ZZ(0);
            y[index + 1] = NTL::ZZ(1);

            NTL::Vec<NTL::ZZ> z;
            z.SetLength(2);
            z[index] = NTL::ZZ(1);
            z[index + 1] = NTL::ZZ(0);
            NTL::ZZ q;
            index = 1;

            while ((x[index] > big_n) != 0)
            {
                x.SetLength(x.length() + 2);
                y.SetLength(y.length() + 2);
                z.SetLength(z.length() + 2);
                q = x[index - 1] / x[index];
                x[index + 1] = x[index - 1] - q * x[index];
                y[index + 1] = y[index - 1] - q * y[index];
                z[index + 1] = z[index - 1] - q * z[index];
                index++;
            }

            NTL::ZZ c = NTL::power(NTL::ZZ(-1), (index + 1)) * x[index];
            NTL::ZZ d = NTL::power(NTL::ZZ(-1), (index + 1)) * y[index];

            Rational m = Rational(c, d);

            return m;
        }

        static NTL::ZZX PolyEncode(const long &b, const long &n, Rational m)
        {
            NTL::ZZ g = NTL::ZZ(NTL::power_ZZ(b, n) + 1);
            HenselCode hc = Encode(g, 1, m);

            NTL::ZZ h = hc.code;

            NTL::ZZX hx;
            hx.SetLength(n);

            for (long i = 0; i < n; i++)
            {
                hx[i] = tools::SymMod(h / NTL::power_ZZ(b, i), NTL::ZZ(b));
            }

            return NTL::reverse(hx);
        }

        static Rational PolyDecode(const long &b, const long &n, NTL::ZZX hx)
        {
            NTL::ZZ g = NTL::ZZ(NTL::power_ZZ(b, n) + 1);
            NTL::ZZ h = NTL::ZZ(0);

            hx = NTL::reverse(hx);

            for (long i = 0; i < n; i++)
            {
                h += (hx[i] % b) * NTL::power_ZZ(b, i);
            }

            HenselCode hc = HenselCode(g, 1, h);

            return Decode(g, 1, hc);
        }
    };
}
