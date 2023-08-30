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

#define sigma 3.19
#define B 6 * sigma

using namespace NTL;
using namespace std;

/**
 * @brief All the preudo random number generation
 * will be handled in this part of the code
 *
 */
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
    ZZX randomPolynomial(int length)
    {
        ZZX r;
        int coeff = getNumber();
        r.SetLength(length);
        for (long i = 0; i < length; i++)
        {
            SetCoeff(r, i, coeff);
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
    ZZX randomPolynomial(ZZ from, ZZ to, int length)
    {
        ZZX r;
        r.SetLength(length);
        for (long i = 0; i < length; i++)
        {
            r[i] = getNumber(from, to);
        }
        return r;
    }

    ZZX randomPolynomial(int from, int to, int length)
    {
        return randomness::randomPolynomial(conv<ZZ>(from), conv<ZZ>(to), length);
    }

    void clear()
    {
        gsl_rng_free(r);
    }
}

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

using namespace tools;

class FanVercScheme
{

public:
    int n, w, l;

    ZZX sk; // Secret key of the class
    ZZX irreduciblePolynomial;
    ZZX delta;
    ZZ from;
    ZZ to;
    ZZ b;
    ZZ q;               // Modulus in the ciphertext space
    ZZ t;               // Modulus in the plaintext space
    tuple<ZZX, ZZX> pk; // Public key tuple of the class
    Vec<tuple<ZZX, ZZX>> evk;

    /**
     * @brief Construct a new Fan-Vercauteren Scheme object initializing the modulo, the secret key, the public key
     * and the evaluation keys.
     */
    FanVercScheme(int n, int b, ZZ q, int w = 2)
    {
        randomness::setSeed();
        this->from = -(q - 1) / 2;
        this->to = (q - 1) / 2;
        this->n = n;
        this->b = b;
        this->w = w;
        this->q = q;
        this->t = t;
        this->l = floor(log(this->q) / log(this->w));
        ZZ_p::init(this->q);

        printf("Setting n as %i and b as %i and q as %i\n", n, b, conv<int>(this->q));

        // Calculating ∆b
        RR delta_multiplicative = -(conv<RR>(this->q) / conv<RR>(power(conv<RR>(this->b), this->n) + 1));
        this->delta.SetLength(this->n);
        for (int i = 0; i < this->n; i++)
        {
            this->delta[i] = conv<ZZ>(round(delta_multiplicative * power(conv<RR>(this->b), this->n - (i + 1))));
        }

        this->irreduciblePolynomial.SetLength(n + 1);
        this->irreduciblePolynomial[0] = 1;
        this->irreduciblePolynomial[n] = 1;

        secretKeyGen();
        publicKeyGen();
        evaluationKeyGen();
    }

    /**
     * @brief Converts a fraction to a polynomial form based on the fractional encoder method
     * assuming n and d are numerator and denominator computes n*d^-1 mod q and than converts
     * the result binary representation in a polynomial
     *
     * @param numerator Numerator of the fraction
     * @param denominator Denominator of the fraction
     * @return ZZX The polynomial
     */
    ZZX fractionalEncoder(int numerator, int denominator)
    {
        ZZ tempModulo = power(this->b, this->n) + 1;
        printf("fraction %i/%i\n", numerator, denominator);
        ZZ_pPush push(tempModulo);
        try
        {
            ZZ_p x;
            power(x, (ZZ_p)denominator, (long)-1);
            cout << "inverse of " << denominator << " = " << x << " in modulo " << tempModulo;
            x = x * numerator;
            cout << " - Frac to int = " << x << "\n";
            return FanVercScheme::polyConversion(conv<ZZ>(x));
        }
        catch (const std::exception &e)
        {
            printf("Error: No inverse for the number %i in modulo ", denominator);
            cout << this->q << "\n";
            exit(1);
        }
    }

    // TODO Check if encrypts correctly
    tuple<ZZX, ZZX> encrypt(int numerator, int denominator)
    {
        ZZX message = fractionalEncoder(numerator, denominator);
        cout << "message = " << logPoly(message) << "\n";
        cout << "message = " << message << "\n";
        ZZX u = randomness::randomPolynomial(-1, 1, this->n);
        ZZX e0, e1, c0, c1, pk0, pk1;
        tie(pk0, pk1) = this->pk;
        e0 = randomness::randomPolynomial(this->n);
        e1 = randomness::randomPolynomial(this->n);
        cout << "u = " << logPoly(u) << "\n";
        cout << "e0 = " << logPoly(e0) << "\n";
        cout << "e1 = " << logPoly(e1) << "\n";
        cout << "delta = " << logPoly(this->delta) << "\n";

        c0 = MulMod(message, this->delta, this->irreduciblePolynomial) + MulMod(pk0, u, this->irreduciblePolynomial) + e0;
        c1 = MulMod(pk1, u, this->irreduciblePolynomial) + e1;
        resizePolynomialsCoeff(&c0, this->q);
        resizePolynomialsCoeff(&c1, this->q);
        cout << "c0 = " << logPoly(c0) << "\n";
        cout << "c1 = " << logPoly(c1) << "\n";
        return make_tuple(c0, c1);
    }

    // TODO Fix decryption function
    ZZX decrypt(tuple<ZZX, ZZX> ct)
    {
        ZZX c0, c1, M;
        ZZX temp;
        temp.SetLength(2);
        temp[0] = -this->b;
        temp[1] = 1;
        tie(c0, c1) = ct;
        M = MulMod(c1, this->sk, this->irreduciblePolynomial);
        M += c0;
        resizePolynomialsCoeff(&M, this->q);
        M = MulMod(M, temp, this->irreduciblePolynomial);
        cout << "M before " << logPoly(M) << "\n";
        RR q, coeff;
        conv(q, this->q);
        for (int i = 0; i < this->n; i++)
        {
            conv(coeff, M[i]);
            M[i] = conv<ZZ>(round(coeff / q));
        }
        cout << "M = " << M << "\n";
        resizePolynomialsCoeff(&M, this->b);
        cout << "M = " << M << "\n";
        return M;
    }

    // ZZ add(ZZ ct0, ZZ ct1) {

    // }

    // ZZ multiply(ZZ ct0, ZZ ct1, ZZ evk) {
    //     relinearize(mul(ct0, ct1), evk);
    // }

    // ZZ relinearize(ZZ ct, ZZ evk) {

    // }

    // private:

    /**
     * @brief Converts an integer value to a polynomial based on his binary representation
     *
     * @param value The integer to convert
     * @return The polynomial of type ZZX
     */
    ZZX polyConversion(ZZ value)
    {
        ZZ temp = value;
        ZZX f;
        f.SetLength(this->n);
        for (long i = 0; i < this->n; i++)
        {
            f[i] = (value >> i) & 1;
            printf("%i ", conv<int>((value >> i) & 1));
        }
        cout << "\nencoding = " << logPoly(f) << "\n";
        cout << "encoding = " << f << "\n";
        printf("Conversion completed...\n\n");
        return f;
    }
    /**
     * @brief Generates a random polynomial with random coeficients (-1, 0, 1)
     *
     */
    void secretKeyGen()
    {
        this->sk = randomness::randomPolynomial(-1, 1, this->n);
        cout << "sk = " << logPoly(sk) << "\n";
        printf("Private key generated....\n\n");
    }

    /**
     * @brief Generates and stores the public key of the cipher
     *
     * @param sk The secret key required to generate the pubblic key
     */
    void publicKeyGen()
    {
        ZZX pk1 = randomness::randomPolynomial(this->from, this->to, this->n), pk0;
        ZZX e = randomness::randomPolynomial(this->n);

        cout << "e = " << logPoly(e) << "\n";
        cout << "pk1 = " << logPoly(pk1) << "\n";
        pk0 = -(MulMod(pk1, this->sk, this->irreduciblePolynomial) + e);
        resizePolynomialsCoeff(&pk0, this->q);
        cout << "pk0 = " << logPoly(pk0) << "\n";

        this->pk = make_tuple(pk0, pk1);
        printf("Public keys generated....\n\n");
    }

    void evaluationKeyGen()
    {
        Vec<tuple<ZZX, ZZX>> pairs;
        pairs.SetLength(l);
        int w = 2;
        for (int i = 0; i < this->l; i++)
        {
            ZZX x = randomness::randomPolynomial(this->from, this->to, this->n), y;
            ZZX s, e;
            s = randomness::randomPolynomial(this->n);
            e = randomness::randomPolynomial(this->n);
            y = -(MulMod(x, this->sk, this->irreduciblePolynomial) + e) + pow(w, i) * (s * s);

            // cout << "i = " << i << "\n";
            // cout << "y = " << y << "\n";
            // cout << "x = " << x << "\n";
            pairs[i] = make_tuple(y, x);
        }
        this->evk = pairs;
        printf("Evaluation keys generated....\n\n");
    }

    void resizePolynomialsCoeff(ZZX *polynomial, ZZ mod)
    {
        tools::resizePolynomialsCoeff(polynomial, mod);
    }

    // ZZ mul(ZZ ct0, ZZ ct1) {
    // }
};

int main(int argc, char *argv[])
{
    if (argc >= 3)
    {
        int numerator = atoi(argv[1]), denominator = atoi(argv[2]);
        ZZ q = power2_ZZ(60);
        tuple<ZZX, ZZX> ct;
        int n = 2048, b = 2;
        FanVercScheme f(n, b, q);
        ct = f.encrypt(numerator, denominator);
        f.decrypt(ct);
        randomness::clear();
    }
    else
    {
        printf("Too few arguments\n");
    }

    return 0;
}
