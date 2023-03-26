#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
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
namespace randomness {
    gsl_rng* r;

    void setSeed() {
        r = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(r, time(0));
    }

    int getNumber() {
        int sample;
        do {
            sample = round(gsl_ran_gaussian(r, 3.19));
        } while (sample < 0 || sample > B);
        return sample;
    }
    
    ZZ_pX randomPolynomial(int length) {
        ZZ_pX r;
        int coeff = getNumber();
        r.SetLength(length);
        for (long i = 0; i < length; i++) {
            SetCoeff(r, i, coeff);
        }
        return r;
    }

    void clear() {
        gsl_rng_free(r);
    }
}


class FanVercScheme {
    
    public:
        int n = 4, b = 2, w = 2, l;

        ZZ_pX sk; // Secret key of the class
        ZZ_pX irreduciblePolynomial;
        ZZ_pX delta;
        ZZ q; // Modulus in the ciphertext space
        ZZ t; // Modulus in the plaintext space
        tuple<ZZ_pX, ZZ_pX> pk; // Public key tuple of the class
        Vec<tuple<ZZ_pX, ZZ_pX>> evk;

        /**
         * @brief Construct a new Fan-Vercauteren Scheme object initializing the modulo, the secret key, the public key
         * and the evaluation keys.
         */
        FanVercScheme(int n, int b, int q, int t, int w = 2) {
            randomness::setSeed();
            printf("Setting n as %i and b as %i\n", n, b);
            this->n = n; this->b = b; this->w = w; this->q = q; this->t = t; this->l = floor(log(this->q) / log(this->w));
            ZZ_p::init(this->q);
            
            // Calculating ∆b
            double delta_multiplicative = -(conv<double>(this->q) / (pow(this->b, this->n) + 1));
            this->delta.SetLength(this->n);
            for (int i = 0; i < this->n; i++) {
                SetCoeff(this->delta, i, conv<long>(delta_multiplicative * pow(this->b, this->n - (i + 1))));
            }
            
            irreduciblePolynomial.SetLength(n);
            SetCoeff(irreduciblePolynomial, 0);
            SetCoeff(irreduciblePolynomial, n);
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
        ZZ_pX fractionalEncoder (int numerator, int denominator, int sign) {
            try {
                ZZ_p x;
                power(x, (ZZ_p) denominator, (long) -1);
                x = x * numerator;
                return FanVercScheme::polyConversion(x, sign);
            } catch(const std::exception& e) {
                printf("Error: No inverse for the number %i in modulo ", denominator);
                cout << this->q << "\n";
                exit(1);
            }
        }

        // TODO Check if encrypts correctly
        tuple<ZZ_pX, ZZ_pX> encrypt(int numerator, int denominator, int sign) {
            ZZ_pX message = fractionalEncoder(numerator, denominator, sign);
            ZZ_pX u = random_ZZ_pX(this->n);
            ZZ_pX e0, e1;
            ZZ_pX c0, c1, pk0, pk1;
            tie(pk0, pk1) = this->pk;
            e0 = randomness::randomPolynomial(this->n);
            e1 = randomness::randomPolynomial(this->n);
            cout << "e0 = " << e0 << "\n";
            cout << "e1 = " << e1 << "\n";
            cout << "delta = " << this->delta << "\n";
            
            c0 = reduce(message * this->delta) + reduce(pk0 * u) + e0;
            c1 = reduce(pk1 * u) + e1;
            cout << "c0 = " << c0 << "\n";
            cout << "c1 = " << c1 << "\n";
            return make_tuple(c0, c1);
        }

        // TODO Fix decryption function
        ZZ_pX decrypt(tuple<ZZ_pX, ZZ_pX> ct) {
            ZZ_pX c0, c1, M, temp;
            tie(c0, c1) = ct;
            temp.SetLength(2);
            printf("b = %i q = %i\n", this->b, conv<int>(this->q));
            SetCoeff(temp, 0, b);
            SetCoeff(temp, 1);
            M = reduce((c0 + reduce(c1 * this->sk)) * temp);
            
            cout << "M before " << M << "\n";
            double messageCoefficent;
            for(int i = 0; i <= deg(M); i++) {
                messageCoefficent = conv<double>(conv<ZZ>(coeff(M, i))) / conv<double>(this->q);
                printf("%f ", messageCoefficent);
                SetCoeff(M, i, round(messageCoefficent));
            }
            printf("\n");

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
    
    private:
    
        /**
         * @brief Converts an integer value to a polynomial based on his binary representation
         * 
         * @param value The integer to convert
         * @return The polynomial of type ZZX
         */
        static ZZ_pX polyConversion(ZZ_p value, int sign) {
            ZZ temp = conv<ZZ>(value);
            int polyLength = 0;
            while (temp != 0) {
                temp = temp >> 1;
                polyLength++;
            }
            temp = conv<ZZ>(value);
            ZZ_pX f;
            f.SetLength(polyLength);
            for (long i = 0; i < polyLength; i++) {
                SetCoeff(f, i, (long)(conv<int>((temp >> i) & 1)));
            }
            cout << "encoding = " << f << "\n";
            return f * sign;
        }
        
        /**
         * @brief Generates a random polynomial with random coeficients (-1, 0, 1)
         * 
         */
        void secretKeyGen() {
            this->sk = random_ZZ_pX(this->n);
            cout << "sk = " << sk << "\n";
            printf("Private key generated....\n");
        }

        /**
         * @brief Generates and stores the public key of the cipher
         * 
         * @param sk The secret key required to generate the pubblic key
         */
        void publicKeyGen() {
            ZZ_pX pk1 = random_ZZ_pX(this->n), pk0;
            ZZ_pX e = randomness::randomPolynomial(this->n);

            cout << "e = " << e << "\n";
            cout << "b = " << pk1 << "\n";
            pk0 = -((pk1 * this->sk) + e);
            reduce(&pk0);
            cout << "a = " << pk0 << "\n";

            this->pk = make_tuple(pk0, pk1);
            printf("Public keys generated....\n");
        }
        
        
        void evaluationKeyGen() {
            Vec<tuple<ZZ_pX, ZZ_pX>> pairs;
            pairs.SetLength(l);
            int w = 2;
            for (int i = 0; i < this->l; i++) {
                ZZ_pX x = random_ZZ_pX(this->n), y;
                ZZ_pX s, e;
                s = randomness::randomPolynomial(this->n);
                e = randomness::randomPolynomial(this->n);
                y = -((x * this->sk) + e) + pow(w, i) * (s * s);
                reduce(&y);
                
                // cout << "i = " << i << "\n";
                // cout << "y = " << y << "\n";
                // cout << "x = " << x << "\n";
                pairs[i] = make_tuple(y, x);
            }
            this->evk = pairs;
            printf("Evaluation keys generated....\n");
        }


        void reduce(ZZ_pX* polynomial) {
            while (deg(*polynomial) > n) {
                *polynomial /= this->irreduciblePolynomial;
            }
        }

        ZZ_pX reduce(ZZ_pX polynomial) {
            while (deg(polynomial) > n) {
                polynomial /= this->irreduciblePolynomial;
            }
            return polynomial;
        }

        // ZZ mul(ZZ ct0, ZZ ct1) {
        // }

};


 void testing() {
    ZZ_pX x, y, z;
    cin >> x;
    cin >> y;
    z = x * y;
    cout << z << "\n";
 }

 int main() {
    int numerator, denominator;
    printf("Insert numerator and denominator: ");
    int dataReceived = scanf("%i %i", &numerator, &denominator);
    int q = 73, n = 4, b = 2;
    FanVercScheme f(n, b, q, (q / 8) - 1);
    f.decrypt(f.encrypt(numerator, denominator, -1));
    randomness::clear();
    return 0;
}
