#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vector.h>
#include <stdio.h>
#include <stdlib.h>
#include <tuple>
#include "DiscreteGaussianNumber.cpp"


// A power of 2 tipically 1024
#define n 4
#define b 2

using namespace NTL;
using namespace std;

class FanVercScheme {
    
    public:
        ZZ_pX sk; // Secret key of the class
        ZZ q = power_ZZ(b, n) + 1; // Modulus in the ciphertext space
        tuple<ZZ_pX, ZZ_pX> pk; // Public key tuple of the class
        Vec<tuple<ZZ_pX, ZZ_pX>> evk;
        DiscreteGaussianNumber dgn;

        /**
         * @brief Construct a new Fan-Vercauteren Scheme object initializing the modulo, the secret key, the public key
         * and the evaluation keys.
         */
        FanVercScheme() {
            ZZ_p::init(this->q);
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
            try { // TODO: Handle the not initialized modulo possibility
                ZZ_p x;
                power(x, (ZZ_p) denominator, (long) -1);
                printf("%i ^ -1 = ", denominator);
                cout << x << "\n";
                x = x * numerator;
                printf("%i/%i mod ", numerator, denominator);
                cout << this->q << " = " << x << "\n";
                return FanVercScheme::polyConversion(x, sign);
            } catch(const std::exception& e) {
                printf("Error: No inverse for the number %i in modulo ", denominator);
                cout << this->q << "\n";
                exit(1);
            }
        }

        // FIXME: change from usage of modulo numbers to modulo polynomials
        tuple<ZZ_pX, ZZ_pX> encrypt(int numerator, int denominator, int sign) {
            ZZ_pX message = fractionalEncoder(numerator, denominator, sign);
            ZZ_pX u = random_ZZ_pX(n);
            double e0, e1;
            ZZ_pX c0, c1, pk0, pk1, delta;
            tie(pk0, pk1) = this->pk;
            e0 = gaussianRandomNumber();
            e1 = gaussianRandomNumber();

            ZZ delta_coefficient = -(q / (pow(b, n) + 1));
            delta.SetLength(n);
            for (long i = 0; i < n; i++) {
                SetCoeff(delta, i, conv<long>(delta_coefficient * pow(b, n - (i + 1))));
            }
            cout << "delta = " << delta << "\n";
            
            c0 = message * delta + pk0 * u + e0;
            c1 = pk1 * u + e1;
            return make_tuple(c0, c1);
        }

        ZZ_pX decrypt(tuple<ZZ_pX, ZZ_pX> ct) {
            ZZ_pX c0, c1, M, temp;
            tie(c0, c1) = ct;
            temp.SetLength(2);
            SetCoeff(temp, 0, conv<long>(-(conv<ZZ>(b) / this->q)));
            SetCoeff(temp, 1);
            M = (c0 + c1 * this->sk) * temp;
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
        
        double gaussianRandomNumber() {
            return this->dgn.getNumber();
        }
    
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
        // TODO: add change of sign for each element
        void secretKeyGen() {
            ZZX s;
            s.SetLength(n);
            srandom(time(0));
            for (int i = 0; i < n; i++) {
                SetCoeff(s, i, random() * ((random() % 3) - 1));
            }
            cout << "sk ZZX = " << s << "\n";
            this->sk = random_ZZ_pX(n);
            cout << "sk = " << sk << "\n";
            printf("Private key generated....\n");
        }

        /**
         * @brief Generates and stores the public key of the cipher
         * 
         * @param sk The secret key required to generate the pubblic key
         */
        void publicKeyGen() {
            ZZ_pX pk1 = random_ZZ_pX(n), pk0;
            double e = gaussianRandomNumber();

            cout << "e = " << e << "\n";
            cout << "b = " << pk1 << "\n";
            pk0 = -((pk1 * this->sk) + e);
            cout << "a = " << pk0 << "\n";

            this->pk = make_tuple(pk0, pk1);
            printf("Public keys generated....\n");
        }
        
        
        void evaluationKeyGen() {
            Vec<tuple<ZZ_pX, ZZ_pX>> pairs;
            double l = floor(log2(conv<double>(this->q)));
            pairs.SetLength(l);
            int w = 2;
            for (int i = 0; i < l; i++) {
                ZZ_pX x = random_ZZ_pX(n), y;
                double s, e;
                s = gaussianRandomNumber();
                e = gaussianRandomNumber();
                y = -((x * this->sk) + e) + pow(w, i) * (s * s);
                
                cout << "i = " << i << "\n";
                cout << "y = " << y << "\n";
                cout << "x = " << x << "\n";
                pairs[i] = make_tuple(y, x);
            }
            this->evk = pairs;
            printf("Evaluation keys generated....\n");
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
    int numerator = 8, denominator = 3;
    ZZX s;
    FanVercScheme f;
    f.decrypt(f.encrypt(numerator, denominator, -1));
}
