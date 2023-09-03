#include "randomness.cpp"
#include "tools.cpp"
#include "fractional_encoder.cpp"
#include "globals.h"

using namespace NTL;
using namespace std;
using namespace randomness;
using namespace tools;
using namespace fractionalEncoder;

class FanVercScheme
{

public:
    int n, l;

    ZZX sk; // Secret key of the class
    ZZX irreduciblePolynomial;
    ZZX delta;
    ZZ from;
    ZZ to;
    ZZ b;
    ZZ w;
    ZZ q;               // Modulus in the ciphertext space
    ZZ t;               // Modulus in the plaintext space
    tuple<ZZX, ZZX> pk; // Public key tuple of the class
    Vec<tuple<ZZX, ZZX>> evk;

    /**
     * @brief Construct a new Fan-Vercauteren Scheme object initializing the modulo, the secret key, the public key
     * and the evaluation keys.
     */
    FanVercScheme(int n, int b_, ZZ q, int w_)
    {
        randomness::setSeed();
        this->from = ZZ(-(q - 1) / 2);
        this->to = ZZ((q - 1) / 2);
        this->n = n;
        this->b = ZZ(b_);
        this->w = ZZ(w_);
        this->q = q;
        this->t = t;
        this->l = (int)floor(log(this->q) / log(this->w));
        ZZ_p::init(this->q);

        printf("Setting n as %i and b as %i and q as %i\n", this->n, conv<int>(this->b), conv<int>(this->q));

        setDelta();
        setIrreduciblePoly();

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
    ZZX fractionalEncoder(Rational message)
    {

        NTL::ZZX encoded = fractionalEncoder::HenselCode::PolyEncode(conv<long>(this->b), this->n, message);
        Rational decoded = fractionalEncoder::HenselCode::PolyDecode(conv<long>(this->b), this->n, encoded);

        return encoded;
    }

    tuple<ZZX, ZZX> encrypt(Rational message_)
    {
        if (DEBUG)
        {
            printf("Encoding %s...\n", message_.ToString().c_str());
        }
        ZZX message = fractionalEncoder(message_);
        if (DEBUG)
        {
            cout << "message = " << logPoly(message) << "\n";
            cout << "message = " << message << "\n";
        }
        ZZX e0, e1, u, c0, c1, pk0, pk1;

        u = randomness::RandomPolynomial(-1, 1, this->n);
        e0 = randomness::NormalPolynomial(this->n);
        e1 = randomness::NormalPolynomial(this->n);

        tie(pk0, pk1) = this->pk;
        if (DEBUG)
        {
            cout << "u = " << logPoly(u) << "\n";
            cout << "e0 = " << logPoly(e0) << "\n";
            cout << "e1 = " << logPoly(e1) << "\n";
            cout << "delta = " << logPoly(this->delta) << "\n";
        }

        c0 = MulMod(message, this->delta, this->irreduciblePolynomial) // ∆ * m
             + MulMod(pk0, u, this->irreduciblePolynomial)             // p0 * u
             + e0;
        c1 = MulMod(pk1, u, this->irreduciblePolynomial) // p1 * u
             + e1;
        resizePolynomialsCoeff(&c0, this->q);
        resizePolynomialsCoeff(&c1, this->q);
        if (DEBUG)
        {
            cout << "c0 = " << logPoly(c0) << "\n";
            cout << "c1 = " << logPoly(c1) << "\n";
        }
        return make_tuple(c0, c1);
    }

    Rational decrypt(tuple<ZZX, ZZX> ct)
    {
        ZZX c0, c1, M;
        ZZX temp;
        temp.SetLength(2);
        temp[0] = -this->b;
        temp[1] = ZZ(1);
        tie(c0, c1) = ct;
        if (DEBUG)
        {
            cout << "Decoding...\n\t" << logPoly(c0) << "\n\t" << logPoly(c1) << "\n";
        }
        M = MulMod(c1, this->sk, this->irreduciblePolynomial) + c0; // c0 + c1 * s
        resizePolynomialsCoeff(&M, this->q);
        M = MulMod(M, temp, this->irreduciblePolynomial);
        if (DEBUG)
        {
            cout << "M before " << logPoly(M) << "\n";
        }
        RR q, coeff;
        conv(q, this->q);
        for (int i = 0; i < this->n; i++)
        {
            conv(coeff, M[i]);
            M[i] = conv<ZZ>(round(coeff / q));
        }
        if (DEBUG)
        {
            cout << "M = " << M << "\n";
        }
        Rational result = fractionalEncoder::HenselCode::PolyDecode(conv<long>(this->b), this->n, M);
        return result;
    }

    tuple<ZZX, ZZX> add(tuple<ZZX, ZZX> ct0, tuple<ZZX, ZZX> ct1)
    {
        ZZX ct00, ct01, ct10, ct11, ct20, ct21;
        tie(ct00, ct01) = ct0;
        tie(ct10, ct11) = ct1;

        ct20 = ct00 + ct10;
        ct21 = ct01 + ct11;
        tools::resizePolynomialsCoeff(&ct20, this->q);
        tools::resizePolynomialsCoeff(&ct21, this->q);

        return make_tuple(ct20, ct21);
    }

    tuple<ZZX, ZZX> multiply(tuple<ZZX, ZZX> c0, tuple<ZZX, ZZX> c1) {
        ZZX c00, c01, c10, c11, c20, c21, c22, a;
        tie(c00, c01) = c0;
        tie(c10, c11) = c1;
        a.SetLength(2);
        a[0] = this->b;
        a[1] = 1;

        c20 = MulMod(MulMod(c00, c10, this->irreduciblePolynomial), a, this->irreduciblePolynomial);
        c21 = MulMod(
            MulMod(c00, c11, this->irreduciblePolynomial) +  // c00 * c11
            MulMod(c01, c10, this->irreduciblePolynomial),   // c01 * c10
            a, this->irreduciblePolynomial);
        c22 = MulMod(MulMod(c01, c11, this->irreduciblePolynomial), a, this->irreduciblePolynomial);

        RR q, coeff;
        conv(q, this->q);
        for (int i = 0; i < this->n; i++)
        {
            conv(coeff, c20[i]);
            c20[i] = conv<ZZ>(round(coeff / q));
            conv(coeff, c21[i]);
            c21[i] = conv<ZZ>(round(coeff / q));
            conv(coeff, c22[i]);
            c22[i] = conv<ZZ>(round(coeff / q));
        }

        tools::resizePolynomialsCoeff(&c20, this->q);
        tools::resizePolynomialsCoeff(&c21, this->q);
        tools::resizePolynomialsCoeff(&c22, this->q);

        return relinearize(c20, c21, c22);

    }

    tuple<ZZX, ZZX> relinearize(ZZX c0, ZZX c1, ZZX c2) {
        ZZX ct0, ct1, evk0, evk1, evkSum0, evkSum1;
        evkSum0.SetLength(this->n);
        evkSum1.SetLength(this->n);
        for (int i = 0; i < this->n; i++)
        {
            evkSum0[i] = 0;
            evkSum1[i] = 0;
        }
        
        ct0 = c0;
        ct1 = c1;
        for (int i = 0; i < this->l + 1; i++) {
            tie(evk0, evk1) = this->evk[i];
            ct0 = ct0 + MulMod(evk0, c2, this->irreduciblePolynomial);
            ct1 = ct1 + MulMod(evk1, c2, this->irreduciblePolynomial);
            tools::resizePolynomialsCoeff(&evkSum0, this->w);
            tools::resizePolynomialsCoeff(&evkSum1, this->w);
        }

        return make_tuple(ct0, ct1);
    }

    // #region Private functions

private:
    /**
     * @brief Generates a random polynomial with random coeficients (-1, 0, 1)
     *
     */
    void secretKeyGen()
    {
        this->sk = randomness::RandomPolynomial(-1, 1, this->n);
        if (DEBUG)
        {
            cout << "sk = " << logPoly(sk) << "\n";
        }
        printf("Private key generated....\n\n");
    }

    /**
     * @brief Generates and stores the public key of the cipher
     *
     * @param sk The secret key required to generate the pubblic key
     */
    void publicKeyGen()
    {
        ZZX pk1 = randomness::RandomPolynomial(this->from, this->to, this->n), pk0;
        ZZX e = randomness::NormalPolynomial(this->n);

        if (DEBUG)
        {
            cout << "e = " << logPoly(e) << "\n";
            cout << "pk1 = " << logPoly(pk1) << "\n";
        }
        pk0 = -(MulMod(pk1, this->sk, this->irreduciblePolynomial) + e);
        resizePolynomialsCoeff(&pk0, this->q);
        if (DEBUG)
        {
            cout << "pk0 = " << logPoly(pk0) << "\n";
        }

        this->pk = make_tuple(pk0, pk1);
        printf("Public keys generated....\n\n");
    }

    void evaluationKeyGen()
    {
        Vec<tuple<ZZX, ZZX>> pairs;
        pairs.SetLength(l + 1);
        const int length = pairs.length();

        NTL::ZZX wi, spow2, e, y;
        wi.SetLength(1);
        spow2 = MulMod(this->sk, this->sk, this->irreduciblePolynomial);
        for (int i = 0; i < length; i++)
        {
            ZZX x = randomness::RandomPolynomial(this->from, this->to, this->n); // a
            wi[0] = NTL::ZZ(NTL::power_ZZ(conv<long>(this->w), i));
            e = randomness::NormalPolynomial(this->n);
            y = -(
                    MulMod(x, this->sk, this->irreduciblePolynomial) + e) // ai * s + e
                + MulMod(wi, spow2, this->irreduciblePolynomial           // w^i * s^2
                  );
            resizePolynomialsCoeff(&y, this->q);

            // cout << "i = " << i << "\n";
            // cout << "y = " << y << "\n";
            // cout << "x = " << x << "\n";
            pairs[(length - i) - 1] = make_tuple(y, x);
        }
        this->evk = pairs;
        printf("Evaluation keys generated....\n\n");
    }

    void resizePolynomialsCoeff(ZZX *polynomial, ZZ mod)
    {
        tools::resizePolynomialsCoeff(polynomial, mod);
    }

    void setIrreduciblePoly()
    {
        this->irreduciblePolynomial.SetLength(n + 1);
        this->irreduciblePolynomial[0] = NTL::ZZ(1);
        this->irreduciblePolynomial[n] = NTL::ZZ(1);
    }

    void setDelta()
    {
        // Calculating ∆b
        RR delta_multiplicative = -conv<RR>(this->q) / conv<RR>(power(conv<RR>(this->b), this->n) + 1);
        this->delta.SetLength(this->n);
        for (int i = 0; i < this->n; i++)
        {
            this->delta[i] = conv<ZZ>(round(delta_multiplicative * power(conv<RR>(this->b), this->n - (i + 1))));
        }
    }
    // #endregion
};

int main(int argc, char *argv[])
{
    ZZ q = power2_ZZ(100);
    tuple<ZZX, ZZX> ct, c0, c1, c2, csum, cmul;
    int n = 4;
    int b = 128;
    int w = 4;
    FanVercScheme f(n, b, q, w);

    if (argc >= 3)
    {
        Rational message = Rational(ZZ(atoi(argv[1])), ZZ(atoi(argv[2])));
        ct = f.encrypt(message);
        Rational dt = f.decrypt(ct);
        printf("Message %s\n", dt.ToString().c_str());
        randomness::clear();
    }
    else
    {
        printf("Launching test run...\n\n");

        Rational m0 = Rational(1, 3);
        Rational m1 = Rational(1, 3);
        Rational m2 = Rational(2, 3);

        c0 = f.encrypt(m0);
        Rational r0 = f.decrypt(c0);
        printf("Decoded m0 in %s\n", r0.ToString().c_str());

        c1 = f.encrypt(m1);
        Rational r1 = f.decrypt(c1);
        printf("Decoded m1 in %s\n", r1.ToString().c_str());

        csum = f.add(c0, c1);
        Rational rsum = f.decrypt(csum);
        printf("Result sum is %s\n", rsum.ToString().c_str());

        cmul = f.multiply(c0, c1);
        Rational rmul = f.decrypt(cmul);
        printf("Result mul is %s\n", rmul.ToString().c_str());

        
        c2 = f.encrypt(m2);
        Rational r2 = f.decrypt(c2);
        printf("Decoded m2 in %s\n", r2.ToString().c_str());

        ZZX cs10, cs11, c20, c21;
        tie(cs10, cs11) = csum;
        tie(c20, c21) = c2;
        cout << "CS" << cs10 << cs11 << "\n";
        cout << "C2" << c20 << c21 << "\n";
    }

    return 0;
}
