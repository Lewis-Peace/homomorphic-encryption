#include <random>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <time.h>

// TODO: implement the discrete gaussian distribution
class DiscreteGaussianNumber {
    gsl_rng* r;

    public:

        DiscreteGaussianNumber() {
            this->r = gsl_rng_alloc(gsl_rng_mt19937);
            gsl_rng_set(this->r, time(0));
        }

        double getNumber() {
            return floor(gsl_ran_gaussian(r, 3.19) + 0.5);
        }
};