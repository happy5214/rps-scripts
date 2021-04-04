/* Aliquot integer power abundance calculator.
 *
 * This program scans an exponent range for a given base for an aliquot sum
 * which can be determined to be abundant based on its factors up to 10^4. It
 * will output any such exponents, along with the factors that establish
 * abundance, to a log file.
 *
 * (C) Alexander Jones, 2021. My code is under the MIT License, which is
 * included in this repository.
 *
 * Adapted from aliqueit by Mikael Klasson, Christian Beer, et al., which can
 * be found at https://github.com/ChristianBeer/aliqueit. My understanding of
 * its rather unique license agreement is that it is public domain with a
 * disclaimer of warranty.
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <gmpxx.h>

using namespace std;

typedef vector<pair<mpz_class, int> > FactorVector;

vector<unsigned int> trial_primes; //precalced primes for trial factoring

void factor(mpz_class n, FactorVector & factors);

void precalc_trial_primes() {
    cout << "Precalcing primes for trial factoring..." << endl;
    mpz_class p(2);
    while (p < 10000) {
        trial_primes.push_back(p.get_ui());
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    }
}

//prints a log msg and adds <factor> to <factors>

void found_factor(mpz_class & factor, FactorVector & factors) {
    factors.push_back(make_pair(factor, 1));
}

//sorts factors and merges any <p,x>,<p,y> into <p,x+y>

void merge_factors(FactorVector & factors) {
    sort(factors.begin(), factors.end());
    for (FactorVector::size_type j = 1; j < factors.size();) {
        if (factors[j].first == factors[j - 1].first) {
            factors[j - 1].second += factors[j].second;
            factors.erase(factors.begin() + j);
        } else {
            ++j;
        }
    }
}

//factors <n> and returns its prime components (with exponents) in <factors>.
//vector<pair<p_i,x_i> >, n = product(p_i^x_i)

void factor(mpz_class n, FactorVector & factors) {
    factors.clear();

    for (vector<unsigned int>::size_type j = 0; j < trial_primes.size(); ++j) {
        while (mpz_divisible_ui_p(n.get_mpz_t(), trial_primes[j])) {
            mpz_class tp(trial_primes[j]);
            found_factor(tp, factors);
            mpz_divexact_ui(n.get_mpz_t(), n.get_mpz_t(), trial_primes[j]);
        }
    }

    merge_factors(factors); //sorts factors and merges any <p,x>,<p,y> into <p,x+y>
}

//calculates <n>=product(<factors>) and <s>=sigma(<n>)

void sigma(FactorVector & factors, mpz_class & s, mpz_class & n) {
    mpz_class t, tmp, tmp2;
    n = 1;
    s = 1;
    for (FactorVector::size_type j = 0; j < factors.size(); ++j) {
        t = 1;
        mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second + 1);
        tmp -= 1;
        tmp2 = factors[j].first - 1;
        mpz_divexact(t.get_mpz_t(), tmp.get_mpz_t(), tmp2.get_mpz_t());

        s *= t;

        mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second);
        n *= tmp;
    }
}

void print_help() {
    cout << "usage: powerAbundance <base> <minExp> <maxExp> [<skip>]" << endl;
}

int main(int argc, char ** argv) {
    if (argc < 4) {
        print_help();
        return 1;
    }

    mpz_class base;
    base.set_str(argv[1], 10);
    int min = atoi(argv[2]);
    int max = atoi(argv[3]);
    int skip;
    if (argc >= 5) {
        skip = atoi(argv[4]);
    } else {
        skip = 2;
    }

    precalc_trial_primes();

    FactorVector factors; //vector<pair<p_i,x_i> >, n = product(p_i^x_i)
    mpz_class n, s, partial;

    FactorVector::size_type j;
    for (int i = min; i <= max; i += skip) {
        // Index 0 -> 1
        factor(base, factors);
        for (j = 0; j < factors.size(); ++j) {
            factors[j].second *= i;
        }
        sigma(factors, s, partial); //calculate sigma(n) and partial = product(factors)
        n = s - partial;

        // Index 1 -> 2
        factor(n, factors);
        sigma(factors, s, partial); //calculate sigma(n) and partial = product(factors)
        n = s - partial;

        // Check for abundance of partial factors.
        if (n > partial) {
            ofstream fff("power_abundant_exponents", ios::app);
            if (!fff.is_open()) {
                cout << "WARNING: couldn't open output file for writing!" << endl;
                return 1;
            }
            string factorization = "";
            for (j = 0; j < factors.size(); ++j) {
                if (j) factorization += " * ";
                factorization += factors[j].first.get_str() + (factors[j].second > 1 ? ("^" + to_string(factors[j].second)) : "");
            }
            fff << base.get_str() << " " << i << " (" << factorization << ")" << endl;
            fff.close();
            cout << base.get_str() << "^" << i << " is abundant!" << endl;
        }
    }
    return 0;
}
