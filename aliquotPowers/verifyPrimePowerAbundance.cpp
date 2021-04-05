/* Aliquot prime power abundance verifier.
 *
 * This program verifies the divisibility and abundance of a set of divisors in
 * relation to a given prime base and exponent pair.
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

// adds <factor> to <factors>

void found_factor(mpz_class & factor, FactorVector & factors, int exponent = 1) {
    factors.push_back(make_pair(factor, exponent));
}

void load_factors(FactorVector & factors) {
    ifstream factorFile("partial_factors");
    if (!factorFile.is_open()) {
        cout << "WARNING: couldn't open input file for reading!" << endl;
        exit(1);
    }
    string line;
    while (getline(factorFile, line)) {
        mpz_class factor;
        int caretIndex = line.find("^");
        int power;
        if (caretIndex == string::npos) {
            factor.set_str(line, 10);
            power = 1;
        } else {
            factor.set_str(line.substr(0, caretIndex), 10);
            power = stoi(line.substr(caretIndex + 1));
        }
        found_factor(factor, factors, power);
    }
    factorFile.close();
}

//calculates <n>=product(<factors>) and <s>=sigma(<n>)

void sigma(FactorVector & factors, mpz_class & s, mpz_class & n) {
    mpz_class t, tmp, tmp2;
    n = 1;
    s = 1;
    for (FactorVector::size_type j = 0; j < factors.size(); ++j) {
        if (factors[j].second == 1) {
            t = factors[j].first + 1;
        } else {
            t = 1;
            mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second + 1);
            tmp -= 1;
            tmp2 = factors[j].first - 1;
            mpz_divexact(t.get_mpz_t(), tmp.get_mpz_t(), tmp2.get_mpz_t());
        }

        s *= t;

        mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second);
        n *= tmp;
    }
}

void print_help() {
    cout << "usage: verifyPrimePowerAbundance <base> <exponent>" << endl
         << "Place partial factorization (one factor per line) in file 'partial_factors'" << endl;
}

int main(int argc, char ** argv) {
    if (argc < 3) {
        print_help();
        return 1;
    }

    mpz_class base;
    base.set_str(argv[1], 10);
    long exponent = atol(argv[2]);

    FactorVector factors; //vector<pair<p_i,x_i> >, n = product(p_i^x_i)
    load_factors(factors);

    // Validate abundance
    mpz_class n, s, partial;
    sigma(factors, s, partial); //calculate sigma(n) and partial = product(factors)
    n = s - partial;
    mpq_class abundance(n, partial);
    if (n > partial) {
        cout << "Index 1 of " << base.get_str() << "^" << exponent << " is abundant! (" << abundance.get_d() << ")" << endl;
    } else {
        cout << "Index 1 of " << base.get_str() << "^" << exponent << " is not abundant. (" << abundance.get_d() << ")" << endl;
    }

    // Validate divisibility
    mpz_class base_modulo = base - 1;
    mpz_class modulo;
    mpz_class result;
    bool allFactorsDivide = true;
    for (FactorVector::size_type i = 0; i < factors.size(); ++i) {
        mpz_class factor;
        mpz_pow_ui(factor.get_mpz_t(), factors[i].first.get_mpz_t(), factors[i].second);
        modulo = base_modulo * factor;
        mpz_powm_ui(result.get_mpz_t(), base.get_mpz_t(), exponent, modulo.get_mpz_t());
        if (mpz_cmp_ui(result.get_mpz_t(), 1)) {
            if (factors[i].second > 1) {
                cout << "Does not divide: " << factors[i].first.get_str() << "^" << factors[i].second << endl;
            } else {
                cout << "Does not divide: " << factor.get_str() << endl;
            }
            allFactorsDivide = false;
        }
    }
    if (allFactorsDivide) {
        cout << "All factors divide!" << endl;
    }
    return 0;
}
