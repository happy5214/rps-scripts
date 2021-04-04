/* Aliquot prime power trial factoring program.
 *
 * This program verifies the divisibility and abundance of a set of divisors in
 * relation to a given prime base and exponent pair.
 *
 * (C) Alexander Jones, 2021. My own code is under the MIT License, which is
 * included in this repository.
 *
 * Adapted from code by Oliver Kruse, available at
 * https://www.mersenneforum.org/showpost.php?p=574965&postcount=1052.
 */

#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <gmpxx.h>
#include <primesieve.hpp>

using namespace std;

// Factor typedef
typedef pair<mpz_class, uint64_t> Factor;

// Convert Factor typedef to string
string factorToString(Factor & factor) {
    if (factor.second == 1) {
        return factor.first.get_str();
    } else {
        return factor.first.get_str() + "^" + to_string(factor.second);
    }
}

// Write string for factor vector
string getBaseFactorString(vector<Factor> & factors) {
    ostringstream ss;
    for (vector<Factor>::size_type i = 0; i < factors.size(); i++) {
        if (i) {
            ss << " * ";
        }
        ss << factorToString(factors[i]);
    }
    return ss.str();
}

// Push new factor to vector
void foundFactor(vector<Factor> & factors, mpz_class factor, uint64_t exponent = 1) {
    factors.push_back(make_pair(factor, exponent));
}

// Validate a number string
bool isnumber(string s) {
    for (string::size_type j = 0; j < s.size(); ++j) {
        if (s[j] < '0' || s[j] > '9') return false;
    }
    return true;
}

// Parse a factor (with exponent) string and add it to the factor vector
void parseExponent(vector<Factor> & factors, string exponentString) {
    factors.clear();
    string s;
    stringstream ss(exponentString);
    while (ss >> ws >> s) {
        if (s == "*") continue;
        string::size_type o;
        if ((o = s.find('^')) != s.npos) {
            string p = s.substr(0, o);
            string e = s.substr(o + 1);
            if (!isnumber(p) || !isnumber(e)) {
                cerr << "not a number: " << s << endl;
                exit(1);
            }
            foundFactor(factors, mpz_class(p), stoi(e));
        } else {
            if (!isnumber(s)) {
                cerr << "not a number: " << s << endl;
                exit(1);
            }
            foundFactor(factors, mpz_class(s));
        }
    }
}

// Sort the factor vector and merge common factors
void merge_factors(vector<Factor> & factors) {
    sort(factors.begin(), factors.end());
    for (vector<Factor>::size_type j = 1; j < factors.size();) {
        if (factors[j].first == factors[j - 1].first) {
            factors[j - 1].second += factors[j].second;
            factors.erase(factors.begin() + j);
        } else {
            ++j;
        }
    }
}

// Perform simple trial factoring
void simpleFactor(mpz_class n, vector<Factor> & factors, uint64_t factoringLimit) {
    factors.clear();
    if (mpz_probab_prime_p(n.get_mpz_t(), 25)) {
        foundFactor(factors, n);
        return;
    }

    primesieve::iterator it;
    uint64_t prime = it.next_prime();

    bool fullyFactored = false;
    for (; prime < factoringLimit; prime = it.next_prime()) {
        while (mpz_divisible_ui_p(n.get_mpz_t(), prime)) {
            mpz_class tp(prime);
            foundFactor(factors, tp);
            mpz_divexact_ui(n.get_mpz_t(), n.get_mpz_t(), prime);
        }
        if (mpz_cmp_ui(n.get_mpz_t(), 1) == 0) {
            break;
        } else if (mpz_probab_prime_p(n.get_mpz_t(), 25)) {
            foundFactor(factors, n);
            break;
        }
    }

    merge_factors(factors); //sorts factors and merges any <p,x>,<p,y> into <p,x+y>
}

// Perform simple trial factoring
uint64_t fullFactor(mpz_class & base, vector<Factor> & baseFactors, mpz_class & exponent, uint64_t factoringLimit, vector<Factor> & factors) {
    factors.clear();

    mpz_class divisor = 1;
    mpz_class exponentPlusOne = exponent + 1;
    mpz_class tmp;

    uint64_t totalFactorCount = 0;

    vector<Factor>::size_type i;
    uint64_t j;

    for (i = 0; i < baseFactors.size(); i++) {
        mpz_class factorMinusOne = baseFactors[i].first - 1;
        for (j = 0; j < baseFactors[i].second; j++) {
            divisor *= factorMinusOne;
        }
    }

    primesieve::iterator it;
    uint64_t prime = it.next_prime();

    bool fullyFactored = false;
    for (; prime < factoringLimit; prime = it.next_prime()) {
        mpz_class candidate(prime);
        bool doesDivideThisPrime;
        int divideAmount = -1;
        do {
            mpz_class firstAddend = 1;
            mpz_class modulus = divisor * candidate;
            for (i = 0; i < baseFactors.size(); i++) {
                mpz_class div = baseFactors[i].first;
                for (j = 0; j < baseFactors[i].second; j++) {
                    mpz_powm(tmp.get_mpz_t(), div.get_mpz_t(), exponentPlusOne.get_mpz_t(), modulus.get_mpz_t());
                    if (mpz_cmp_ui(tmp.get_mpz_t(), 0) == 0) { // we want to avoid negative numbers
                        tmp = modulus - 1;
                    } else if (mpz_cmp_ui(tmp.get_mpz_t(), 1) == 0) { // special case: in this case, we can break out early, since the product will be 0 if a single factor is 0
                        firstAddend = 0;
                        break;
                    } else {
                        tmp--;
                    }
                    firstAddend *= tmp;
                    // for bases with tons of factors, we could do
                    //     firstAddend %= modulus;
                    // here, at least sometimes
                }
            }
            mpz_powm(tmp.get_mpz_t(), base.get_mpz_t(), exponent.get_mpz_t(), modulus.get_mpz_t());
            mpz_class secondAddend = modulus - tmp * divisor;
            mpz_class sum = firstAddend + secondAddend;
            sum %= modulus;
            doesDivideThisPrime = mpz_cmp_ui(sum.get_mpz_t(), 0) == 0;
            divideAmount++; // this will also be increased if the division is unsuccessful, that's why we start with -1
            candidate *= prime;
        } while (doesDivideThisPrime); // the number could divide n and n² and n³...

        if (divideAmount > 0) {
            foundFactor(factors, mpz_class(prime), divideAmount);
            totalFactorCount += divideAmount;
        }
    }

    return totalFactorCount;
}

// Multiply out the factor vector
void multiply(vector<Factor> & factors, mpz_class & n) {
    mpz_class tmp;
    n = 1;
    for (vector<Factor>::size_type j = 0; j < factors.size(); ++j) {
        mpz_pow_ui(tmp.get_mpz_t(), factors[j].first.get_mpz_t(), factors[j].second);
        n *= tmp;
    }
}

// Calculate the sum of divisors and product of the factor vector
void sigma(vector<Factor> & factors, mpz_class & s, mpz_class & n) {
    mpz_class t, tmp, tmp2;
    n = 1;
    s = 1;
    for (vector<Factor>::size_type j = 0; j < factors.size(); ++j) {
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

// Print help
void print_help() {
    cout << "usage: verifyPrimePowerAbundance <base> <exponent> [<limit>]" << endl;
}

#define DEFAULT_TF_LIMIT 100000

int main(int argc, char ** argv) {
    // Validate argument count
    if (argc < 3) {
        print_help();
        return 1;
    }

    // Load command-line arguments
    mpz_class base;
    base.set_str(argv[1], 10);

    vector<Factor> exponentFactors;
    parseExponent(exponentFactors, argv[2]);
    mpz_class exponent;
    multiply(exponentFactors, exponent);

    uint64_t factoringLimit;
    if (argc > 3) {
        factoringLimit = atol(argv[3]);
    } else {
        factoringLimit = DEFAULT_TF_LIMIT;
    }

    vector<Factor> baseFactors;
    simpleFactor(base, baseFactors, factoringLimit);

    vector<Factor> resultFactors;
    uint64_t totalFactorCount = fullFactor(base, baseFactors, exponent, factoringLimit, resultFactors);
    vector<Factor>::size_type uniqueFactorCount = resultFactors.size();

    if (resultFactors.empty()) {
        cout << "No factors found up to given limit." << endl;
    } else {
        string resultFactorString = getBaseFactorString(resultFactors);
        cout << "d = " << resultFactorString << " * remainder up to limit=" << factoringLimit << endl;
    }

    return 0;
}
