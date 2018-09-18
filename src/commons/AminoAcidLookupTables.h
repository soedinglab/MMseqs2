#ifndef AMINOACIDLOOKUP_H
#define AMINOACIDLOOKUP_H

#include <unordered_map>


class Doolittle {
public:
    std::unordered_map<char, float> values;

    Doolittle() {
        values['a'] = 6.3;
        values['r'] = 0.0;
        values['n'] = 1.0;
        values['d'] = 1.0;
        values['c'] = 7.0;
        values['q'] = 1.0;
        values['e'] = 1.0;
        values['g'] = 4.1;
        values['h'] = 1.3;
        values['i'] = 9.0;
        values['l'] = 5.2;
        values['k'] = 0.6;
        values['m'] = 6.4;
        values['f'] = 7.2;
        values['p'] = 2.9;
        values['s'] = 3.6;
        values['t'] = 3.8;
        values['w'] = 3.6;
        values['y'] = 3.2;
        values['v'] = 8.7;
        values['x'] = 0.0;
        values['0'] = 0.0; // N-ter
        values['1'] = 0.0; // C-ter
    }
};



class Charges {
public:
    std::unordered_map<char, float> values;

    Charges() {
        const float pH = 7.0;

        std::unordered_map<char, float> pKs, chargeSign;

        // pKs values:
        // Bjellqvist Dawson EMBOSS Lehninger Murray Rodwell Sillero Solomon Stryer
        pKs['c'] = 9.00;//    8.3    8.5      8.18   8.33    8.33     9.0     8.3    8.5
        pKs['d'] = 4.05;//    3.9    3.9      3.65   3.68    3.86     4.0     3.9    4.4
        pKs['e'] = 4.45;//    4.3    4.1      4.25   4.25    4.25     4.5     4.3    4.4
        pKs['h'] = 5.98;//    6.0    6.5      6.00   6.00    6.00     6.4     6.0    6.5
        pKs['k'] = 10.00;//   10.5   10.8     10.53  11.50   11.50    10.4    10.5   10.0
        pKs['r'] = 12.00;//   12.0   12.5     12.48  11.50   11.50    12.0    12.5   12.0
        pKs['y'] = 10.00;//   10.1   10.1     10.07  10.07   10.70    10.0    10.1   10.0
        pKs['1'] = 3.55;//    3.2    3.6      2.34   2.15    3.10     3.2     3.2    3.2 // C ter
        pKs['0'] = 7.50;//    8.2    8.6      9.69   9.52    8.00     8.2     8.2    8.2 // N ter

        chargeSign['c'] = -1.0f;
        chargeSign['d'] = -1.0f;
        chargeSign['e'] = -1.0f;
        chargeSign['y'] = -1.0f;
        chargeSign['h'] = 1.0f;
        chargeSign['k'] = 1.0f;
        chargeSign['r'] = 1.0f;
        chargeSign['1'] = -1.0f; // C ter
        chargeSign['0'] = 1.0f; // N ter

        for (std::unordered_map<char, float>::iterator k = pKs.begin(); k != pKs.end(); k++) {
            values[k->first] = chargeSign[k->first] / (1 + pow(10, (chargeSign[k->first] * (pH - pKs[k->first]))));
        }
    }
};


#endif