#ifndef BYTE_PARSER_H
#define BYTE_PARSER_H

#include "Debug.h"
#include <string>

class ByteParser {
public:    
    static size_t parse(const std::string& sizeAndUnit) {
        // default unit is M
        size_t unitFactor = TWO_POW_10 * TWO_POW_10;
        size_t size = 0;

        size_t strLen = sizeAndUnit.size();
        char lastChar = sizeAndUnit[strLen - 1];
        std::string digitsString;

        if (std::isdigit(lastChar)) {
            // all digits - default unit
            digitsString = sizeAndUnit.substr(0);
        } else {
            digitsString = sizeAndUnit.substr(0,(strLen - 1));

            if ((lastChar == 't') || (lastChar == 'T')) {
                unitFactor = TWO_POW_10 * TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
            } else if ((lastChar == 'g') || (lastChar == 'G')) {
                unitFactor = TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
            } else if ((lastChar == 'm') || (lastChar == 'M')) {
                unitFactor = TWO_POW_10 * TWO_POW_10;
            } else if ((lastChar == 'k') || (lastChar == 'K')) {
                unitFactor = TWO_POW_10;
            } else if ((lastChar == 'b') || (lastChar == 'B')) {
                unitFactor = 1;
            } else {
                // unrecognized unit...
                return INVALID_SIZE;
            }
        }
        
        // convert to size_t
        if (1 == sscanf(digitsString.c_str(), "%zu", &size)) {
            size_t sizeBits = highestOneBitPosition(size);
            size_t unitFactorBits = highestOneBitPosition(unitFactor);

            if ((sizeBits + unitFactorBits) > 64) {
                // cannot store (size * unitFactor) in a size_t
                return INVALID_SIZE;
            }
            size_t numBytes = (size * unitFactor);
            return numBytes;
        } else {
            // conversion failed
            return INVALID_SIZE;
        }
    };
    
    static std::string format(size_t numBytes, char unit='a') {
        if (unit == 'a') {
            // auto-detect the unit to use:
            if ((numBytes / (TWO_POW_10 * TWO_POW_10 * TWO_POW_10 * TWO_POW_10)) > 0) {
                unit = 'T';
            } else if ((numBytes / (TWO_POW_10 * TWO_POW_10 * TWO_POW_10)) > 0) {
                unit = 'G';
            } else if ((numBytes / (TWO_POW_10 * TWO_POW_10)) > 0) {
                unit = 'M';
            } else if ((numBytes / TWO_POW_10) > 0) {
                unit = 'K';
            } else {
                unit = 'B';
            }
        }

        size_t unitFactor = 1;
        if ((unit == 't') || (unit == 'T')) {
            unitFactor = TWO_POW_10 * TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
        } else if ((unit == 'g') || (unit == 'G')) {
            unitFactor = TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
        } else if ((unit == 'm') || (unit == 'M')) {
            unitFactor = TWO_POW_10 * TWO_POW_10;
        } else if ((unit == 'k') || (unit == 'K')) {
            unitFactor = TWO_POW_10;
        } else if ((unit == 'b') || (unit == 'B')) {
            unitFactor = 1;
        } else {
            // unrecognized unit
            Debug(Debug::ERROR) << "Invalid unit " << unit << " for format conversion given\n";
            EXIT(EXIT_FAILURE);
        }

        size_t value = (size_t)(numBytes / unitFactor);
        std::string str(SSTR(value));
        if (value > 0) {
            str.append(1, unit);
        }
        return str;
    };

    static const size_t INVALID_SIZE = SIZE_MAX;
private:
    static size_t highestOneBitPosition(size_t number) {
        size_t bits = 0;
        while (number != 0) {
            bits++;
            number >>= 1;
        };
        return bits;
    };

    static const size_t TWO_POW_10 = 1024;
};

#endif
