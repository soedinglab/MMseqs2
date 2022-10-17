#ifndef BYTE_PARSER_H
#define BYTE_PARSER_H

#include <cinttypes>
#include <cerrno>
#include <string>
#include "Debug.h"

class ByteParser {
public:    
    static uint64_t parse(const std::string& sizeAndUnit) {
        // default unit is M
        uint64_t unitFactor = TWO_POW_10 * TWO_POW_10;

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
        
        char* rest;
        errno = 0;
        uint64_t size = strtoull(digitsString.c_str(), &rest, 10);
        if ((rest != digitsString.c_str() && *rest != '\0') || errno == ERANGE) {
            // conversion to uint64_t failed
            return INVALID_SIZE;
        }
        uint64_t sizeBits = highestOneBitPosition(size);
        uint64_t unitFactorBits = highestOneBitPosition(unitFactor);

        if ((sizeBits + unitFactorBits) > 64) {
            // cannot store (size * unitFactor) in a uint64_t
            return INVALID_SIZE;
        }
        uint64_t numBytes = (size * unitFactor);
        return numBytes;
    };
    
    static std::string format(uint64_t numBytes, char unit='a', char accuracy='l') {
        uint64_t unitT = TWO_POW_10 * TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
        uint64_t unitG = TWO_POW_10 * TWO_POW_10 * TWO_POW_10;
        uint64_t unitM = TWO_POW_10 * TWO_POW_10;
        uint64_t unitK = TWO_POW_10;

        // in default mode (l), 1,433,600 will be rounded to 1M.
        // in more informative mode (h), 1,433,600 will be formatted to 1400K.
        uint64_t valForModCheck = 0;
        if (accuracy != 'l') {
            valForModCheck = numBytes;
        }
        
        if (unit == 'a') {
            // auto-detect the unit to use:
            if ((numBytes / unitT > 0) && (valForModCheck % unitT == 0)) {
                unit = 'T';
            } else if ((numBytes / unitG > 0) && (valForModCheck % unitG == 0)) {
                unit = 'G';
            } else if ((numBytes / unitM > 0) && (valForModCheck % unitM == 0)) {
                unit = 'M';
            } else if ((numBytes / unitK > 0) && (valForModCheck % unitK == 0)) {
                unit = 'K';
            } else {
                unit = 'B';
            }
        }

        uint64_t unitFactor = 1;
        if ((unit == 't') || (unit == 'T')) {
            unitFactor = unitT;
        } else if ((unit == 'g') || (unit == 'G')) {
            unitFactor = unitG;
        } else if ((unit == 'm') || (unit == 'M')) {
            unitFactor = unitM;
        } else if ((unit == 'k') || (unit == 'K')) {
            unitFactor = unitK;
        } else if ((unit == 'b') || (unit == 'B')) {
            unitFactor = 1;
        } else {
            // unrecognized unit
            Debug(Debug::ERROR) << "Invalid unit " << unit << " for format conversion given\n";
            EXIT(EXIT_FAILURE);
        }

        uint64_t value = (uint64_t)(numBytes / unitFactor);
        std::string str(SSTR(value));
        if (value > 0) {
            str.append(1, unit);
        }
        return str;
    };

    static const uint64_t INVALID_SIZE = UINT32_MAX;
private:
    static uint64_t highestOneBitPosition(uint64_t number) {
        uint64_t bits = 0;
        while (number != 0) {
            bits++;
            number >>= 1;
        }
        return bits;
    };

    static const uint64_t TWO_POW_10 = 1024;
};

#endif
