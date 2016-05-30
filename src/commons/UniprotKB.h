//
// Created by Milot on 23/05/16.
//

#ifndef MMSEQS_UNIPROTKB_H
#define MMSEQS_UNIPROTKB_H

class UniprotKB {
public:
    UniprotKB() : dbColumns(17), isInEntry(false), hasEntry(false) {
        streams = new std::ostringstream[dbColumns];
    };

    ~UniprotKB() {
        delete[] streams;
    }

    size_t getColumnCount() const {
        return dbColumns;
    }

    bool readLine (const char* line);
    std::string getColumn (size_t column);

    enum {
        COL_KB_ID = 0,
        COL_KB_AC,
        COL_KB_DT,
        COL_KB_DE,
        COL_KB_GN,
        COL_KB_OS,
        COL_KB_OG,
        COL_KB_OC,
        COL_KB_OX,
        COL_KB_OH,
        COL_KB_REF,
        COL_KB_CC,
        COL_KB_DR,
        COL_KB_PE,
        COL_KB_KW,
        COL_KB_FT,
        COL_KB_SEQ
    };

    static const std::string columnNames[];

private:
    const size_t dbColumns;
    bool isInEntry;
    bool hasEntry;
    std::ostringstream* streams;

};



#endif //MMSEQS_UNIPROTKB_H
