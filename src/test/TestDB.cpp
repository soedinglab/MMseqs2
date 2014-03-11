#include <iostream>

#include "DBReader.h"
#include "DBWriter.h"

int main(int argc, char **argv)
{
    // DBReader test
    // argv[1] = ffindex_data_file, argv[2] = ffindex_index_file
    DBReader* dbr = new DBReader(argv[1], argv[2]);
    dbr->open(DBReader::NOSORT);

    char* up1 = dbr->getDbKey(0);
    std::cout << "first entry UpID: " << up1 << "\n";

    char* d1 = dbr->getData(0);
    std::cout << "data:\n" << d1;

    int dbsize = dbr->getSize();
    std::cout << "DB size: " << dbsize << "\n";

    char* up2 = dbr->getDbKey(dbsize-1);
    std::cout << "last UpID: " << up2 << "\n";

    char* d2 = dbr->getData(dbsize-1);
    std::cout << "last data:\n" << d2;

    dbr->close();
/*
    // DBWriter test
    DBWriter* dbw = new DBWriter(argv[1], argv[2], 8);
    dbw->open();

    char* key1 = "1";
    char data1[] = "abc\n";
    dbw->write(&data1[0], 4, key1, 0);

    char* key2 = "2";
    char data2[] = "defg\n";
    dbw->write(&data2[0], 5, key2, 1);

    char* key3 = "6";
    char data3[] = "1234567\n";
    dbw->write(&data3[0], 8, key3, 6);

    char* key4 = "4";
    char data4[] = "xyzxxx\n";
    dbw->write(&data4[0], 7, key4, 4);


    dbw->close();*/
}
