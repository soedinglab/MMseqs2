//
// Created by Martin Steinegger on 2019-02-05.
//

#include "FileUtil.h"
#include "Util.h"
#include "Parameters.h"

int mvdb(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 2);
    std::vector<std::string> files = FileUtil::findDatafiles(par.db1.c_str());
    for(size_t i = 0; i < files.size(); i++){
        //Util::file
        std::string extention = files[i].substr(files[i].find_last_of(".") + 1);
        if(Util::isNumber(extention)){
            std::string dst = (par.db2 + "." + extention);
            FileUtil::move(files[i].c_str(), dst.c_str());
        } else{
            FileUtil::move(files[i].c_str(), par.db2.c_str());
        }
    }

    FileUtil::move(par.db1Index.c_str(), par.db2Index.c_str());
    FileUtil::move(par.db1dbtype.c_str(), par.db2dbtype.c_str());

    return EXIT_SUCCESS;
}
