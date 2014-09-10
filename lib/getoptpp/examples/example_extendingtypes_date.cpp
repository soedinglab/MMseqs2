/*
GetOpt_pp: Yet another C++ version of getopt.
    This file is part of GetOpt_pp.

    Copyright (C) Daniel Gutson, FuDePAN 2007-2008
    Distributed under the Boost Software License, Version 1.0.
    (See accompanying file LICENSE_1_0.txt in the root directory or 
    copy at http://www.boost.org/LICENSE_1_0.txt)
    
    GetOpt_pp IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
    SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
    FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    Example of extending a type.
    Usage:
        short option: -d mm/dd/yyyy
        long option:  --date mm/dd/yyyy
*/

#include <iostream>
#include "getoptpp/getopt_pp.h"

using namespace GetOpt;

struct Date
{
    unsigned int year;
    unsigned int month;
    unsigned int day;

    bool valid() const
    {
        return (month >= 1 && month <= 12 && day >= 1 && day <= 31);
    }
    Date() {}
    Date(unsigned int y, unsigned int m, unsigned int d) : year(y), month(m), day(d) {}


};

namespace GetOpt
{
template <> _Option::Result convert<Date>(const std::string& s, Date& d, std::ios::fmtflags)
{
    _Option::Result ret = _Option::BadType;
    Date tmp;
    char slash;
    std::stringstream ss(s);
    if (ss >> tmp.month)
    {
        ss >> slash;
        if (ss >> tmp.day)
        {
            ss >> slash;
            if (ss >> tmp.year)
            {
                if (tmp.valid())
                {
                    ret = _Option::OK;
                    d = tmp;
                }
            }
        }
    }

    return ret;
}
}

int main(int argc, char* argv[])
{
    Date date;
    const Date myBirthday(1977, 7, 31);

    GetOpt_pp ops(argc, argv);

    if (ops >> Option('d', "date", date, myBirthday))
        std::cout << "Date: " << date.month << "-" << date.day << "-" << date.year << std::endl;
    else
        std::cerr << "Invalid date. Enter mm/dd/yyyy" << std::endl;

    return 0;
}

