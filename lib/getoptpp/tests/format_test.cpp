/*
    This file is part of GetOpt_pp.

    Copyright (C) Hugo Arregui, FuDePAN 2011
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
*/

#include <string>
#include <vector>
#include <iostream>
#include <mili/mili.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "getoptpp/getopt_pp.h"

using namespace GetOpt;
using namespace std;
using namespace mili;

TEST(GetOptPPFormatTest, hexa)
{
    const char* argv[] = {"test", "-i", "15"}; //15 hex
    GetOpt_pp ops(3, argv);
    int i;
    ops >> hex >> Option('i', "number", i);
    ASSERT_EQ(21, i);
}

struct Date
{
    unsigned int year;
    unsigned int month;
    unsigned int day;

    Date() :
        year(), 
        month(), 
        day() 
    {}

    Date(unsigned int y, unsigned int m, unsigned int d) : 
        year(y), 
        month(m), 
        day(d) 
    {}

    bool valid() const
    {
        return in_range(month, 1U, 12U) && in_range(day, 1U, 31U);
    }
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

TEST(GetOptPPFormatTest, date)
{
    const char* argv[] = {"test", "-d", "12/01/1990"}; // mm/dd/yyyy

    GetOpt_pp ops(3, argv);

    Date date;
    ops >> Option('d', "date", date);

    ASSERT_EQ(12, date.month);
    ASSERT_EQ(1, date.day);
    ASSERT_EQ(1990, date.year);
}
