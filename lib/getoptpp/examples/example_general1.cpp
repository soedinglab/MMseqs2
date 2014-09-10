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
    
    This file is a sample usage (and test case).
    TODO:
        * Support for arrays (of any type):
            opt >> vector<T>
        * More validations
        * Pick good exceptions types
        * Support the 'empty option' arguments
        * fill a structure at once
*/

#include <iostream>
#include "getoptpp/getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
    int test1 = 10;
    float test2 = 3.14f;
    std::string test3 = "hello";
    bool flag;
    std::vector<int> vec_int;

    try
    {
        GetOpt_pp ops(argc, argv);

        ops.exceptions(std::ios::failbit | std::ios::eofbit);

        ops
                >> Option('i', "test1", test1)
                >> Option('f', test2)
                >> Option('s', "string", test3)
                >> OptionPresent('x', "flag", flag)
                >> Option('v', "vec", vec_int)
                ;

        if (!ops.options_remain())
        {
            std::cout << test1 << "\n" << test2 << "\n" << test3 << "\n" << flag << "\n";

            for (std::vector<int>::const_iterator it = vec_int.begin(); it != vec_int.end(); ++it)
                std::cout << *it << " ";

            std::cout << std::endl;

            return 0;
        }
        else
        {
            std::cerr << "too many options" << std::endl;
            return 1;
        }
    }
    catch (const GetOptEx& e)
    {
        std::cerr << "Invalid options\n";
    }
}

