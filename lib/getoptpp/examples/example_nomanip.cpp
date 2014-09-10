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

    Example of a Hello World, without using manipulators.
    Usage:
        short option: -n Name
        long option:  --name Name
*/

#include <iostream>
#include "getoptpp/getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
    GetOpt_pp ops(argc, argv);

    std::cout << "Hello " << ops.getopt<std::string>('n', "name", "world") << "!" << std::endl;

    return 0;
}

