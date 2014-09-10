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

    Example testing global options.
    Usage:
        pass a list of arguments to the program, without specifying any option.
        For example: ./program one two three
*/

#include <iostream>
#include "getoptpp/getopt_pp.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
    std::vector<std::string> args;

    GetOpt_pp ops(argc, argv);

    ops >> GlobalOption(args);

    std::cout << "RAN: " << ops.app_name() << " ";
    for (std::vector<std::string>::const_iterator it = args.begin(); it != args.end(); ++it)
        std::cout << *it << " ";

    std::cout << std::endl;

    return 0;
}

