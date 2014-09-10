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

#define ARGC    (sizeof(argv)/sizeof(const char*))

TEST(GetOptPPTest, just_long_option)
{
    const char* argv[] = {"test", "--name", "Hugo"};
    GetOpt_pp ops(ARGC, argv);
    string name;
    ops >> Option("name", name);
    ASSERT_EQ("Hugo", name);
}

TEST(GetOptPPTest, just_long_option_default)
{
    const char* argv[] = {"test"};
    GetOpt_pp ops(ARGC, argv);
    string name;
    ops >> Option("name", name, "world");
    ASSERT_EQ("world", name);
}

TEST(GetOptPPTest, both_options_long)
{
    const char* argv[] = {"test", "--name", "Hugo"};
    GetOpt_pp ops(ARGC, argv);
    string name;
    ops >> Option('n', "name", name);
    ASSERT_EQ("Hugo", name);
}

TEST(GetOptPPTest, both_options_short)
{
    const char* argv[] = {"test", "-n", "Hugo"};
    GetOpt_pp ops(ARGC, argv);
    string name;
    ops >> Option('n', "name", name);
    ASSERT_EQ("Hugo", name);
}

TEST(GetOptPPTest, option_not_found)
{
    const char* argv[] = {"test"};
    GetOpt_pp ops(ARGC, argv);
    ops.exceptions(ios::eofbit);
    string name;

    ASSERT_THROW(ops >> Option('n', "name", name), OptionNotFoundEx);
}

TEST(GetOptPPTest, no_manipulators)
{
    const char* argv[] = {"test", "-n", "Hugo"};
    GetOpt_pp ops(ARGC, argv);
    string name;

    ASSERT_EQ("Hugo", ops.getopt<std::string>('n', "name"));
}

TEST(GetOptPPTest, no_manipulators_option_not_found)
{
    const char* argv[] = {"test"};
    GetOpt_pp ops(ARGC, argv);
    ops.exceptions(ios::eofbit);

    ASSERT_THROW(ops.getopt<std::string>('n', "name"), OptionNotFoundEx);
}


TEST(GetOptPPTest, global_options)
{
    const char* argv[] = {"test", "arg1", "arg2"};

    GetOpt_pp ops(ARGC, argv);

    std::vector<std::string> args;
    ops >> GlobalOption(args);

    ASSERT_EQ("test", ops.app_name());

    for (std::vector<std::string>::const_iterator it = args.begin(); it != args.end(); ++it)
    {
        const unsigned int index = it - args.begin() + 1;
        ASSERT_EQ(argv[index], *it);
    }

}

TEST(GetOptPPTest, negative_integer)
{
    const char* argv[] = {"app", "--test", "-1", "-12"};

    GetOpt_pp ops(ARGC, argv);

    int value = 0;
    ops >> Option("test", value);

    ASSERT_EQ(-1, value);
    ASSERT_FALSE(ops >> OptionPresent('1'));
    ASSERT_FALSE(ops >> OptionPresent('2'));
}

TEST(GetOptPPTest, negative_float)
{
    const char* argv[] = {"app", "--test", "-.23"};

    GetOpt_pp ops(ARGC, argv);

    float value = 0;
    ops >> Option("test", value);

    ASSERT_EQ(float(-.23), value);
    ASSERT_FALSE(ops >> OptionPresent('2'));
    ASSERT_FALSE(ops >> OptionPresent('3'));
}

TEST(GetOptPPTest, negative_integer_as_string_token)
{
    const char* argv[] = {"app", "--test", "-1"};

    GetOpt_pp ops(ARGC, argv);

    ASSERT_TRUE(ops >> OptionPresent('1'));
}

TEST(GetOptPPTest, negative_integers_vector)
{
    const char* argv[] = {"app", "--test", "1", "-32", "4"};

    GetOpt_pp ops(ARGC, argv);

    std::vector<int> args;
    ops >> Option("test", args);

    ASSERT_EQ(3, args.size());
}

TEST(GetOptPPTest, negative_integers_as_string_tokens)
{
    const char* argv[] = {"app", "--test", "-3","-32", "4", "-1a"};

    GetOpt_pp ops(ARGC, argv);

    std::vector<string> test_args;
    ops >> Option('t', "test", test_args);
    ASSERT_EQ(3, test_args.size());
    ASSERT_TRUE(ops >> OptionPresent('1'));
    ASSERT_TRUE(ops >> OptionPresent('a'));

    ASSERT_FALSE(ops >> OptionPresent('2'));
    ASSERT_FALSE(ops >> OptionPresent('3'));
    ASSERT_FALSE(ops >> OptionPresent('4'));
}

TEST(GetOptPPTest, suboptions_good)
{
    const char* argv[] = {"app", "-x", "@ref/options_file.opt", "-y"};

    GetOpt_pp ops(ARGC, argv);

    EXPECT_TRUE(ops >> OptionPresent('x'));
    EXPECT_TRUE(ops >> OptionPresent('y'));

    int arg;
    EXPECT_TRUE(ops >> Option("first", arg));
    EXPECT_EQ(1, arg);

    EXPECT_TRUE(ops >> Option("second", arg));
    EXPECT_EQ(2, arg);

    EXPECT_TRUE(ops >> Option("third", arg));
    EXPECT_EQ(3, arg);

    EXPECT_TRUE(ops >> Option("fourth", arg));
    EXPECT_EQ(4, arg);

    EXPECT_TRUE(ops >> Option("fifth", arg));
    EXPECT_EQ(5, arg);

    EXPECT_TRUE(ops >> Option("sixth", arg));
    EXPECT_EQ(6, arg);

}

TEST(GetOptPPTest, suboptions_bad)
{
    const char* argv[] = {"app", "-x", "@inexistent.opt", "-y"};
    bool exception_thrown = false;

    try
    {
        GetOpt_pp ops(ARGC, argv);
    }
    catch (const OptionsFileNotFoundEx& e)
    {
        EXPECT_EQ("inexistent.opt", e.targetFile);
        exception_thrown = true;
    }

    EXPECT_TRUE(exception_thrown);
}

TEST(GetOptPPTest, globals_last)
{
    const char* argv[] = {"app", "one", "two", "-a", "-b", "1", "2", "-c", "3"};

    GetOpt_pp ops(ARGC, argv);

    vector<int> ivec;

    EXPECT_TRUE( ops >> OptionPresent('a') );

    EXPECT_TRUE( ops >> Option('b', ivec ) );
    EXPECT_EQ(2, ivec.size());

    ivec.clear();
    EXPECT_TRUE( ops >> Option('c', ivec ) );
    EXPECT_EQ(1, ivec.size());

    vector<string> svec;
    EXPECT_TRUE( ops >> GlobalOption(svec) );
    EXPECT_EQ(2, svec.size() );
}

