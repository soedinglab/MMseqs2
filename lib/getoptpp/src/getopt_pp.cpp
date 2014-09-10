/*
GetOpt_pp: Yet another C++ version of getopt.
    This file is part of GetOpt_pp.

    Copyright (C) Daniel Gutson, FuDePAN 2007-2010
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

#include <fstream>

#if __APPLE__
#include <crt_externs.h>
#define environ (*_NSGetEnviron())
#elif _WIN32
#include <Stdio.h>
#define environ _environ
#else
#include <unistd.h>
#endif

#include "getoptpp/getopt_pp.h"

namespace GetOpt
{

GETOPT_INLINE Token* GetOpt_pp::_add_token(const std::string& value, Token::Type type)
{
    Token* const ret = new Token(value, type);
    if (_first_token == NULL)
        _first_token = ret;
    else
        _last_token->link_to(ret);
    _last_token = ret;
    return ret;
}

GETOPT_INLINE void GetOpt_pp::_init_flags()
{
    std::stringstream ss;
    _flags = ss.flags();
}

GETOPT_INLINE void GetOpt_pp::_parse_sub_file(const std::string& file)
{
    std::ifstream ifile(file.c_str());
    if (!ifile)
        throw OptionsFileNotFoundEx(file);

    std::vector<std::string> args;
    std::string arg;

    while (ifile >> arg)
        args.push_back(arg);

    _parse(args);
}

GETOPT_INLINE void GetOpt_pp::_parse(const std::vector<std::string>& args)
{
    bool any_option_processed = false;
    const size_t argc = args.size();

    size_t start = 0;
    if ( _app_name.empty() )
    {
        _app_name = args[0];
        start = 1;
    }

    // parse arguments by their '-' or '--':
    //   (this will be a state machine soon)
    for (size_t i = start; i < argc; i++)
    {
        const std::string& currentArg = args[i];

        if (currentArg[0] == '-' && currentArg.size() > 1)
        {
            // see what's next, differentiate whether it's short or long:
            if (currentArg[1] == '-')
            {
                if ( currentArg.size() > 2 )
                {
                    // long option
                    _longOps[currentArg.substr(2)].token = _add_token(currentArg.substr(2), Token::LongOption);
                }
                else
                {
                    // it's the -- option alone
                    _longOps[currentArg].token = _add_token(currentArg, Token::GlobalArgument);
                }

                any_option_processed = true;
            }
            else
            {
                // check if it is a negative number: rules
                //  * floating point negative numbers are straight classified as 'arguments'
                //  * integer negative numbers of more than 1 digit length are also 'arguments'
                //  * integer negatives of 1 digit length can be either arguments or short options.
                //  * anything else: short options.
                int anInt;
                float aFloat;
                std::stringstream dummy;
                if ( convert(currentArg, anInt, dummy.flags()) == _Option::OK )
                {
                    if ( currentArg.size() > 2 ) // if it's larger than -d (d=digit), then assume it's a negative number:
                        _add_token(currentArg, any_option_processed ? Token::UnknownYet : Token::GlobalArgument);
                    else // size == 2: it's a 1 digit negative number
                        _shortOps[currentArg[1]].token = _add_token(currentArg, Token::PossibleNegativeArgument);
                }
                else if ( convert(currentArg, aFloat, dummy.flags()) == _Option::OK )
                    _add_token(currentArg, any_option_processed ? Token::UnknownYet : Token::GlobalArgument);
                else
                {
                    // short option
                    // iterate over all of them, keeping the last one in currentData
                    // (so the intermediates will generate 'existent' arguments, as of '-abc')
                    for( size_t j = 1; j < currentArg.size(); j++ )
                        _shortOps[currentArg[j]].token = _add_token(std::string(currentArg, j, 1), Token::ShortOption);
                }

                any_option_processed = true;
            }
        }
        else if ( currentArg[0] == '@' && currentArg.size() > 1 )
        {
            // suboptions file
            _parse_sub_file(currentArg.substr(1));
        }
        else
        {
            _add_token(currentArg, any_option_processed ? Token::UnknownYet : Token::GlobalArgument);
        }
    }

    _last = _Option::OK;    // TODO: IMPROVE!!
}

GETOPT_INLINE void GetOpt_pp::_parse_env()
{
    // this will be optimized in version 3
    std::string var_name;
    std::string var_value;
    size_t var = 0;
    std::string::size_type pos;
    OptionData* data;

    while (environ[var] != NULL)
    {
        var_name = environ[var];
        pos = var_name.find('=');

        if (pos != std::string::npos)
        {
            var_value = var_name.substr(pos + 1);
            var_name = var_name.substr(0, pos);

            if (_longOps.find(var_name) == _longOps.end())
            {
                data = &_longOps[var_name];
                data->token = _add_token(var_name, Token::LongOption);
                data->flags = OptionData::Envir;
                _add_token(var_value, Token::OptionArgument);
            }
        }
        else
            (data = &_longOps[var_name])->flags = OptionData::Envir;

        var++;
    }
}


GETOPT_INLINE void GetOpt_pp::_argc_argv_to_vector(int argc, const char* const* const argv, std::vector<std::string>& args)
{
    for (int i = 0; i < argc; i++)
        args.push_back(argv[i]);
}

GETOPT_INLINE GetOpt_pp::TokensDeleter::~TokensDeleter()
{
    Token* next;
    Token* current(_first);
    while (current != NULL)
    {
        next = current->next;
        delete current;
        current = next;
    }
}

GETOPT_INLINE GetOpt_pp::GetOpt_pp(int argc, const char* const* const argv)
    : _exc(std::ios_base::goodbit), _first_token(NULL), _last_token(NULL), _tokens_deleter(_first_token)
{
    _init_flags();
    std::vector<std::string> args;
    _argc_argv_to_vector(argc, argv, args);
    _parse(args);
}

GETOPT_INLINE GetOpt_pp::GetOpt_pp(int argc, const char* const* const argv, _EnvTag)
    : _first_token(NULL), _last_token(NULL), _tokens_deleter(_first_token)
{
    _init_flags();
    std::vector<std::string> args;
    _argc_argv_to_vector(argc, argv, args);
    _parse(args);
    _parse_env();
}

GETOPT_INLINE GetOpt_pp& GetOpt_pp::operator >> (const _Option& opt) throw(GetOptEx)
{
    if (_last != _Option::ParsingError)
    {
        _last = opt(_shortOps, _longOps, _first_token, _flags);

        switch (_last)
        {
            case _Option::OK:
                break;

            case _Option::OptionNotFound:
                if (_exc & std::ios_base::eofbit)
                    throw OptionNotFoundEx();
                break;

            case _Option::BadType:
                if (_exc & std::ios_base::failbit)
                    throw InvalidFormatEx();
                break;

            case _Option::NoArgs:
                if (_exc & std::ios_base::eofbit)
                    throw ArgumentNotFoundEx();
                break;

            case _Option::TooManyArgs:
                if (_exc & std::ios_base::failbit)
                    throw TooManyArgumentsEx();
                break;

            case _Option::OptionNotFound_NoEx:
                break;  // Ok, it will be read by casting to bool

            case _Option::ParsingError:
                break;  // just to disable warning
        }
    }
    else if (_exc & std::ios_base::failbit)
        throw ParsingErrorEx();

    return *this;
}

GETOPT_INLINE GetOpt_pp& GetOpt_pp::operator >> (std::ios_base & (*iomanip)(std::ios_base&))
{
    std::stringstream ss;
    ss.flags(_flags);
    _flags = (ss << iomanip).flags();
    return *this;
}

GETOPT_INLINE bool GetOpt_pp::options_remain() const
{
    bool remain = false;
    ShortOptions::const_iterator it = _shortOps.begin();
    while (it != _shortOps.end() && !remain)
    {
        remain = (it->second.flags == OptionData::CmdLine_NotExtracted);
        ++it;
    }

    if (!remain)
    {
        LongOptions::const_iterator it = _longOps.begin();
        while (it != _longOps.end() && !remain)
        {
            remain = (it->second.flags == OptionData::CmdLine_NotExtracted);
            ++it;
        }
    }

    if (!remain)
    {
        // check for the global arguments:
        Token* token = _first_token;
        while (!remain && token != NULL)
        {
            remain = (token->type == Token::GlobalArgument || token->type == Token::UnknownYet);
            token = token->next;
        }
    }

    return remain;
}

}
