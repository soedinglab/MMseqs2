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

#ifndef GETOPT_PP_H
#define GETOPT_PP_H

#include <string>
#include <vector>   // candidate to be removed
#include <map>
#include <sstream>
#include <list>

/*
    DESIGN GOALS:
        - EASY to use
        - EASY to learn
        - mimc STL's streams
        - EASY to extend
*/

#ifndef GETOPT_INLINE
#   define GETOPT_INLINE
#endif

namespace GetOpt
{

struct Token
{
    enum Type
    {
        ShortOption,
        LongOption,
        GlobalArgument,
        GlobalArgumentUsed, // already read, skip in the next read
        OptionArgument,
        PossibleNegativeArgument,
        UnknownYet  // can be a global option, or an option of the previous one
    };

    Type type;
    std::string value;
    Token* next;

    Token(const std::string& value, Type type = UnknownYet)
        : type(type), value(value), next(NULL)
    {}

    bool is_last() const
    {
        return next == NULL;
    }

    void link_to(Token* new_next)
    {
        next = new_next;
    }

    Token* get_next_option_argument() const
    {
        if (is_last())
            return NULL;
        else
        {
            if (next->type == UnknownYet || next->type == OptionArgument || next->type == PossibleNegativeArgument)
                return next;
            else
                return NULL;
        }
    }
};

struct OptionData
{
    enum _Flags
    {
        CmdLine_NotExtracted,
        CmdLine_Extracted,
        Envir
    };

    _Flags flags;
    Token* token;
    OptionData() : flags(CmdLine_NotExtracted) {}
};

typedef std::map<std::string, OptionData> LongOptions;
typedef std::map<char, OptionData> ShortOptions;


struct _Option
{
    enum Result
    {
        OK,
        ParsingError,
        OptionNotFound,
        BadType,
        NoArgs,
        TooManyArgs,
        OptionNotFound_NoEx
    };

    virtual Result operator()(ShortOptions& short_ops, LongOptions& long_ops, Token* first, std::ios::fmtflags flags) const = 0;
    virtual ~_Option() {}

    static const char NO_SHORT_OPT = 0;
protected:
    static void setTokenAsUsed(Token* token, ShortOptions& short_ops, Token::Type usedAs)
    {
        if (token->type == Token::PossibleNegativeArgument)
            short_ops.erase(token->value[1]);

        token->type = usedAs;
    }
};

template <class T> inline _Option::Result convert(const std::string& s, T& result, std::ios::fmtflags flags)
{
    std::stringstream ss;
    ss.clear();
    ss.flags(flags);
    ss << s;
    ss >> result;
    if (ss.fail() || !ss.eof())
        return _Option::BadType;
    else
        return _Option::OK;
}

template <> inline _Option::Result convert<std::string>(const std::string& s, std::string& result, std::ios::fmtflags /*flags*/)
{
    result = s;
    return _Option::OK;
}


template <class T> class _OptionTBase : public _Option
{
    const char short_opt;
    const std::string long_opt;
protected:
    T& target;
    virtual Result _assign(Token* token, std::ios::fmtflags flags, ShortOptions& short_ops) const = 0;

public:
    _OptionTBase(const _OptionTBase<T>& other)
        : _Option(), short_opt(other.short_opt), long_opt(other.long_opt), target(other.target)
    {}

    _OptionTBase(char short_opt, const std::string& long_opt, T& target)
        : short_opt(short_opt), long_opt(long_opt), target(target)
    {}

    virtual Result operator()(ShortOptions& short_ops, LongOptions& long_ops, Token* /*first*/, std::ios::fmtflags flags) const
    {
        Result ret = OptionNotFound;
        ShortOptions::iterator it;
        if (short_opt == _Option::NO_SHORT_OPT)
            it = short_ops.end();
        else
            it = short_ops.find(short_opt);

        if (it != short_ops.end())
        {
            it->second.flags = OptionData::CmdLine_Extracted;
            ret = _assign(it->second.token, flags, short_ops);
        }
        else if (!long_opt.empty())
        {
            LongOptions::iterator it = long_ops.find(long_opt);
            if (it != long_ops.end())
            {
                it->second.flags = OptionData::CmdLine_Extracted;
                ret = _assign(it->second.token, flags, short_ops);
            }
        }

        return ret;
    }
};


template <class T> class _OptionT : public _OptionTBase<T>
{
protected:
    virtual _Option::Result _assign(Token* token, std::ios::fmtflags flags, ShortOptions& short_ops) const
    {
        Token* const option_token = token->get_next_option_argument();
        if (option_token == NULL)
            return _Option::NoArgs;
        else
        {
            this->setTokenAsUsed(option_token, short_ops, Token::OptionArgument);
            return convert<T>(option_token->value, this->target, flags);
        }
    }
public:
    _OptionT(const _OptionT<T>& other)
        : _OptionTBase<T>(other)
    {}

    _OptionT(char short_opt, const std::string& long_opt, T& target)
        : _OptionTBase<T>(short_opt, long_opt, target)
    {}

};

template <class T> class _OptionT<std::vector<T> > : public _OptionTBase<std::vector<T> >
{
protected:
    virtual _Option::Result _assign(Token* token, std::ios::fmtflags flags, ShortOptions& short_ops) const
    {
        Token* option_token = token->get_next_option_argument();
        if (option_token != NULL)
        {
            _Option::Result result;
            //OptionArgs::const_iterator it = args.begin();
            T temp;

            do
            {
                this->setTokenAsUsed(option_token, short_ops, Token::OptionArgument);
                result = convert<T>(option_token->value, temp, flags);
                if (result == _Option::OK)
                    this->target.push_back(temp);

                option_token = option_token->get_next_option_argument();
            }
            while (option_token != NULL && result == _Option::OK);

            return result;
        }
        else
            return _Option::NoArgs;
    }

public:
    _OptionT(const _OptionT<std::vector<T> >& other)
        : _OptionTBase<std::vector<T> >(other)
    {}

    _OptionT(char short_opt, const std::string& long_opt, std::vector<T>& target)
        : _OptionTBase<std::vector<T> >(short_opt, long_opt, target)
    {}
};


template <class T, class BaseOption>
class _DefValOption : public BaseOption
{
    const T default_value;
public:

    _DefValOption(const _DefValOption<T, BaseOption>& other)
        : BaseOption(other), default_value(other.default_value)
    {}

    _DefValOption(char short_opt, const std::string& long_opt, T& target, const T& default_value)
        : BaseOption(short_opt, long_opt, target), default_value(default_value)
    {}

    virtual _Option::Result operator()(ShortOptions& short_ops, LongOptions& long_ops, Token* first, std::ios::fmtflags flags) const
    {
        _Option::Result ret = BaseOption::operator()(short_ops, long_ops, first, flags);

        if (ret == _Option::OptionNotFound)
        {
            this->target = default_value;
            ret = _Option::OK;
        }

        return ret;
    }
};

template <class T>
class _GlobalOption : public _Option
{
    T& target;
    virtual Result operator()(ShortOptions& short_ops, LongOptions& long_ops, Token* first, std::ios::fmtflags flags) const
    {
        // find first token GlobalArgument or UnknownYet (candidate) or PossibleNegativeArgument (candidate too)
        Token* token(first);
        bool found(false);
        while (token != NULL && !found)
        {
            found = (token->type == Token::GlobalArgument || token->type == Token::UnknownYet || token->type == Token::PossibleNegativeArgument);
            if (!found)
                token = token->next;
        }
        if (found)
        {
            this->setTokenAsUsed(token, short_ops, Token::GlobalArgumentUsed);
            return convert<T>(token->value, target, flags);
        }
        else
            return OptionNotFound;
    }
public:
    _GlobalOption(const _GlobalOption<T>& other)
        : target(other.target)
    {}

    _GlobalOption(T& target)
        : target(target)
    {}
};

template <class T>
class _GlobalOption<std::vector<T> > : public _Option
{
    std::vector<T>& target;
    virtual Result operator()(ShortOptions& short_ops, LongOptions& /*long_ops*/, Token* first, std::ios::fmtflags flags) const
    {
        // find first token GlobalArgument or UnknownYet (candidate) or PossibleNegativeArgument (candidate too)
        Token* token(first);
        bool found_any(false);
        T tmp;
        Result res(OK);

        while (token != NULL && res == OK)
        {
            if (token->type == Token::GlobalArgument || token->type == Token::UnknownYet || token->type == Token::PossibleNegativeArgument)
            {
                res = convert<T>(token->value, tmp, flags);
                if (res == OK)
                {
                    this->setTokenAsUsed(token, short_ops, Token::GlobalArgumentUsed);
                    found_any = true;
                    target.push_back(tmp);
                }
            }
            token = token->next;
        }
        if (res == OK)
        {
            if (found_any)
                return res;
            else
                return OptionNotFound;
        }
        else
            return res;
    }
public:
    _GlobalOption(const _GlobalOption<std::vector<T> >& other)
        : target(other.target)
    {}

    _GlobalOption(std::vector<T>& target)
        : target(target)
    {}
};

template <class T>
inline _OptionT<T> Option(char short_opt, const std::string& long_opt, T& target)
{
    return _OptionT<T>(short_opt, long_opt, target);
}

template <class T>
inline _OptionT<T> Option(char short_opt, T& target)
{
    return _OptionT<T>(short_opt, std::string(), target);
}

// LongOpt only
template <class T>
inline _OptionT<T> Option(const std::string& long_opt, T& target)
{
    return _OptionT<T>(_Option::NO_SHORT_OPT, long_opt, target);
}

// Defaulted version
template <class T>
inline _DefValOption<T, _OptionT<T> >
Option(char short_opt, const std::string& long_opt, T& target, const T& def)
{
    return _DefValOption<T, _OptionT<T> >(short_opt, long_opt, target, def);
}

template <class T>
inline _DefValOption<T, _OptionT<T> > Option(char short_opt, T& target, const T& def)
{
    return _DefValOption<T, _OptionT<T> >(short_opt, std::string(), target, def);
}

//  no short opt.
template <class T>
inline _DefValOption<T, _OptionT<T> >
Option(const std::string& long_opt, T& target, const T& def)
{
    return _DefValOption<T, _OptionT<T> >(_Option::NO_SHORT_OPT, long_opt, target, def);
}

// Defaults for strings:
inline _DefValOption<std::string, _OptionT<std::string> >
Option(char short_opt, const std::string& long_opt, std::string& target, const char* def)
{
    return _DefValOption<std::string, _OptionT<std::string> >(short_opt, long_opt, target, def);
}

inline _OptionT<std::string> Option(char short_opt, std::string& target, const char* def)
{
    return _DefValOption<std::string, _OptionT<std::string> >(short_opt, std::string(), target, def);
}

//  no short opt.
inline _DefValOption<std::string, _OptionT<std::string> >
Option(const std::string& long_opt, std::string& target, const char* def)
{
    return _DefValOption<std::string, _OptionT<std::string> >(_Option::NO_SHORT_OPT, long_opt, target, def);
}

// Global Option:
template <class T>
inline _GlobalOption<T> GlobalOption(T& target)
{
    return _GlobalOption<T>(target);
}

class OptionPresent : public _Option
{
    const char short_opt;
    const std::string long_opt;
    bool* const present;
public:
    // two combinations: with/without target, and with/without long opt.

    // WITH long_opt:
    OptionPresent(char short_opt, const std::string& long_opt, bool& present)
        : short_opt(short_opt), long_opt(long_opt), present(&present)
    {}

    OptionPresent(char short_opt, const std::string& long_opt)
        : short_opt(short_opt), long_opt(long_opt), present(NULL)
    {}

    // WITHOUT long_opt:
    OptionPresent(char short_opt, bool& present)
        : short_opt(short_opt), present(&present)
    {}

    OptionPresent(char short_opt)
        : short_opt(short_opt), present(NULL)
    {}

    // WITHOUT short_opt
    OptionPresent(const std::string& long_opt, bool& present)
        : short_opt(_Option::NO_SHORT_OPT), long_opt(long_opt), present(&present)
    {}

    OptionPresent(const std::string& long_opt)
        : short_opt(_Option::NO_SHORT_OPT), long_opt(long_opt), present(NULL)
    {}
protected:
    virtual Result operator()(ShortOptions& short_ops, LongOptions& long_ops, Token* /*first*/, std::ios::fmtflags /*flags*/) const
    {
        bool found;
        ShortOptions::iterator it = short_ops.find(short_opt);

        found = (it != short_ops.end());
        if (found)
        {
            it->second.flags = OptionData::CmdLine_Extracted;
        }
        else if (!long_opt.empty())
        {
            LongOptions::iterator it = long_ops.find(long_opt);
            found = (it != long_ops.end());
            if (found)
            {
                it->second.flags = OptionData::CmdLine_Extracted;
            }
        }

        if (present != NULL)
        {
            *present = found;
            return OK;
        }
        else
        {
            return found ? OK : OptionNotFound_NoEx;
        }
    }
};

class GetOptEx : public std::exception {};
struct ParsingErrorEx : GetOptEx {};
struct InvalidFormatEx : GetOptEx {};
struct ArgumentNotFoundEx : GetOptEx {};
struct TooManyArgumentsEx : GetOptEx {};
struct OptionNotFoundEx : GetOptEx {};
struct TooManyOptionsEx : GetOptEx {};
struct OptionsFileNotFoundEx : GetOptEx
{
    const std::string targetFile;
    OptionsFileNotFoundEx(const std::string& file) : targetFile(file) {}
    ~OptionsFileNotFoundEx() throw() {}
};

enum _EnvTag
{
    Include_Environment
};

class GetOpt_pp
{
    ShortOptions _shortOps;
    LongOptions _longOps;
    std::ios_base::iostate _exc;
    _Option::Result _last;
    std::ios::fmtflags _flags;
    std::string _app_name;
    Token* _first_token;
    Token* _last_token;

    class TokensDeleter
    {
        Token*& _first;
    public:
        TokensDeleter(Token*& first) : _first(first) {}

        GETOPT_INLINE ~TokensDeleter();
    };

    TokensDeleter _tokens_deleter;

    GETOPT_INLINE Token* _add_token(const std::string& value, Token::Type type);
    GETOPT_INLINE void _init_flags();
    GETOPT_INLINE void _parse(const std::vector<std::string>& args);
    GETOPT_INLINE void _parse_env();
    static GETOPT_INLINE void _argc_argv_to_vector(int argc, const char* const* const argv, std::vector<std::string>& args);
    GETOPT_INLINE void _parse_sub_file(const std::string& file);
public:
    GETOPT_INLINE GetOpt_pp(int argc, const char* const* const argv);
    GETOPT_INLINE GetOpt_pp(int argc, const char* const* const argv, _EnvTag);

    std::ios_base::iostate exceptions() const
    {
        return _exc;
    }
    void exceptions(std::ios_base::iostate except)
    {
        _exc = except;
    }
    void exceptions_all()
    {
        _exc = std::ios_base::failbit | std::ios_base::eofbit;
    }

    operator bool() const
    {
        return _last == _Option::OK;
    }

    GETOPT_INLINE bool options_remain() const;

    void end_of_options() const throw(GetOptEx)
    {
        if (options_remain() && (_exc & std::ios_base::eofbit))
            throw TooManyOptionsEx();
    }

    std::ios::fmtflags flags() const
    {
        return _flags;
    }
    void flags(std::ios::fmtflags flags)
    {
        _flags = flags;
    }

    const std::string& app_name() const
    {
        return _app_name;
    }

    GETOPT_INLINE GetOpt_pp& operator >> (const _Option& opt) throw(GetOptEx);

    GETOPT_INLINE GetOpt_pp& operator >> (std::ios_base& (*iomanip)(std::ios_base&));

    // Alternative to manipulators, for those who don't like them: the 'getopt' method :)
    // Combination 1: with long option:
    template <class T> inline T getopt(char short_opt, const std::string& long_opt) throw(GetOptEx)
    {
        T result;
        operator >> (Option(short_opt, long_opt, result));
        return result;
    }

    template <class T> inline T getopt(char short_opt, const std::string& long_opt, const T& def_value)
    {
        T result;
        operator >> (Option(short_opt, long_opt, result, def_value));
        return result;
    }

    // Combination 2: without long option:
    template <class T> inline T getopt(char short_opt) throw(GetOptEx)
    {
        T result;
        operator >> (Option(short_opt, result));
        return result;
    }

    template <class T> inline T getopt(char short_opt, const T& def_value)
    {
        T result;
        operator >> (Option(short_opt, result, def_value));
        return result;
    }

    struct ItCtorData
    {
        ShortOptions::const_iterator short_iter;
        LongOptions::const_iterator  long_iter;
        GetOpt_pp* getopt_pp;
    };

    template <class Container, class Adapter, class OptionType>
    class _iterator
    {
        typename Container::const_iterator _it;
        GetOpt_pp* _getopt_pp;
    public:
        _iterator(const ItCtorData& ctor_data)
        {
            _it = Adapter::adapt(ctor_data);
            _getopt_pp = ctor_data.getopt_pp;
        }

        _iterator() : _getopt_pp(NULL)
        {}

        _iterator<Container, Adapter, OptionType>& operator = (const _iterator<Container, Adapter, OptionType>& other)
        {
            _it = other._it;
            _getopt_pp = other._getopt_pp;
            return *this;
        }

        bool operator != (const _iterator<Container, Adapter, OptionType>& other) const
        {
            return _it != other._it;
        }

        OptionType option() const
        {
            return _it->first;
        }
        OptionType operator*() const
        {
            return option();
        }

        _iterator<Container, Adapter, OptionType>& operator ++()
        {
            ++_it;
            return *this;
        }

        template <class T>
        GetOpt_pp& operator >> (T& t)
        {
            Adapter::extract(t, *_getopt_pp, option());
            return *_getopt_pp;
        }
    };

    ItCtorData begin()
    {
        ItCtorData ret;
        ret.short_iter = _shortOps.begin();
        ret.long_iter  = _longOps.begin();
        ret.getopt_pp  = this;
        return ret;
    }

    ItCtorData end()
    {
        ItCtorData ret;
        ret.short_iter = _shortOps.end();
        ret.long_iter  = _longOps.end();
        ret.getopt_pp  = this;
        return ret;
    }

    struct ShortAdapter
    {
        static ShortOptions::const_iterator adapt(const ItCtorData& data)
        {
            return data.short_iter;
        }

        template <class T>
        static void extract(T& t, GetOpt_pp& getopt_pp, char option)
        {
            getopt_pp >> Option(option, t);
        }
    };

    struct LongAdapter
    {
        static LongOptions::const_iterator adapt(const ItCtorData& data)
        {
            return data.long_iter;
        }

        template <class T>
        static void extract(T& t, GetOpt_pp& getopt_pp, const std::string& option)
        {
            getopt_pp >> Option('\0', option, t);
        }
    };

    typedef _iterator<ShortOptions, ShortAdapter, char> short_iterator;
    typedef _iterator<LongOptions, LongAdapter, const std::string&> long_iterator;
};

class Environment
{
    // Coming soon!
};

}

#endif
