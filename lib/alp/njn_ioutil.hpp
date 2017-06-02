#ifndef INCLUDED_NJN_IOUTIL
#define INCLUDED_NJN_IOUTIL

/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: njn_ioutil.hpp

Author: John Spouge

Contents: 

******************************************************************************/


#include <sstream>
//#include <iostream>
//#include <fstream>
#include <string.h>
#include <cstdlib>

namespace Njn {
	namespace IoUtil {


        char getTerminator (); // gets the terminator character; default '!'
        void setTerminator (char t_); // sets the terminator character
        char clearTerminator (); // resets the terminator character to and returns default

        enum Format {HUMAN, MACHINE, FORMAT}; // format type

        Format getFormat (); // gets the format type
        void setFormat (Format format_); // sets the format type
        Format clearFormat (); // resets the format type to and returns default

        void abort (); // abort 
        void abort (const std::string &s_);
        inline void abort (const char *s_);

        std::istream &getLine ( 
        std::istream &in_, // input filestream
        std::string &str_, // string before '!' from next line
        const char t_ = getTerminator ()); // termination character (could be '\n'; default '!')
        // behaves like std::getline (in_, str_) but throws away the string after first character t_ && skips lines of white space

        template <typename T>
        std::istream &getString ( 
        std::istream &in_, // input filestream
        T &value_, // string before '!' from next line
        std::string &line_, // string before '!' from next line
        const char t_ = getTerminator ()); // termination character (could be '\n'; default '!')
        // extracts the first non-white-space string from the result of getLine (in_, str_, t_) and throws everything else away

        std::istream &getLine ( 
        std::istream &in_, // input filestream
        std::stringstream &sstr_, // string before '!' from next line
        const char t_ = getTerminator ()); // termination character (could be '\n'; default '!')
        // behaves like std::getline (in_, str_) but throws away the string after first character t_ && skips lines of white space

        std::istream &in (
        std::istream &in_,
        double &x_);

		}
	}

inline std::ostream &operator << (std::ostream &ostr_, Njn::IoUtil::Format format_);
inline std::istream &operator >> (std::istream &istr_, Njn::IoUtil::Format format_);

//
// There are no more declarations beyond this point.
//

namespace Njn {
	namespace IoUtil {

        void abort (const char *s_)
        {
            abort (static_cast <const std::string> (s_));
        }

        template <typename T>
        std::istream &getString ( 
        std::istream &in_, // input filestream
        T &value_, // string before '!' from next line
        std::string &line_, // string before '!' from next line
        const char t_) // termination character (could be '\n')
        { // behaves like getline (in_, str_) but throws away the string after first character t_ && skips lines of white space

            std::string line;
            std::stringstream sstream;
            std::string str;

            getLine (in_, line_, t_);
            sstream << line_;
            sstream >> std::skipws >> str;
            sstream.str (str);
            sstream.clear ();
            sstream >> value_;

            return in_; 
        }

		}
	}


std::ostream &operator << (std::ostream &ostr_, Njn::IoUtil::Format format_)
{
    Njn::IoUtil::setFormat (format_);
    return ostr_;
}

std::istream &operator >> (std::istream &istr_, Njn::IoUtil::Format format_)
{
    Njn::IoUtil::setFormat (format_);
    return istr_;
}


#endif //! INCLUDED

