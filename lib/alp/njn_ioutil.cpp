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

File name: njn_ioutil.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include <math.h>

#include <assert.h>
#include "njn_ioutil.hpp"
#include <iostream>
#include <fstream>


using namespace std;
using namespace Njn;
using namespace IoUtil;



   const Format FORMAT = HUMAN;
   Format format1 = HUMAN;

   const char TERMINATOR = '!';
   char terminator = TERMINATOR;


Format IoUtil::getFormat ()
{
   return format1;
}

void IoUtil::setFormat (Format format_)
{
   format1 = format_;
}

Format IoUtil::clearFormat ()
{
   return format1 = FORMAT;
}

char IoUtil::getTerminator ()
{
   return terminator;
}

void IoUtil::setTerminator (char t_)
{
   terminator = t_;
}

char IoUtil::clearTerminator ()
{
   return terminator = TERMINATOR;
}

void IoUtil::abort ()
{
   exit (1);
}

void IoUtil::abort (const string &s_)
{
   cerr << s_ << endl;
   IoUtil::abort ();
}

istream &IoUtil::getLine ( 
istream &in_, // input filestream
string &str_, // string before '!' from next line
const char t_) // termination character (could be '\n')
{ // behaves like getline (in_, str_) but throws away the string after first character t_ && skips lines of white space

   assert (t_ != '\0');

   if (! in_) return in_;

   const char *pbuf = 0;

   while (getline (in_, str_)) {
      for (pbuf = str_.c_str (); *pbuf != '\0' && isspace (*pbuf); pbuf++); // advance past white space
      if (*pbuf != '\0' && *pbuf != t_) break; // skip lines full of white space
   }

   if (t_ != '\n') {
      size_t pos = str_.find (t_, 0);
      if (pos < str_.size ()) str_.erase (pos, str_.size ());
   }

   return in_; // buf is empty and eof is reached 
}

istream &IoUtil::getLine ( 
istream &in_, // input filestream
stringstream &sstr_, // sstr_ contains the string before '!' from next line
const char t_) // termination character (could be '\n')
{ // behaves like getline (in_, str_) but throws away the string after first character t_ && skips lines of white space

   string s;

   IoUtil::getLine (in_, s, t_);
   sstr_.clear ();
   sstr_.str ("");
   sstr_ << s;
   sstr_.clear ();

   return in_;
}

std::istream &IoUtil::in (
std::istream &in_,
double &x_)
{
    string s;
    in_ >> s;

    for (string::iterator i = s.begin (); i != s.end (); i++) *i = /*sls added (char)*/(char) tolower (*i);

    if (s == "1.#inf")
    {
        x_ = HUGE_VAL;   
    }
    else if (s == "nan")
    {
        x_ = HUGE_VAL;   
    }
    else
    {
        stringstream str (s);
        str >> x_;

        if (str.fail ()) in_.setstate (ios_base::failbit);
    }

    return in_;
}

