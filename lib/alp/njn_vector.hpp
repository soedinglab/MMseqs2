#ifndef INCLUDED_NJN_VECTOR
#define INCLUDED_NJN_VECTOR

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

File name: njn_vector.hpp

Author: John Spouge

Contents: Matrix routines

******************************************************************************/

#include <vector>
#include "njn_approx.hpp"

namespace Njn {

   template <typename T> 
   class Vector {

   public:

      inline Vector (); 
      inline Vector (const Vector &vector_); 
      inline Vector (const T *begin_, const T *end_); // vector_ [0...m-1]
      inline Vector (size_t m_, const T &a_); // copies a_ into m_ vector elements
      virtual inline ~Vector ();

      inline Vector &operator= (const Vector &vector_);
      virtual inline void *create (bool isCopy_ = false) const; 
      virtual inline void copy (const T *begin_, const T *end_); 
      virtual inline void copy (size_t m_, const T &a_); 

      virtual inline operator bool () const;

      virtual inline std::ostream &out (std::ostream &ostr_) const;
      virtual inline std::istream &in (std::istream &istr_);

      virtual inline Vector &operator= (const T &a_); // copies a_ into all vector elements

      virtual inline size_t size () const; // vector size
      virtual inline T &operator [] (size_t i_); // i-th element
      virtual inline const T &operator [] (size_t i_) const; // i-th element

      virtual inline size_t getM () const {return d_m;} // m
      virtual inline const T *getVector () const {return d_vector_p;} // the vector

   private:

      size_t d_m; // m
      T *d_vector_p; // the vector

      virtual inline void init (size_t m_);
      virtual inline void free2 ();

   };

	}

template <typename T> 
bool operator== (const Njn::Vector <T> &vector_, const Njn::Vector <T> &vector0_);

template <typename T> 
std::ostream &operator<< (std::ostream &ostr_, const Njn::Vector <T> &vector_);

template <typename T> 
std::istream &operator>> (std::istream &istr_, Njn::Vector <T> &vector_);

//
// There are no more declarations beyond this point.
//

namespace Njn {

   template <typename T> 
   Vector <T>::Vector () 
      : d_m (0), d_vector_p (0)
   {}

   template <typename T> 
   Vector <T>::Vector (
   const T *begin_,
   const T *end_)
      : d_m (0), d_vector_p (0)
   {
      copy (begin_, end_);
   }

   template <typename T> 
   Vector <T>::Vector ( // copies a_ into m_ vector elements 
   size_t m_, // m
   const T &a_) // the vector elements
      : d_m (0), d_vector_p (0)
   {
      copy (m_, a_);
   }

   template <typename T> 
   Vector <T>::~Vector () {free2 ();}

   template <typename T> 
   Vector <T> &Vector <T>::operator= (const Vector <T> &vector_)
   {
      if (this != &vector_) copy (vector_.getVector (), vector_.getVector () + vector_.getM ());
      return *this;
   }

   template <typename T> 
   void *Vector <T>::create (bool isCopy_) const 
   {
      Vector *vector = new Vector;
      if (isCopy_) vector->operator= (*this);
      return vector; 
   }

   template <typename T> 
   void Vector <T>::copy ( 
   const T *begin_,
   const T *end_)
   {
      assert (begin_ <= end_);

      size_t m = end_ - begin_;
      if (m != getM ()) {
         free2 ();
         init (m);
      }

      for (size_t i = 0; i < getM (); i++) {
         d_vector_p [i] = begin_ [i];
      }
   }

   template <typename T> 
   void Vector <T>::copy ( 
   size_t m_, // m
   const T &a_) // the vector elements
   {
      if (m_ != getM ()) {
         free2 ();
         init (m_);
      }

      for (size_t i = 0; i < getM (); i++) {
         d_vector_p [i] = a_;
      }
   }

   template <typename T> 
   Vector <T>::operator bool () const {return getM () > 0;}

   template <typename T> 
   Vector <T> &Vector <T>::operator= (const T &a_) // copies a_ into all vector elements
   {
      copy (getM (), a_);
      return *this;
   }


   template <typename T> 
   std::ostream &Vector <T>::out (std::ostream &ostr_) const
   {
      assert (ostr_);
      
	  using namespace Njn;
	  using namespace IoUtil;
	  using namespace std;

      size_t i = 0;

      switch (getFormat ()) {

      case MACHINE :

         for (i = 0; i < getM (); i++) { 
   
            if (i != 0) ostr_ << endl;
            ostr_ << getVector () [i];
         }

         break;

      case HUMAN :
      default :

         ostr_ << getM () << "\t! dimension of vector\n";

         for (i = 0; i < getM (); i++) {

            if (i != 0) ostr_ << '\t';
            ostr_ << getVector () [i];
         }

         ostr_ << "\t! vector elements";

         break;
      } 

      clearFormat ();

      return ostr_;
   }

   template <typename T> 
   std::istream &Vector <T>::in (std::istream &istr_)
   {
      assert (istr_);

      size_t i = 0;
      std::vector <T> v;
      T value;

      //USING_SCOPE(Njn);
      //USING_SCOPE(IoUtil);

	  using namespace Njn;
	  using namespace IoUtil;
	  using namespace std;


      string s;
      stringstream sstream;
      size_t m = 0;

      switch (getFormat ()) {

      case MACHINE :

         while (istr_ >> value) v.push_back (value);

         if (v.size () != getM ()) {

            free2 ();
            init (v.size ());
         }

         std::copy (v.begin (), v.end (), d_vector_p);

         break;

      case HUMAN :
      default :

         getLine (istr_, s);
         sstream.str ("");
         sstream << s;
         sstream >> m;
         if (sstream.fail ()) IoUtil::abort ("Njn::Vector::in : bad m");

         if (m != getM ()) {

            free2 ();
            init (m);
         }

         sstream.str ("");
         getLine (istr_, s);
         sstream << s;
         for (size_t i = 0; i < getM (); i++) {

            sstream >> d_vector_p [i];

            if (sstream.fail ()) {
               ostringstream sstream;
               sstream << i;
               IoUtil::abort ("Njn::Vector::in : bad d_vector_p [" + sstream.str () + "]");
            }
         }

         break;
      }

      clearFormat ();

      return istr_;
   }

   template <typename T> 
   size_t Vector <T>::size () const
   {
      return getM ();
   }

   template <typename T> 
   T &Vector <T>::operator [] (size_t i_)
   {
      return d_vector_p [i_];
   }

   template <typename T> 
   const T &Vector <T>::operator [] (size_t i_) const
   {
      return d_vector_p [i_];
   }

   template <typename T> 
   void Vector <T>::init (size_t m_)
   {
      if (m_ > 0) d_vector_p = new T [m_];

      d_m = m_;
   }

   template <typename T> 
   void Vector <T>::free2 ()
   {
      if (getM () > 0) delete [] d_vector_p; d_vector_p = 0;

      d_m = 0;
   }

}

template <typename T> 
bool operator== (const Njn::Vector <T> &vector_, const Njn::Vector <T> &vector0_)
{
   if (vector_.getM () != vector0_.getM ()) return false;
   for (size_t i = 0; i < vector_.getM (); i++) {
      if (vector_.getVector () [i] != vector0_.getVector () [i]) return false;
   }

   return true;
}

template <typename T> 
std::ostream &operator<< (std::ostream &ostr_, const Njn::Vector <T> &vector_)
{return vector_.out (ostr_);}

template <typename T> 
std::istream &operator>> (std::istream &istr_, Njn::Vector <T> &vector_)
{return vector_.in (istr_);}



#endif //! INCLUDED

