#ifndef INCLUDED_NJN_MATRIX
#define INCLUDED_NJN_MATRIX

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

File name: njn_matrix.hpp

Author: John Spouge

Contents: Matrix routines

******************************************************************************/


#include "njn_approx.hpp"
#include "njn_doubletype.hpp"
#include "njn_ioutil.hpp"
#include "njn_vector.hpp"


namespace Njn {

    
    template <typename T>
    class Matrix {

    public:

        static inline Matrix <T> *matrix (size_t k_, size_t m_, size_t n_, T a_ = static_cast <T> (0)); 
        // allocates & returns an  a_-initialized pointer to k_ matrices of dimension m_ x n_ 
        // to deallocate : 
        //     delete [] matrix; matrix = 0;

        static inline bool approx (const Matrix <T> &x_, const Matrix <T> &y_, T eps_);
        static inline bool relApprox (const Matrix <T> &x_, const Matrix <T> &y_, T eps_);
        static inline bool absRelApprox (const Matrix <T> &x_, const Matrix <T> &y_, T tol_, T rtol_);
        static inline bool isSymmetric (const Matrix <T> &x_);

        inline Matrix (); 
        inline Matrix (const Matrix &matrix_);
        inline Matrix (size_t m_, size_t n_, const T *vector_); // assumes vector_ contains m_ * n_ elements
        inline Matrix (size_t m_, size_t n_, T a_ = static_cast <T> (0)); // copies a_ into m_ x n_ matrix elements
        virtual inline ~Matrix ();

        inline Matrix &operator= (const Matrix &matrix_);
        virtual inline void *create (bool isCopy_ = false) const; 
        virtual inline void copy (size_t m_, size_t n_, const T *const *matrix_); 
        virtual inline void copy (size_t m_, size_t n_, const T *vector_); // assumes vector_ contains m_ * n_ elements
        virtual inline void copy (size_t m_, size_t n_, T a_ = static_cast <T> (0)); // copies a_ into m_ x n_ matrix elements

        virtual inline operator bool () const;

        virtual inline std::ostream &out (std::ostream &ostr_) const;
        virtual inline std::istream &in (std::istream &istr_);

        virtual inline Matrix &operator= (const T &a_); // copies a_ into all vector elements

        virtual inline T *operator [] (size_t i_); // i-th row
        virtual inline const T *operator [] (size_t i_) const; // i-th row

        virtual inline T &setValue () {return d_value;} // pointer to a typical matrix value for input

        virtual inline size_t getM () const {return d_m;} // m
        virtual inline size_t getN () const {return d_n;} // n
        virtual inline const T *const *getMatrix () const {return d_matrix_p;} // the matrix
        virtual inline const T getValue () const {return d_value;} // pointer to a typical matrix value for input

    private:

        size_t d_m; // m
        size_t d_n; // n
        T **d_matrix_p; // the matrix
        T d_value; // pointer to a typical matrix value for input

        virtual inline void init (size_t m_, size_t n_);
        virtual inline void free2 ();
    };

	namespace MatrixIO {


        // The following routines handle specialized I/O for Matrix.
        // They require I/O for the underlying values to be available.

        enum Format {GENERAL, SYMMETRIC, FORMAT}; 
        // GENERAL    : IoUtil sets the format 
        // SYMMETRIC : IO the upper half of the matrix after checking for symmetry

        Format getFormat (); // gets the format type

        void setFormat (Format format_); // sets the format type
        
        Format clearFormat (); // resets the format type to and returns default

		}
	}

template <typename T> 
bool operator== (const Njn::Matrix <T> &matrix_, const Njn::Matrix <T> &matrix0_);

template <typename S, typename T> 
void copy (Njn::Matrix <S> *matrix_, const Njn::Matrix <T> &matrix0_);

template <typename T> 
std::ostream &operator<< (std::ostream &ostr_, const Njn::Matrix <T> &matrix_);

template <typename T> 
std::istream &operator>> (std::istream &istr_, Njn::Matrix <T> &matrix_);

inline std::ostream &operator << (std::ostream &ostr_, Njn::MatrixIO::Format format_);
inline std::istream &operator >> (std::istream &istr_, Njn::MatrixIO::Format format_);

//
// There are no more declarations beyond this point.
//

namespace Njn {

    template <typename T> 
    Matrix <T> *Matrix <T>::matrix (size_t k_, size_t m_, size_t n_, T a_) 
    // allocates & returns an a_-initialized pointer to k_ matrices of dimension m_ x n_ 
    // to deallocate : 
    //     delete [] matrix; matrix = 0;
    {
        Matrix <T> *matrix = new Matrix <T> [k_];

        for (size_t i = 0; i < k_; i++) 
            matrix [i].copy (m_, n_, a_);

        return matrix;
    }

    template <typename T> 
    bool Matrix <T>::approx (const Matrix <T> &x_, const Matrix <T> &y_, T eps_)
    {
        assert (x_.getM () == y_.getM ());
        assert (x_.getN () == y_.getN ());

        for (size_t i = 0; i < x_.getM (); i++) {
            for (size_t j = 0; j < x_.getN (); j++) {
                if (! Approx::approx (x_ [i][j], y_ [i][j], eps_)) return false;
            }
        }

        return true;
    }

    template <typename T> 
    bool Matrix <T>::relApprox (const Matrix <T> &x_, const Matrix <T> &y_, T eps_)
    {
        assert (x_.getM () == y_.getM ());
        assert (x_.getN () == y_.getN ());

        for (size_t i = 0; i < x_.getM (); i++) {
            for (size_t j = 0; j < x_.getN (); j++) {
                if (! Approx::relApprox (x_ [i][j], y_ [i][j], eps_)) return false;
            }
        }

        return true;
    }

    template <typename T> 
    bool Matrix <T>::absRelApprox (const Matrix <T> &x_, const Matrix <T> &y_, T tol_, T rtol_)
    {
        assert (x_.getM () == y_.getM ());
        assert (x_.getN () == y_.getN ());

        for (size_t i = 0; i < x_.getM (); i++) {
            for (size_t j = 0; j < x_.getN (); j++) {
                if (! Approx::absRelApprox (x_ [i][j], y_ [i][j], tol_, rtol_)) return false;
            }
        }

        return true;
    }

    template <typename T> 
    bool Matrix <T>::isSymmetric (const Matrix <T> &x_)
    {
        if (x_.getM () != x_.getN ()) return false;

        for (size_t i = 0; i < x_.getM (); i++) {
            for (size_t j = 0; j < i; j++) {
                if (x_ [i][j] != x_ [j][i]) return false;
            }
        }

        return true;
    }

    template <typename T> 
    Matrix <T>::Matrix () 
        : d_m (0), d_n (0), d_matrix_p (0), d_value (0)
    {}

    template <typename T> 
    Matrix <T>::Matrix (const Matrix &matrix_) 
        : d_m (0), d_n (0), d_matrix_p (0), d_value (0)
    {
        operator= (matrix_);
    }

    template <typename T> 
    Matrix <T>::Matrix ( 
    size_t m_, // m
    size_t n_, // n
    const T *vector_) // the matrix
        : d_m (0), d_n (0), d_matrix_p (0), d_value (0)
    {
        copy (m_, n_, vector_);
    }

    template <typename T> 
    Matrix <T>::Matrix ( 
    size_t m_, // m
    size_t n_, // n
    T a_) // the matrix elements
        : d_m (0), d_n (0), d_matrix_p (0), d_value (0)
    {
        copy (m_, n_, a_);
    }

    template <typename T> 
    Matrix <T>::~Matrix () {free2 ();}

    template <typename T> 
    Matrix <T> &Matrix <T>::operator= (const Matrix <T> &matrix_)
    {
        if (this != &matrix_) 
        {
            copy (matrix_.getM (), matrix_.getN (), matrix_.getMatrix ());
            this->setValue () = matrix_.getValue ();
        }
        return *this;
    }

    template <typename T> 
    void *Matrix <T>::create (bool isCopy_) const 
    {
        Matrix *matrix = new Matrix;
        if (isCopy_) matrix->operator= (*this);
        return matrix; 
    }

    template <typename T> 
    Matrix <T>::operator bool () const {return getM () > 0 && getN () > 0;}

    template <typename T> 
    std::ostream &Matrix <T>::out (std::ostream &ostr_) const
    {
        assert (ostr_);

	using namespace Njn;
	using namespace Njn::IoUtil;
	using namespace std;

        size_t i = 0;
        size_t j = 0;

        if (MatrixIO::getFormat () == MatrixIO::GENERAL)
        {
            switch (IoUtil::getFormat ()) {

            case IoUtil::MACHINE :

                    for (i = 0; i < getM (); i++) { 
             
						if (i != 0) ostr_ << endl;

                        for (size_t j = 0; j < getN (); j++) {

                        if (j != 0) ostr_ << '\t';
                        ostr_ << getMatrix () [i][j];
                        }
                    }

                    break;

            case IoUtil::HUMAN :
            default :

                    ostr_ << getM () << "\t! row dimension of matrix\n";
                    ostr_ << getN () << "\t! column dimension of matrix\n";

                    for (i = 0; i < getM (); i++) {

                        if (i != 0) ostr_ << "\t! matrix elements\n";

                        for (j = 0; j < getN (); j++) {

                        if (j != 0) ostr_ << '\t';
                        ostr_ << getMatrix () [i][j];
                        }
                    }

                    ostr_ << "\t! matrix elements";

                    break;
            } 
        }
        else if (MatrixIO::getFormat () == MatrixIO::SYMMETRIC)
        {
            if (! isSymmetric (*this)) IoUtil::abort ("Matrix::out : matrix is not symmetric");

            for (i = 0; i < getM (); i++) 
            {
                if (i != 0) ostr_ << endl;

                for (size_t j = i; j < getN (); j++) 
                {
                    if (j != i) ostr_ << '\t';
                    ostr_ << getMatrix () [i][j];
                }
            }
        }
        else IoUtil::abort ("Matrix::out : impossible MatrixIO::getFormat ()");

        IoUtil::clearFormat ();
        MatrixIO::clearFormat ();

        return ostr_;
    }

    template <typename T> 
    std::istream &Matrix <T>::in (std::istream &istr_)
    {
        assert (istr_);

		using namespace Njn;
		using namespace Njn::IoUtil;
		using namespace std;

        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        size_t m = 0;
        size_t n = 0;

        string s;
        stringstream sstream;

        std::vector <T> v;
        T value = this->getValue ();

        if (MatrixIO::getFormat () == MatrixIO::GENERAL)
        {
            switch (IoUtil::getFormat ()) {

            case IoUtil::MACHINE :

                getline (istr_, s);
                sstream.str ("");
                sstream << s;

                while (sstream >> value) v.push_back (value);

                if (! sstream.eof ()) IoUtil::abort ("Njn::Matrix::in : bad value for the MACHINE");

                n = v.size ();
                if (n == 0) IoUtil::abort ("Njn::Matrix::in : bad n for the MACHINE");

                while (istr_ >> value) v.push_back (value);
                if (! sstream.eof ()) IoUtil::abort ("Njn::Matrix::in : bad value for the MACHINE");

                m = v.size () / n;
                if (m * n != v.size ()) IoUtil::abort ("Njn::Matrix::in : rows for the MACHINE have different lengths.");

                if (m != this->getM () || n != this->getN ()) {

                    this->free2 ();
                    this->init (m, n);
                }

                for (i = 0; i < m; i++) {

                    for (j = 0; j < n; j++) d_matrix_p [i][j] = v [k++];
                }

                break;

            case IoUtil::HUMAN :
            default :

                string s;
                stringstream sstream;

                IoUtil::getLine (istr_, s);
                sstream.str ("");
                sstream.clear ();
                sstream << s;
                sstream >> m;
                if (sstream.fail ()) IoUtil::abort ("Njn::Matrix::in : bad m");

                IoUtil::getLine (istr_, s);
                sstream.str ("");
                sstream.clear ();
                sstream << s;
                sstream >> n;
                if (sstream.fail ()) IoUtil::abort ("Njn::Matrix::in : bad n");

                if (m != this->getM () || n != this->getN ()) {

                    this->free2 ();
                    this->init (m, n);
                }

                for (i = 0; i < this->getM (); i++) {

                    sstream.str ("");
                    IoUtil::getLine (istr_, s);
                    sstream.clear ();
                    sstream << s;

                    for (j = 0; j < this->getN (); j++) {

                        sstream >> d_matrix_p [i][j];

                        if (sstream.fail ()) {

                            ostringstream sistream;
                            sistream << i;
                            ostringstream sjstream;
                            sjstream << j;
                            IoUtil::abort ("Njn::Matrix::in : bad d_matrix_p [" + sistream.str () + "][" + sjstream.str () + "]");
                        }
                    }
                }

                break;
            }
        }
        else if (MatrixIO::getFormat () == MatrixIO::SYMMETRIC)
        {
            getline (istr_, s);
            sstream.str (s);
            sstream.clear ();

            while (sstream >> value) 
            {
                v.push_back (value);
            }
            if (! sstream.eof ()) IoUtil::abort ("Njn::Matrix::in : bad value for the MatrixIO::SYMMETRIC");

            m = n = v.size ();
            if (n == 0) IoUtil::abort ("Njn::Matrix::in : bad n for the MatrixIO::SYMMETRIC");

            while (istr_ >> value) v.push_back (value);
            if (! istr_.eof ()) IoUtil::abort ("Njn::Matrix::in : bad value for the MatrixIO::SYMMETRIC");


            if (m * (m + 1) / 2 != v.size ()) IoUtil::abort ("Njn::Matrix::in : rows for the MatrixIO::SYMMETRIC have incorrect lengths.");

            if (m != this->getM () || n != this->getN ()) {

                this->free2 ();
                this->init (m, n);
            }

            for (i = 0; i < m; i++) 
            {
                for (j = 0; j < i; j++) d_matrix_p [i][j] = d_matrix_p [j][i];
                for (j = i; j < n; j++) d_matrix_p [i][j] = v [k++];
            }
        }
        else IoUtil::abort ("Matrix::in : impossible MatrixIO::getFormat ()");

        IoUtil::clearFormat ();
        MatrixIO::clearFormat ();

        return istr_;
    }

    template <typename T> 
    Matrix <T> &Matrix <T>::operator= (const T &a_)
    {
        copy (this->getM (), this->getN (), a_);
        return *this;
    }

    template <typename T> 
    T *Matrix <T>::operator [] (size_t i_)
    {
        return d_matrix_p [i_];
    }

    template <typename T> 
    const T *Matrix <T>::operator [] (size_t i_) const
    {
        return d_matrix_p [i_];
    }

    template <typename T> 
    void Matrix <T>::copy ( 
    size_t m_, // m
    size_t n_, // n
    const T *const *matrix_) // the matrix
    {
        if (m_ != this->getM () || n_ != this->getN ()) {
            this->free2 ();
            this->init (m_, n_);
        }

        for (size_t i = 0; i < this->getM (); i++)
        {
            for (size_t j = 0; j < this->getN (); j++) 
            {
                d_matrix_p [i][j] = matrix_ [i][j];
            }
        }
    }

    template <typename T> 
    void Matrix <T>::copy ( 
    size_t m_, // m
    size_t n_, // n
    const T *vector_) // the matrix
    {
        if (m_ != this->getM () || n_ != this->getN ()) {
            this->free2 ();
            this->init (m_, n_);
        }

        size_t k = 0;
        for (size_t i = 0; i < this->getM (); i++) 
        {
            for (size_t j = 0; j < this->getN (); j++) 
            {
                d_matrix_p [i][j] = vector_ [k++];
            }
        }
    }

    template <typename T> 
    void Matrix <T>::copy ( 
    size_t m_, // m
    size_t n_, // n
    T a_) // the matrix elements
    {
        if (m_ != this->getM () || n_ != this->getN ()) 
        {
            this->free2 ();
            this->init (m_, n_);
        }

        for (size_t i = 0; i < this->getM (); i++)
        {
            for (size_t j = 0; j < this->getN (); j++) 
            {
                d_matrix_p [i][j] = a_;
            }
        }
   }

    template <typename T> 
    void Matrix <T>::init (size_t m_, size_t n_)
    {
        if (m_ > 0) d_matrix_p = new T* [m_];

        for (size_t i = 0; i < m_; i++) 
        {
            d_matrix_p [i] = new T [n_];
        }

        d_m = m_;
        d_n = n_;
    }

    template <typename T> 
    void Matrix <T>::free2 ()
    {
        for (size_t i = 0; i < this->getM (); i++) {
            delete [] d_matrix_p [i]; d_matrix_p [i] = 0;
        }
        if (this->getM () > 0) delete [] d_matrix_p; d_matrix_p = 0;

        d_m = 0;
        d_n = 0;
    }

}

template <typename S, typename T> 
void copy (Njn::Matrix <S> *matrix_, const Njn::Matrix <T> &matrix0_)
{
    matrix_->copy (matrix0_.getM (), matrix0_.getN ());
    
    for (size_t i = 0; i < matrix0_.getM (); i++) {
        for (size_t j = 0; j < matrix0_.getN (); j++) {
            (*matrix_) [i][j] = static_cast <S> (matrix0_ [i][j]);
        }
    }

}

template <typename T> 
bool operator== (const Njn::Matrix <T> &matrix_, const Njn::Matrix <T> &matrix0_)
{
    if (matrix_.getM () != matrix0_.getM ()) return false;
    if (matrix_.getN () != matrix0_.getN ()) return false;
    for (size_t i = 0; i < matrix_.getM (); i++) {
        for (size_t j = 0; j < matrix_.getN (); j++) {
            if (matrix_.getMatrix () [i][j] != matrix0_.getMatrix () [i][j]) return false;
        }
    }

    return true;
}

template <typename T> 
std::ostream &operator<< (std::ostream &ostr_, const Njn::Matrix <T> &matrix_)
{return matrix_.out (ostr_);}

template <typename T> 
std::istream &operator>> (std::istream &istr_, Njn::Matrix <T> &matrix_)
{return matrix_.in (istr_);}

std::ostream &operator << (std::ostream &ostr_, Njn::MatrixIO::Format format_)
{
   Njn::MatrixIO::setFormat (format_);
   return ostr_;
}

std::istream &operator >> (std::istream &istr_, Njn::MatrixIO::Format format_)
{
   Njn::MatrixIO::setFormat (format_);
   return istr_;
}

#endif //! INCLUDED

