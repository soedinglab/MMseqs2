#ifndef INCLUDED_NJN_MEMUTIL
#define INCLUDED_NJN_MEMUTIL

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

File name: njn_memutil.hpp

Author: John Spouge

Contents: 

******************************************************************************/

#include "njn_ioutil.hpp"


namespace Njn {
	namespace MemUtil {

        inline void *memNew (size_t size);

        inline void *memMore (void *p, size_t size);

        inline void memCpy (void *to, const void *from, size_t size);

        inline void *memSet (void *p, int c, size_t size);

        inline void *memZero (void *p, size_t size);

        template <typename T>
        T **newMatrix (size_t m_, size_t n_); // matrix [0...m_)[0...n_)

        template <typename T>
        void deleteMatrix (T **&matrix, size_t m_, size_t n_); 
        // Assumes matrix was created with the call 
        //         matrix = newMatrix <T> (m_, n_);

//
// There are no more declarations beyond this point.
//


        inline void *memNew (size_t size)
        { 
            return (void *) (size == 0 ? 0 : new char [size]);
        }

        inline void *memMore (void *p, size_t size)
        {
            return (void *) (p == 0 ? new char [size] : realloc (p, size));
        }

        inline void memCpy (void *to, const void *from, size_t size)
        { 
	        if (size != 0) memcpy (to, from, size);
        }

        inline void *memSet (void *p, int c, size_t size)
        {
            return memset (p, c, size);
        }

        inline void *memZero (void *p, size_t size)
        {
            return memset (p, '\0', size);
        }

        template <typename T>
        T **newMatrix (size_t m_, size_t n_) // matrix [iBegin...iEnd)[jBegin...jEnd)
        {
            T **matrix = new T * [m_];
            
            for (size_t index = 0; index != m_; index++) 
            {
                matrix [index] = new T [n_];
            }
        
            return matrix;
        }

        template <typename T>

        void deleteMatrix (T **&matrix, size_t m_ , size_t /*sls deleted n_*/) 
        // Assumes matrix was created with the call 
        //         matrix = newMatrix <T> (iBegin, iEnd, jBegin, jEnd);
        {
            for (size_t index = 0; index != m_; index++) 
            {
                delete [] matrix [index];
            }

            delete [] matrix;
            matrix = 0;
        }

		}
	}



#endif

