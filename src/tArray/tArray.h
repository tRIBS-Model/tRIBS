/*******************************************************************************
 * TIN-based Real-time Integrated Basin Simulator (tRIBS)
 * Distributed Hydrologic Model
 *
 * Copyright (c) 2025. tRIBS Developers
 *
 * See LICENSE file in the project root for full license information.
 ******************************************************************************/

/***************************************************************************
**
**  tArray.h: Header file for tArray Class objects
**
**  A tArray is an object implementation of a one-dimensional array.
**  A template is used so that the array may be of any specified type.
**  Unlike regular arrays, tArray objects provide bounds checking,
**  memberwise equality/inequality comparison, and memberwise copy
**  operations. The size of the array is determined either by an 
**  argument passed to the constructor or by assignment of one array
**  to another. 
**
***************************************************************************/

#ifndef TARRAY_H
#define TARRAY_H

#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;

//=========================================================================
//
//
//                  Section 1: tArray Class Definition
//
//
//=========================================================================

template< class T > 
class tArray
{
public:
    friend class tFlowNet;
    tArray();                      	// default constructor
    tArray( int number ) :
      npts(number), avalue(0) // constructor initialize array size
    {
      assert( number > 0 );
      avalue = new T [npts];
      assert( avalue != 0 );
      for( int i=0; i<npts; i++ )
         avalue[i] = 0;
    }

    tArray( const tArray< T > & ); 	// copy constructor
    ~tArray()                     	// destructor
    { delete [] avalue;}

    const tArray< T > &operator=( const tArray< T > & ); // memberwise assignmt
    int operator==( const tArray< T > & ) const;    	// memberwise comparison
    int operator!=( const tArray< T > & ) const;     
    T &operator[]( int subscript )      // overloaded array index operator 
    {
      checkSubscript(subscript);
      return avalue[subscript];
    }
  
    const T &operator[]( int ) const;


    int getSize() const        // returns the number of elements in the array
    { return npts; } 

    void setSize( int );       // reinitializes and sets array size
    T *getArrayPtr();          
    void addValueTo( int , T );
    void checkSubscript(int) const;

protected:
    int npts;   		// size of array
    T * avalue; 		// an array
};

//Default constructor

template< class T >
tArray< T >::
tArray() :
        npts(1), avalue(0)
{
    avalue = new T [1];
    assert( avalue != 0 );
    avalue[0] = 0;
}

//Copy constructor

template< class T >
tArray< T >::
tArray( const tArray< T > &original ) :
        npts(original.npts), avalue(0)
{
    int i;
    assert( npts > 0 );
    avalue = new T[npts];
    assert( avalue != 0 );
    for( i = 0; i < npts; i++ )
        avalue[i] = original.avalue[i];
}

//=========================================================================
//
//
//                  Section 2: tArray Overloaded Operators
//
//
//=========================================================================

/**************************************************************************
**  Overloaded operators:
**
**    assignment, equality, inequality - memberwise operation
**    index - uses an assertion to check array bounds (assumed to be within
**           bounds at runtime) and returns value
**    left shift - sends array values to output stream, separated by
**                 spaces and with a carriage return after every 10 items
**
**************************************************************************/

//Overloaded assignment operator

template< class T >
const tArray< T > &tArray< T >::operator=( const tArray< T > &right )
{
    //assert( &right != 0 );//WR--09192023:warning: reference cannot be bound to dereferenced null pointer in well-defined C++ code; comparison may be assumed to always evaluate to true
    int i;
    if( &right != this ){
        delete [] avalue; avalue = 0;
        npts = right.npts;
        if( npts>0 ){
            assert( right.avalue != 0 );
            avalue = new T [npts];
            for( i = 0; i < npts; i++ ){
                avalue[i] = right.avalue[i];
            }
        }
    }
    return *this;
}

//Overloaded equality operator:

template< class T >
int tArray< T >::operator==( const tArray< T > &right ) const
{
    if( npts != right.npts ) return 0;
    int i;
    for( i = 0; i < npts; i++ )
        if( avalue[i] != right.avalue[i] )
            return 0;
    return 1;
}

//overloaded inequality operator:

template< class T >
int tArray< T >::operator!=( const tArray< T > &right ) const
{
    if( npts != right.npts ) return 0;
    int i;
    for( i = 0; i < npts; i++ )
        if( avalue[i] != right.avalue[i] )
            return 1;
    return 0;
}

//checkSubscript function
template< class T >                                               //tArray
void tArray< T >::checkSubscript( int subscript ) const
{
    if( 0 > subscript || subscript >= npts ){
        cout<<"subscript is "<<subscript<<" npts is "<<npts<<endl<<flush;
        cout<<"bailing out of tArray[]"<<endl<<flush;
    }
    assert( 0 <= subscript && subscript < npts );
}

//Overloaded subscript operator 2
template< class T >                                               //tArray
const T &tArray< T >::operator[]( int subscript ) const
{
    checkSubscript(subscript);
    return avalue[subscript];
}

//Overloaded left shift operator

template< class T >
ostream &operator<<( ostream &output, const tArray< T > &a )
{
    int i;

    for( i = 0; i < a.npts; i++ ){
        output << a.avalue[i] << " ";
        if( (i + 1) %10 == 0 ) output<<endl;
    }
    if( i % 10 != 0 ) output<<endl;
    return output;
}


//=========================================================================
//
//
//                  Section 3: tArray Get/Set Functions
//
//
//=========================================================================

//Returns a pointer to the head of the array

template< class T > T *tArray< T >::getArrayPtr() {
    return avalue;
}

//Reinitializes and resized array

template< class T >void tArray<T>::setSize( int size )
{
    int i;
    delete [] avalue;
    npts = size;
    avalue = new T [npts];
    assert( avalue!=0 && npts>=0 );
    for( i=0; i<npts; i++ ) avalue[i] = 0;
}

//Adds value to Array member

template< class T > void tArray< T >::addValueTo( int subscript, T value){
    avalue[subscript] = avalue[subscript] + value;
}

#endif

//=========================================================================
//
//
//                           End of tArray.h
//
//
//=========================================================================
