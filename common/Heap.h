
////////////////////////////////////////////////////////////////////////////////
//
// Binary heap using arrays
//
// Element type must have: default constructor, copy ctor, operator=, operator<
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 9/2003

#ifndef _Heap_h
#define _Heap_h

#include "DynArray.h"

template <class T>
class Heap
{

public:

  Heap();
  Heap(const Heap& h);
  ~Heap();

  Heap& operator=(const Heap& h);

  void Allocate(unsigned int size);

  inline void Clear() { m_Elements.Clear(); }

  T ExtractMinimum();

  inline unsigned int GetNumberOfElements() { return m_Elements.GetSize(); }

  void Insert(const T& e);

  bool IsEmpty();

  T* GetElements() { return m_Elements.GetRawArray(); }

  void UpdateElementAt(unsigned int i);

private:

  DynArray<T> m_Elements;

  void PreserveHeapOrder();

};

// Get the first k sorted elements using heap sort
template <class T>
T*
heapFirstK(T* array, unsigned int n, unsigned int k);

// Get the k-th element using heap sort
template <class T>
T
heapKthElement(T* array, unsigned int n, unsigned int k);

// Get median using heap sort
template <class T>
T
heapMedian(T* array, unsigned int n);

#ifndef MU_MANUAL_INSTANTIATION
#include "Heap.txx"
#endif

#endif
