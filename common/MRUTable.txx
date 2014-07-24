
#ifndef _MRUTable_txx
#define _MRUTable_txx

#include "MRUTable.h"

#include "muException.h"

template <class T>
T*
_MRUTable_safeAlloc(unsigned int n)
{ 
  T* array = new T[n];
  if (array == NULL)
    muExceptionMacro(<< "MRUTable: Failed allocating " << n << " elements");
  return array;
}

template <class T>
MRUTable<T>
::MRUTable() : MRUTable(100)
{

}

template <class T>
MRUTable<T>
::MRUTable(unsigned int n)
{

  m_Size = 0;
  m_MaxSize = n;

  //m_Array = new T[m_MaxSize];
  m_Array = _MRUTable_safeAlloc<T>(m_MaxSize);

  m_InsertionIndex = 0;
  m_OldestIndex = m_MaxSize + 1;

  m_NumberRead = 0;

}

template <class T>
MRUTable<T>
::MRUTable(const MRUTable<T>& l)
{
  m_Size = l.m_Size;
  m_MaxSize = l.m_MaxSize;

  //m_Array = new T[m_MaxSize];
  m_Array = _MRUTable_safeAlloc<T>(m_MaxSize);
  for (unsigned int i = 0; i < m_Size; i++)
    m_Array[i] = l.m_Array[i];

  m_InsertionIndex = l.m_InsertionIndex;
  m_OldestIndex = l.m_OldestIndex;

  m_NumberRead = l.m_NumberRead;
}

template <class T>
MRUTable<T>
::~MRUTable()
{
  delete [] m_Array;
}

template <class T>
MRUTable<T>&
MRUTable<T>
::operator=(const MRUTable& l)
{
  delete [] m_Array;

  m_MaxSize = l.m_MaxSize;
  m_Size = l.m_Size;

  //m_Array = new T[m_MaxSize];
  m_Array = _MRUTable_safeAlloc<T>(m_MaxSize);
  for (unsigned int i = 0; i < m_Size; i++)
    m_Array[i] = l.m_Array[i];

  m_InsertionIndex = l.m_InsertionIndex;
  m_OldestIndex = l.m_OldestIndex;

  m_NumberRead = l.m_NumberRead;

  return *this;
}

template <class T>
void
MRUTable<T>
::Clear()
{

  m_Size = 0;

  m_InsertionIndex = 0;
  m_OldestIndex = m_MaxSize+1;

  m_NumberRead = 0;

}

template <class T>
void
MRUTable<T>
::Insert(const T& e)
{
  m_Array[m_InsertionIndex] = e;

  if (m_InsertionIndex == m_OldestIndex) // overwrite oldest value
    m_OldestIndex = (m_OldestIndex + 1) % m_MaxSize;

  if (m_OldestIndex > m_MaxSize) // empty table
    m_OldestIndex = m_InsertionIndex;

  m_InsertionIndex = (m_InsertionIndex + 1) % m_MaxSize;

  if (m_Size < m_MaxSize)
    m_Size++;
}

template <class T>
T&
MRUTable<T>
::ReadNext()
{
  if (m_NumberRead >= m_Size)
    muExceptionMacro(<< "All values in MRU table have been read");

  m_NumberRead++;

  return m_Array[(m_OldestIndex + m_NumberRead - 1) % m_MaxSize];
}

template <class T>
T&
MRUTable<T>
::GetElement(unsigned int i)
{
  if (i >= m_Size)
    muExceptionMacro(<< "Indexing element beyond MRU table");

  return m_Array[(m_OldestIndex + i) % m_MaxSize];
}

#endif
