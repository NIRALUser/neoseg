
////////////////////////////////////////////////////////////////////////////////
//
// Table containing Most Recently Used items
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 9/2005

#ifndef _MRUTable_h
#define _MRUTable_h

template <class T>
class MRUTable
{

public:
  MRUTable();
  MRUTable(unsigned int n);
  MRUTable(const MRUTable& l);
  ~MRUTable();

  MRUTable<T>& operator=(const MRUTable<T>& l);

  void Clear();

  void Insert(const T& e);

  // For data iteration
  inline void GoToBegin() { m_NumberRead = 0; }
  T& ReadNext();

  inline unsigned int GetSize() const { return m_Size; }
  inline unsigned int GetMaxSize() const { return m_MaxSize; }

  inline T* GetOldest() const { return (m_Array+m_OldestIndex); }

  // Retrieve the ith element, indexed starting from the oldest element
  T& GetElement(unsigned int i);

private:

  T* m_Array;

  unsigned int m_MaxSize;
  unsigned int m_Size;

  unsigned int m_InsertionIndex;
  unsigned int m_OldestIndex;

  unsigned int m_NumberRead;

};

#ifndef MU_MANUAL_INSTANTIATION
#include "MRUTable.txx"
#endif

#endif
