
////////////////////////////////////////////////////////////////////////////////
//
// Pair of elements
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 6/2003

#ifndef _Pair_h
#define _Pair_h

template <class T>
class Pair
{

public:

  Pair() { }
  Pair(const Pair& p) { m_First = p.m_First; m_Second = p.m_Second; }
  Pair(const T& a, const T& b) { m_First = a; m_Second = b; }
  ~Pair() { }

  inline Pair& operator=(const Pair& p)
  {
    m_First = p.m_First;
    m_Second = p.m_Second;
    return *this;
  }

  inline T GetFirst() { return m_First; }
  inline T GetSecond() { return m_Second; }

  inline void SetFirst(const T& f) { m_First = f; }
  inline void SetSecond(const T& s) { m_Second = s; }

  inline void Swap()
  {
    T temp;
    temp = m_First;
    m_First = m_Second;
    m_Second = temp;
  }

private:

  T m_First;
  T m_Second;

};

#endif
