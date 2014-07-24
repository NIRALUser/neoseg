
#ifndef _ImageIteratorSet_h
#define _ImageIteratorSet_h

template <class TIterator>
class ImageIteratorSet
{

public:

  typedef TIterator::ImageType ImageType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::IndexType PixelType;

  typedef ImageIteratorSet Self;

  ImageIteratorSet();

  void AddIterator(const TIterator& it);
  void Clear();

  void GoToBegin();
  void GoToEnd();

  bool IsAtEnd();

  DynArray<PixelType> Get();
  void Set(const DynArray<PixelType>& v);

  IndexType GetIndex();

  Self& operator++();
  Self& operator--();

protected:

  DynArray<TIterator> m_Iterators;

};

#endif
