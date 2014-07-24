
//TODO

#ifndef _ImageSetManager_h
#define _ImageSetManager_h

#include "itkImage.h"
#include "itkVector.h"

template <class TImage>
class ImageSetManager: public itk::Object
{

public:

  typedef ImageSetManager Self;

  void UseImageSet(VectorImageType* vecimg);
  void UseImageSet(DynArray<ImagePointer>& imglist);

  DynArray<ImagePointer> GetList();
  VectorImagePointer GetVectorImage();

  void AlignImagesToFirstImage();

protected:

};

Filter(DynArray<ImagePointer> imgset)
{
  typedef Vector<ImagePixelType, n> VectorPixelType;
  typedef itk::Image<VectorPixelType, dim> VectorImageType;

  // Allocate vector image


  // Fill values

  // Filter

  // Assign filtered values to images in the set

}

#endif
