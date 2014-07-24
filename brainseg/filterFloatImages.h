
#ifndef _filterFloatImages_h
#define _filterFloatImages_h

#include "itkImage.h"

#include "DynArray.h"

#include <string>

void
filterFloatImages(
  DynArray<itk::Image<float, 3>::Pointer>& images,
  std::string& method,
  unsigned int iters,
  double dt
);

#endif
