
////////////////////////////////////////////////////////////////////////////////
//
// Tree structured vector quantization
//
// Build tree with centroids in each node, allows for fast nearest point
// search
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 1/2006

#ifndef _TreeStructuredVectorQuantizer_h
#define _TreeStructuredVectorQuantizer_h

#include "DynArray.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class TreeStructuredVectorQuantizer
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_vector<float> VectorType;

  TreeStructuredVectorQuantizer();
  TreeStructuredVectorQuantizer(MatrixType& samples);
  TreeStructuredVectorQuantizer(DynArray<VectorType>& samples);

  ~TreeStructuredVectorQuantizer();

  void ConstructTree(MatrixType& samples);

  void SetCodeWordLimit(unsigned int c);
  inline unsigned int GetCodeWordLimit() const { return m_CodeWordLimit; }

  VectorType FindNearestMatch(const VectorType& v);

  void Update();

protected:

  void _LloydSplit(MatrixType& centroids, VectorType& labels, VectorType& data);

private:

  MatrixType m_Data;

  DynArray<VectorType> m_TreeArray;

  unsigned int m_CodeWordLimit;

  bool m_Modified;

};

#endif
