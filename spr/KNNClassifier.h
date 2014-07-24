
////////////////////////////////////////////////////////////////////////////////
//
// K-nearest neighbors classifier
//
// Speeded up using prototype condensation based on Gabriel graph neighbors.
// For n prototypes, this operation is O(n^3). Do not use this when number
// of prototypes >>> than number of features to be classified.
//
////////////////////////////////////////////////////////////////////////////////

// prastawa@cs.unc.edu 8/2004

#ifndef _KNNClassifier_h
#define _KNNClassifier_h

#include "DynArray.h"

#include "vnl/vnl_matrix.h"
#include "vnl/vnl_vector.h"

class KNNClassifier
{

public:

  typedef vnl_matrix<float> MatrixType;
  typedef vnl_vector<float> VectorType;

  typedef DynArray<unsigned char> LabelArray;

  KNNClassifier();
  ~KNNClassifier();

  unsigned char Classify(const VectorType& x);
  unsigned char ClassifyWithoutCondense(const VectorType& x);

  void SetTrainingData(MatrixType& prototypes, LabelArray& labels);

  void SetDimension(unsigned int d);
  inline unsigned int GetDimension() const { return m_Dimension; }

  void SetKNeighbors(unsigned int k);
  inline unsigned int GetKNeighbors() const { return m_KNeighbors; }

  void Update();

protected:

  unsigned char _ClassifyFunction(
    const MatrixType& protos, const LabelArray& labels, const VectorType& x);

private:

  MatrixType m_Prototypes;
  MatrixType m_CondensedPrototypes;

  LabelArray m_Labels;
  LabelArray m_CondensedLabels;

  unsigned int m_Dimension;
  unsigned int m_KNeighbors;

  // Range of label values [0, m]
  unsigned int m_LabelRange;

  bool m_Modified;

};

#endif
