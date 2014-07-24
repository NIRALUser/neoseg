
/*
  Chained affine transform

  Forward evaluation:
  Transform = Translate * Rotate * Scale * Skew * Centering

  Backward evaluation:
  Transform = Centering * Skew * Scale * Rotate * Translate

  ----

  For matrix inverse, the "inverse" of each parameter is computed and then
  the composition order is set to the opposite of the current one.

  For translation and rotation, the inverse is the negative values

  For scale, the inverse is (1 / scale)

  The skew matrix is K =
    [1, k1, k2]
    [0,  1, k3]
    [0,  0,  1]
  with inverse
    [1, -k1, k1*k3-k2]
    [0,   1,      -k3]
    [0,   0,        1]
  so inverse k1 is -k1, inverse k2 is (k1*k3-k2) and inverse k3 is -k3

*/

#ifndef _ChainedAffineTransform3D_h
#define _ChainedAffineTransform3D_h

#include "itkAffineTransform.h"

#include "vnl/vnl_matrix.h"

class ChainedAffineTransform3D:
  public itk::AffineTransform<double, 3>
{

public:

  typedef ChainedAffineTransform3D Self;
  typedef itk::AffineTransform<double, 3> Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  itkTypeMacro(ChainedAffineTransform3D, SuperClass);

  itkNewMacro(Self);

  /** Parameters Type   */
  typedef Superclass::ParametersType         ParametersType;
  typedef Superclass::JacobianType           JacobianType;
  typedef Superclass::ScalarType             ScalarType;
  typedef Superclass::InputPointType         InputPointType;
  typedef Superclass::OutputPointType        OutputPointType;
  typedef Superclass::InputVectorType        InputVectorType;
  typedef Superclass::OutputVectorType       OutputVectorType;
  typedef Superclass::InputVnlVectorType     InputVnlVectorType;
  typedef Superclass::OutputVnlVectorType    OutputVnlVectorType;
  typedef Superclass::InputCovariantVectorType InputCovariantVectorType;
  typedef Superclass::OutputCovariantVectorType OutputCovariantVectorType;
  typedef Superclass::MatrixType             MatrixType;
  typedef Superclass::InverseMatrixType      InverseMatrixType;
  typedef Superclass::CenterType             CenterType;
  typedef Superclass::OffsetType             OffsetType;
  typedef Superclass::TranslationType        TranslationType;

  typedef vnl_matrix<double> VNLMatrixType;

  void SetIdentity();

  const ParametersType& GetParameters() const
  { return m_Parameters; }

  void ComputeOffset();

  void SetParameters(const ParametersType& parameters);

  void SetSourceCenter(double x, double y, double z);
  void SetTargetCenter(double x, double y, double z);

  const CenterType& GetSourceCenter() const { return m_SourceCenter; }
  const CenterType& GetTargetCenter() const { return m_TargetCenter; }

  // Update all parameters at once, to save on updates
  void SetAllParameters(const ParametersType& params,
    double cx, double cy, double cz,
    double qx, double qy, double qz,
    bool forward);

  inline void EvaluateForward()
  { m_ForwardEvaluation = true; this->Update(); }
  inline void EvaluateBackward()
  { m_ForwardEvaluation = false; this->Update(); }

  inline bool IsForwardEvaluation() const { return m_ForwardEvaluation; }

  Pointer GetInverse() const;

private:

  ChainedAffineTransform3D();
  ChainedAffineTransform3D(const Self& other);

  ~ChainedAffineTransform3D();

  const Self& operator=(const Self& other);

  // Update the actual transformation matrix using the parameters
  void Update();

  CenterType m_SourceCenter;
  CenterType m_TargetCenter;

  bool m_ForwardEvaluation;

};

#endif
