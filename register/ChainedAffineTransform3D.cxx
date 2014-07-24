
#include "ChainedAffineTransform3D.h"

#include <cmath>

ChainedAffineTransform3D
::ChainedAffineTransform3D():
  Superclass()
{
  m_SourceCenter.Fill(0);
  m_TargetCenter.Fill(0);

  this->SetIdentity();
}

ChainedAffineTransform3D
::ChainedAffineTransform3D(const Self& other)
{
  this->SetMatrix(other.GetMatrix());
  this->SetOffset(other.GetOffset());
  this->SetTranslation(other.GetTranslation());
  this->m_SourceCenter = other.m_SourceCenter;
  this->m_TargetCenter = other.m_TargetCenter;
  this->m_Parameters = other.m_Parameters;
  this->m_ForwardEvaluation = other.m_ForwardEvaluation;
}

ChainedAffineTransform3D
::~ChainedAffineTransform3D()
{

}

const ChainedAffineTransform3D&
ChainedAffineTransform3D
::operator=(const Self& other)
{
  this->SetMatrix(other.GetMatrix());
  this->SetOffset(other.GetOffset());
  this->SetTranslation(other.GetTranslation());
  this->m_SourceCenter = other.m_SourceCenter;
  this->m_TargetCenter = other.m_TargetCenter;
  this->m_Parameters = other.m_Parameters;
  this->m_ForwardEvaluation = other.m_ForwardEvaluation;

  return *this;
}

void
ChainedAffineTransform3D
::ComputeOffset()
{
  // Overload superclass method to do nothing
  // this should be handled in Update()
}

void
ChainedAffineTransform3D
::SetIdentity()
{
  TranslationType trans;
  trans.Fill(0);
  this->SetTranslation(trans);

  MatrixType id;
  id.SetIdentity();
  this->SetMatrix(id);

  OffsetType offset;
  for(unsigned int i = 0; i < 3; i++)
    offset[i] = 0;
  this->SetOffset(offset);

  m_ForwardEvaluation = true;

  m_Parameters.SetSize(12);
  m_Parameters.Fill(0);
  m_Parameters[6] = 1.0;
  m_Parameters[7] = 1.0;
  m_Parameters[8] = 1.0;

  this->Modified();
}

void
ChainedAffineTransform3D
::SetParameters(const ParametersType& parameters)
{
  this->m_Parameters = parameters;

  this->Update();
}

void
ChainedAffineTransform3D
::SetSourceCenter(double x, double y, double z)
{
  m_SourceCenter[0] = x;
  m_SourceCenter[1] = y;
  m_SourceCenter[2] = z;

  this->Update();
}

void
ChainedAffineTransform3D
::SetTargetCenter(double x, double y, double z)
{
  m_TargetCenter[0] = x;
  m_TargetCenter[1] = y;
  m_TargetCenter[2] = z;

  this->Update();
}

void
ChainedAffineTransform3D
::SetAllParameters(const ParametersType& parameters,
  double cx, double cy, double cz,
  double qx, double qy, double qz,
  bool forward)
{
  this->m_Parameters = parameters;

  m_SourceCenter[0] = cx;
  m_SourceCenter[1] = cy;
  m_SourceCenter[2] = cz;

  m_TargetCenter[0] = qx;
  m_TargetCenter[1] = qy;
  m_TargetCenter[2] = qz;

  m_ForwardEvaluation = forward;

  this->Update();
}

void
ChainedAffineTransform3D
::Update()
{
  // Translations
  double tx, ty, tz;
  // Rotations
  double rx, ry, rz;
  // Scalings
  double sx, sy, sz;
  // Skews
  double k1, k2, k3;

  // Extract transform parameters
  tx = m_Parameters[0];
  ty = m_Parameters[1];
  tz = m_Parameters[2];

  rx = m_Parameters[3];
  ry = m_Parameters[4];
  rz = m_Parameters[5];

  sx = m_Parameters[6];
  sy = m_Parameters[7];
  sz = m_Parameters[8];

  k1 = m_Parameters[9];
  k2 = m_Parameters[10];
  k3 = m_Parameters[11];

  double cosx = cos(rx);
  double sinx = sin(rx);

  double cosy = cos(ry);
  double siny = sin(ry);

  double cosz = cos(rz);
  double sinz = sin(rz);

  // Translation matrix
  VNLMatrixType T(4, 4, 0);
  T.set_identity();
  T(0, 3) = tx;
  T(1, 3) = ty;
  T(2, 3) = tz;

  // Rotation matrix forward composition, R = Rx * Ry * Rz
  VNLMatrixType R(4, 4, 0);
  R.set_identity();

  if (m_ForwardEvaluation)
  {
    // Rx * Ry * Rz
    // 1st row
    R(0, 0) = cosy*cosz;
    R(0, 1) = cosy*sinz;
    R(0, 2) = -siny;
    // 2nd row
    R(1, 0) = sinx*siny*cosz - cosx*sinz;
    R(1, 1) = sinx*siny*sinz + cosx*cosz;
    R(1, 2) = sinx*cosy;
    // 3rd row
    R(2, 0) = cosx*siny*cosz + sinx*sinz;
    R(2, 1) = cosx*siny*sinz - sinx*cosz;
    R(2, 2) = cosx*cosy;
  }
  else
  {
    // Rz * Ry * Rx
    // 1st row
    R(0, 0) = cosy*cosz;
    R(0, 1) = cosx*sinz + sinx*siny*cosz;
    R(0, 2) = sinx*sinz - cosx*siny*cosz;
    // 2nd row
    R(1, 0) = -cosy*sinz;
    R(1, 1) = cosx*cosz - sinx*siny*sinz;
    R(1, 2) = sinx*cosz + cosx*siny*sinz;
    // 3rd row
    R(2, 0) = siny;
    R(2, 1) = -sinx*cosy;
    R(2, 2) = cosx*cosy;
  }

  // Scale matrix
  VNLMatrixType S(4, 4, 0);
  S.set_identity();
  S(0, 0) = sx;
  S(1, 1) = sy;
  S(2, 2) = sz;

  // Skew matrix
  VNLMatrixType K(4, 4, 0);
  K.set_identity();
  K(0, 1) = k1;
  K(0, 2) = k2;
  K(1, 2) = k3;

  // Centering matrices
  VNLMatrixType C_s(4, 4, 0);
  C_s.set_identity();

  VNLMatrixType C_t(4, 4, 0);
  C_t.set_identity();

  C_s(0, 3) = -m_SourceCenter[0];
  C_s(1, 3) = -m_SourceCenter[1];
  C_s(2, 3) = -m_SourceCenter[2];

  C_t(0, 3) = m_TargetCenter[0];
  C_t(1, 3) = m_TargetCenter[1];
  C_t(2, 3) = m_TargetCenter[2];

  // Compose matrix
  VNLMatrixType M;
  if (m_ForwardEvaluation)
    M = C_t * T * R * K * S * C_s;
  else
    M = C_s * S * K * R * T * C_t;

  // Update translation with last column in M
  TranslationType trans;
  for (unsigned int i = 0; i < 3; i++)
    trans[i] = M(i, 3);
  this->SetTranslation(trans);

  // Update matrix with the 3x3 component in M
  MatrixType matrix;
  for (unsigned int r = 0; r < 3; r++)
    for (unsigned int c = 0; c < 3; c++)
      matrix[r][c] = M(r, c);
  this->SetMatrix(matrix);

  // Update offset
  OffsetType offset;
  for(unsigned int i = 0; i < 3; i++)
    offset[i] = trans[i];
  this->SetOffset(offset);

  this->Modified();
}

ChainedAffineTransform3D::Pointer
ChainedAffineTransform3D
::GetInverse() const
{

  // Extract scale
  double sx = m_Parameters[6];
  double sy = m_Parameters[7];
  double sz = m_Parameters[8];

  if (fabs(sx) < 1e-20 || fabs(sy) < 1e-20 || fabs(sz) < 1e-20)
    itkExceptionMacro(<< "Transform not invertible");

  // Extract skew
  double k1 = m_Parameters[9];
  double k2 = m_Parameters[10];
  double k3 = m_Parameters[11];

  ParametersType inverseParameters(12);

  // Negative translation and rotation
  for (unsigned int i = 0; i < 6; i++)
    inverseParameters[i] = -m_Parameters[i];

  // Inverse scale
  inverseParameters[6] = 1.0 / sx;
  inverseParameters[7] = 1.0 / sy;
  inverseParameters[8] = 1.0 / sz;

  // Inverse skew
  inverseParameters[9] = -k1;
  inverseParameters[10] = k1*k3 - k2;
  inverseParameters[11] = -k3;

  // Build the inverse transform
  ChainedAffineTransform3D::Pointer inverse = ChainedAffineTransform3D::New();

  inverse->SetAllParameters(
    inverseParameters, 
    -m_SourceCenter[0], -m_SourceCenter[1], -m_SourceCenter[2],
    -m_TargetCenter[0], -m_TargetCenter[1], -m_TargetCenter[2],
    !m_ForwardEvaluation);

  return inverse;
}
