
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkExceptionObject.h"

#include "vnl/vnl_math.h"

#include "PowellOptimizer.h"

#define POWELL_GOLD 1.618033988749894848 // (1+sqrt(5))/2
#define POWELL_CONJUGATE_GOLD 0.3819660112501051 // 2 - (1-sqrt(5))/2
#define POWELL_EPS 1e-20

#define LIMIT_BRACKET_ITERS 0

PowellOptimizer
::PowellOptimizer()
{
  itkDebugMacro("Constructor");

  m_Order = OrderType(0);

  m_CurrentDirection = ParametersType(0);

  m_InitialSteps = ParametersType(0);
  m_CurrentSteps = ParametersType(0);

  m_SpaceDimension = 0;

  m_Value = 0.0;

  m_MaximumIterations = 20;
  m_CurrentIteration = 0;

  m_FracTol = 1e-5;

  m_BracketMaxIters = 100000;
  m_BracketMaxStep = 10;

  m_BrentMaxIters = 100;
  m_BrentAbsTol = 1e-3;
  m_BrentFracTol = 1e-2;

  m_UseNewDirections = false;
}

void
PowellOptimizer
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "MaximumIterations: "
     << m_MaximumIterations << std::endl;
  os << indent << "CurrentIteration: "
     << m_CurrentIteration;
  os << indent << "Value: "
     << m_Value;
  if (m_CostFunction)
  {
    os << indent << "CostFunction: "
       << m_CostFunction;
  }
  os << indent << "StopCondition: "
     << m_StopCondition;
  os << std::endl;

}

void
PowellOptimizer
::StartOptimization( void )
{

  itkDebugMacro("StartOptimization");
   
  m_CurrentIteration = 0;

  m_SpaceDimension = m_CostFunction->GetNumberOfParameters();

  // Initialize order if necessary
  if (m_Order.GetNumberOfElements() < m_SpaceDimension)
  {
    m_Order = OrderType(m_SpaceDimension);
    for (unsigned i = 0; i < m_SpaceDimension; i++)
      m_Order[i] = i;
  }

  // Build default list of directions
  m_DirectionArray = MatrixType(m_SpaceDimension, m_SpaceDimension);
  m_DirectionArray.fill(0.0);
  for (unsigned int i = 0; i < m_SpaceDimension; i++)
  {
    m_DirectionArray[i][m_Order[i]] = 1.0;
  }

  // Use unit steps if not specified
  if (m_InitialSteps.GetNumberOfElements() != m_SpaceDimension)
  {
    m_InitialSteps = ParametersType(m_SpaceDimension);
    m_InitialSteps.Fill(1.0);
  }

  m_CurrentSteps = m_InitialSteps;

  if (m_CurrentDirection.GetSize() != m_SpaceDimension)
  {
    m_CurrentDirection = ParametersType(m_SpaceDimension);
    m_CurrentDirection.Fill(0.0);
    m_CurrentDirection[m_Order[0]] = 1.0;
  }

  this->SetCurrentPosition( this->GetInitialPosition() );

  this->ResumeOptimization();

}

void
PowellOptimizer
::ResumeOptimization( void )
{
  
  itkDebugMacro("ResumeOptimization");

  m_Stop = false;

  this->InvokeEvent( itk::StartEvent() );

  try
  {
    m_Value = m_CostFunction->GetValue(this->GetCurrentPosition());
  }
  catch( itk::ExceptionObject& err )
  {
     m_StopCondition = EvaluationError;
     StopOptimization();
     // Pass exception to caller
     throw err;
  }

  while( !m_Stop ) 
  {

    try
    {
  
      AdvanceOneStep();

    }
    catch( itk::ExceptionObject& err )
    {
       m_StopCondition = EvaluationError;
       StopOptimization();
       // Pass exception to caller
       throw err;
    }


    if( m_Stop )
    {
      break;
    }

    m_CurrentIteration++;

    if( m_CurrentIteration >= m_MaximumIterations )
    {
       m_StopCondition = MaximumNumberOfIterations;
       StopOptimization();
       break;
    }
    
  }

}

void
PowellOptimizer
::StopOptimization( void )
{

  itkDebugMacro("StopOptimization");

  m_Stop = true;
  this->InvokeEvent( itk::EndEvent() );
}

void
PowellOptimizer
::AdvanceOneStep( void )
{ 

  itkDebugMacro("AdvanceOneStep");

  ParametersType startPosition = this->GetCurrentPosition();

  double startValue = m_Value;

  unsigned int ibig = 0;
  double del = 0;

  // Line search along each parameter
  for (unsigned int dir = 0; dir < m_SpaceDimension; dir++)
  {
    // Get direction
    for (unsigned int j = 0; j < m_SpaceDimension; j++)
      m_CurrentDirection[j] = m_DirectionArray[dir][j];

    double prev = m_Value;

    // Minimizes along the current direction
    // Updates value, current position, and direction
    this->LineSearch(m_CurrentSteps[dir]);

    double diff = vnl_math_abs(m_Value - prev);
    if (diff > del)
    {
      del = diff;
      ibig = dir;
    }
  }

  m_CurrentSteps /= 2.0;

  if ( (2.0*vnl_math_abs(m_Value - startValue))
       <=
       (m_FracTol*(vnl_math_abs(m_Value)+vnl_math_abs(startValue))) )
  {
    this->InvokeEvent( itk::IterationEvent() );
    StopOptimization();
    return;
  }

  ParametersType p = this->GetCurrentPosition();

  ParametersType tempPosition(m_SpaceDimension);
  for(unsigned int j = 0; j < m_SpaceDimension; j++)
    tempPosition[j] = 2.0*p[j] - startPosition[j];

  double tempValue = m_CostFunction->GetValue(tempPosition);

  if (tempValue < startValue)
  {
    double w = (startValue - m_Value - del);
    double v = (startValue - tempValue);
    double t = 2.0*(startValue-2.0*m_Value+tempValue)*(w*w) - del*(v*v);

    if (t < 0.0)
    {
      for(unsigned int j = 0; j < m_SpaceDimension; j++)
        m_CurrentDirection[j] = p[j] - startPosition[j];

      this->LineSearch(1.0);

      // Update list of directions?
      if (m_UseNewDirections)
      {
        // Replace last direction entry with current direction
        for (unsigned int j = 0; j < m_SpaceDimension; j++)
        {
          m_DirectionArray[ibig][j] = m_DirectionArray[m_SpaceDimension-1][j];
          m_DirectionArray[m_SpaceDimension-1][j] = m_CurrentDirection[j];
        }

        // Also update associated step sizes
        m_CurrentSteps[ibig] = m_CurrentSteps[m_SpaceDimension-1];
        m_CurrentSteps[m_SpaceDimension-1] = 1.0;
      }
    }
  }

  this->InvokeEvent( itk::IterationEvent() );

}

void 
PowellOptimizer
::LineSearch(double step)
{
  itkDebugMacro(<< "LineSearch");

  // Bracket
  itkDebugMacro(<< "Bracket");

  double ax, bx, cx;
  double fa, fb, fc;

  double u;
  double fu;

  ax = 0.0;
  bx = step;

  fa = this->EvaluateLineAt(ax);
  fb = this->EvaluateLineAt(bx);

  if (fb > fa)
  {
    double t;
    t = ax; ax = bx; bx = t;
    t = fa; fa = fb; fb = fa;
  }

  cx = bx + POWELL_GOLD*(bx-ax);

  fc = this->EvaluateLineAt(cx);

  unsigned int bracketIters = 0;

  while (fb > fc)
  {

#if LIMIT_BRACKET_ITERS
    if (bracketIters >= m_BracketMaxIters)
    {
      // Use initial bracket values
      ax = 0.0;
      bx = step;
      cx = bx + POWELL_GOLD*(bx-ax); 

      fa = this->EvaluateLineAt(ax);
      fb = this->EvaluateLineAt(bx);
      fc = this->EvaluateLineAt(cx);

      break;
    }
#endif

    bracketIters++;

    double r = (bx-ax)*(fb-fc);
    double q = (bx-cx)*(fb-fa);

    double dq = q-r;
    if (vnl_math_abs(dq) < POWELL_EPS)
      dq = vnl_math_sgn(dq) * POWELL_EPS;

    double ulim = bx + m_BracketMaxStep*(cx-bx);

    u = bx - ((bx-cx)*q-(bx-ax)*r) / (2.0*dq);

    if (((bx-u)*(u-cx)) > 0.0)
    {
      fu = this->EvaluateLineAt(u);
      if (fu < fc)
      {
        ax = bx; bx = u;
        fa = fb; fb = fu;
        break;
      }
      else if (fu > fb)
      {
        cx = u;
        fc = fu;
        break;
      }
      u = cx + POWELL_GOLD*(cx-bx);
      fu = this->EvaluateLineAt(u);
    }
    else if (((cx-u)*(u-ulim)) > 0.0)
    {
      fu = this->EvaluateLineAt(u);
      if (fu < fc)
      {
        bx = cx; cx = u;
        u = cx+POWELL_GOLD*(cx-bx);
        fb = fc; fc = fu; 
        fu = this->EvaluateLineAt(u);
      }
    }
    else if (((u-ulim)*(ulim-cx)) >= 0.0)
    {
      u = ulim;
      fu = this->EvaluateLineAt(u);
    }
    else
    {
      u = cx+POWELL_GOLD*(cx-bx);
      fu = this->EvaluateLineAt(u);
    }

    ax = bx; bx = cx; cx = u;
    fa = fb; fb = fc; fc = fu;

  }

  // Brent
  itkDebugMacro(<< "Brent");

  double a = (ax < cx) ? ax : cx;
  double b = (ax > cx) ? ax : cx;
  double d = 0.0;
  double e = 0.0;

  double x, w, v;
  double fx, fw, fv;

  x = bx; w = bx; v = bx;
  fx = fb; fw = fb; fv = fb;


  for (unsigned iter = 0; iter < m_BrentMaxIters; iter++)
  {

    double xm = 0.5*(a+b);
    double tol1 = m_BrentFracTol*vnl_math_abs(x) + POWELL_EPS;
    double tol2 = 2.0*tol1;
    double dt = 0.5*(b-a);
    double t = vnl_math_abs(x-xm) + dt;
    if (t <= tol2 || dt < m_BrentAbsTol)
      break;

    if (vnl_math_abs(e) > tol1)
    {
       double p, q, r;
       r = (x-w)*(fx-fv);
       q = (x-v)*(fx-fw);
       p = (x-v)*q - (x-w)*r;
       q = 2.0*(q-r);
       if (q > 0.0)
         p = -p;
       q = vnl_math_abs(q);
       double etemp = e;
       e = d;
       if (vnl_math_abs(p) >= vnl_math_abs(0.5*q*etemp)
           || p <= (q*(a-x)) || p >= (q*(b-x)))
       {
         if (x >= xm)
           e = a - x;
         else
           e = b - x;
         d = POWELL_CONJUGATE_GOLD*e;
       }
       else
       {
         d = p/q;
         u = x + d;
         if ((u-a) < tol2 || (b-u) < tol2)
         {
           d = tol1 * vnl_math_sgn(xm-x);
         }
       }
    }
    else
    {
      if (x >= xm)
        e = a-x;
      else
        e = b-x;
      d = POWELL_CONJUGATE_GOLD*e;
    }
    if (vnl_math_abs(d) >= tol1)
      u = x + d;
    else
    {
      u = x + tol1*vnl_math_sgn(d);
    }
    fu = this->EvaluateLineAt(u);
    if (fu <= fx)
    {
      if (u >= x)
        a = x;
      else
        b = x;
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    }
    else
    {
      if (u < x)
        a = u;
      else
        b = u;
      if (fu <= fw || w == x)
      {
        v = w; w = u;
        fv = fw; fw = fu;
      }
      else if (fu <= fv || v == x || v == w)
      {
        v = u;
        fv = fu;
      }
    }

  }

  // Done with line search, update variables
  // fx is the min value along the line
  // x is the min location along the line

  m_Value = fx;

  ParametersType pos = this->GetCurrentPosition();

  ParametersType newPosition(m_SpaceDimension);
  for (unsigned int j = 0; j < m_SpaceDimension; j++)
    newPosition[j] = pos[j] + x*m_CurrentDirection[j];

  this->SetCurrentPosition(newPosition);

}

double
PowellOptimizer
::EvaluateLineAt(double step)
{
  ParametersType orig = this->GetCurrentPosition();

  PowellOptimizer::ParametersType p(m_SpaceDimension);
  for (unsigned int j = 0; j < m_SpaceDimension; j++)
    p[j] = orig[j] + step*m_CurrentDirection[j];

  return m_CostFunction->GetValue(p);
}
