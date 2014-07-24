
////////////////////////////////////////////////////////////////////////////////
//
// ITK wrapped Powell optimizer from Numerical Recipes in C
//
////////////////////////////////////////////////////////////////////////////////

#ifndef _PowellOptimizer_h
#define _PowellOptimizer_h

#include "itkArray.h"
#include "itkMacro.h"
#include "itkSingleValuedNonLinearOptimizer.h"

#include "vnl/vnl_matrix.h"

class PowellOptimizer : public itk::SingleValuedNonLinearOptimizer
{

public:

  // Standard class typedefs
  typedef PowellOptimizer          Self;
  typedef itk::SingleValuedNonLinearOptimizer    Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  typedef itk::SmartPointer<const Self>          ConstPointer;
  
  // Method for creation through the object factory
  itkNewMacro(Self);

  // Run-time type information (and related methods)
  itkTypeMacro(PowellOptimizer, itk::SingleValuedNonLinearOptimizer);

  // Codes of stopping conditions
  typedef enum {
    MaximumNumberOfIterations,
    EvaluationError
  } StopConditionType;

  // Order of optimization for each parameter
  typedef itk::Array<unsigned int> OrderType;

  typedef vnl_matrix<double> MatrixType;

  // Start optimization
  void StartOptimization( void );

  // Resume previously stopped optimization with current parameters
  void ResumeOptimization( void );

  // Stop optimization */
  void StopOptimization( void );

  // Get / set the number of iterations
  itkGetMacro( MaximumIterations, unsigned long );
  itkSetMacro( MaximumIterations, unsigned long );

  // Get the current iteration number. */
  itkGetConstMacro( CurrentIteration, unsigned int );

  // Get the current value
  itkGetConstMacro(Value, double);

  // Get the value given the current position and direction along with a step
  double EvaluateLineAt(double step);

  // Get stop condition
  itkGetConstMacro( StopCondition, StopConditionType );

  itkGetConstMacro(Order, OrderType);
  itkSetMacro(Order, OrderType);

  itkGetConstMacro(InitialSteps, ParametersType);
  itkSetMacro(InitialSteps, ParametersType);

  itkGetConstMacro(FracTol, double);
  itkSetMacro(FracTol, double);

  itkGetConstMacro(BracketMaxIters, unsigned int);
  itkSetMacro(BracketMaxIters, unsigned int);

  itkGetConstMacro(BracketMaxStep, double);
  itkSetMacro(BracketMaxStep, double);

  itkGetConstMacro(BrentMaxIters, unsigned int);
  itkSetMacro(BrentMaxIters, unsigned int);

  itkGetConstMacro(BrentAbsTol, double);
  itkSetMacro(BrentAbsTol, double);

  itkGetConstMacro(BrentFracTol, double);
  itkSetMacro(BrentFracTol, double);

  itkGetConstMacro(UseNewDirections, bool);
  itkSetMacro(UseNewDirections, bool);
  itkBooleanMacro(UseNewDirections);

protected:

  PowellOptimizer();
  virtual ~PowellOptimizer() {};

  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  void AdvanceOneStep();
  void LineSearch(double step);

  unsigned int m_SpaceDimension;

  ParametersType m_CurrentDirection;

  ParametersType m_InitialSteps;
  ParametersType m_CurrentSteps;

  OrderType m_Order;

  bool m_Stop;

  double m_Value;

  StopConditionType m_StopCondition;
  unsigned long m_MaximumIterations;
  unsigned long m_CurrentIteration;

  double m_FracTol;

  unsigned int m_BracketMaxIters;
  double m_BracketMaxStep;

  unsigned int m_BrentMaxIters;
  double m_BrentAbsTol;
  double m_BrentFracTol;

  MatrixType m_DirectionArray;

  bool m_UseNewDirections;

private:

  PowellOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

#endif
