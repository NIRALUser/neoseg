
// Modified from ITK's original version

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentOptimizer.cxx,v $
  Language:  C++
  Date:      $Date: 2003/12/18 21:21:42 $
  Version:   $Revision: 1.28 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkExceptionObject.h"

#include "GradientDescentOptimizer.h"

#include "MersenneTwisterRNG.h"

/**
 * Constructor
 */
GradientDescentOptimizer
::GradientDescentOptimizer()
{
  itkDebugMacro("Constructor");

  m_LearningRate = 0.1;
  m_NumberOfIterations = 100;
  m_CurrentIteration = 0;
  m_Maximize = false;
  m_Value = 0.0;
  m_StopCondition = MaximumNumberOfIterations;

  m_GradientMagnitudeTolerance = 1e-5;
}



void
GradientDescentOptimizer
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "LearningRate: "
     << m_LearningRate << std::endl;
  os << indent << "NunberOfIterations: "
     << m_NumberOfIterations << std::endl;
  os << indent << "Maximize: "
     << m_Maximize << std::endl;
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
  os << indent << "Gradient: "
     << m_Gradient;
  os << std::endl;
}


/**
 * Start the optimization
 */
void
GradientDescentOptimizer
::StartOptimization( void )
{

  itkDebugMacro("StartOptimization");
   
  m_CurrentIteration   = 0;

  this->SetCurrentPosition( this->GetInitialPosition() );
  this->ResumeOptimization();

}



/**
 * Resume the optimization
 */
void
GradientDescentOptimizer
::ResumeOptimization( void )
{
  
  itkDebugMacro("ResumeOptimization");

  m_Stop = false;

  MersenneTwisterRNG* rng = MersenneTwisterRNG::GetGlobalInstance();

  InvokeEvent( itk::StartEvent() );
  while( !m_Stop ) 
    {

    try
      {
      m_CostFunction->GetValueAndDerivative( 
        this->GetCurrentPosition(), m_Value, m_Gradient );
      }
    catch( itk::ExceptionObject& err )
      {
      // An exception has occurred. 
      // Terminate immediately.
      m_StopCondition = MetricError;
      StopOptimization();

      // Pass exception to caller
      throw err;
      }

    if( m_Stop )
      {
      break;
      }

    double gradMag = 0;
    for (int i = 0; i < m_Gradient.Size(); i++)
      gradMag += m_Gradient[i]*m_Gradient[i];
    gradMag = sqrt(gradMag);

    if (gradMag < m_GradientMagnitudeTolerance)
    {
      m_StopCondition = Converged;
      StopOptimization();
      break;
    }

    for (int i = 0; i < m_Gradient.Size(); i++)
      m_Gradient[i] += 0.1*gradMag*rng->GenerateNormal(0, 1.0);

    AdvanceOneStep();

    m_CurrentIteration++;

    if( m_CurrentIteration >= m_NumberOfIterations )
      {
      m_StopCondition = MaximumNumberOfIterations;
      StopOptimization();
      break;
      }
    
    }
    

}


/**
 * Stop optimization
 */
void
GradientDescentOptimizer
::StopOptimization( void )
{

  itkDebugMacro("StopOptimization");

  m_Stop = true;
  InvokeEvent( itk::EndEvent() );
}





/**
 * Advance one Step following the gradient direction
 */
void
GradientDescentOptimizer
::AdvanceOneStep( void )
{ 

  itkDebugMacro("AdvanceOneStep");

  double direction;
  if( this->m_Maximize ) 
    {
    direction = 1.0;
    }
  else 
    {
    direction = -1.0;
    }

  const unsigned int spaceDimension = 
    m_CostFunction->GetNumberOfParameters();

  const ParametersType & currentPosition = this->GetCurrentPosition();

  ScalesType scales = this->GetScales();

  // Make sure the scales have been set properly
  if (scales.size() != spaceDimension)
    {
    itkExceptionMacro(<< "The size of Scales is "
                      << scales.size()
                      << ", but the NumberOfParameters for the CostFunction is "
                      << spaceDimension
                      << ".");
    }

  DerivativeType transformedGradient( spaceDimension ); 

  for(unsigned int j = 0; j < spaceDimension; j++)
    {
    transformedGradient[j] = m_Gradient[j] / scales[j];
    }

  ParametersType newPosition( spaceDimension );
  for(unsigned int j = 0; j < spaceDimension; j++)
    {
    newPosition[j] = currentPosition[j] + 
      direction * m_LearningRate * transformedGradient[j];
    }

  this->SetCurrentPosition( newPosition );

  this->InvokeEvent( itk::IterationEvent() );

}

