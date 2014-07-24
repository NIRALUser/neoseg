
//
// Modified from ITK's gradient descent optimizer
// Added stop criteria for gradient magnitude convergence
//

/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGradientDescentOptimizer.h,v $
  Language:  C++
  Date:      $Date: 2004/11/04 20:40:41 $
  Version:   $Revision: 1.29 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __GradientDescentOptimizer_h
#define __GradientDescentOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"

class GradientDescentOptimizer : 
    public itk::SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef GradientDescentOptimizer          Self;
  typedef itk::SingleValuedNonLinearOptimizer    Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  typedef itk::SmartPointer<const Self>          ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( GradientDescentOptimizer, itk::SingleValuedNonLinearOptimizer );


  /** Codes of stopping conditions */
  typedef enum {
    Converged,
    MaximumNumberOfIterations,
    MetricError
  } StopConditionType;

  /** Methods to configure the cost function. */
  itkGetConstReferenceMacro( Maximize, bool );
  itkSetMacro( Maximize, bool );
  itkBooleanMacro( Maximize );
  bool GetMinimize( ) const
  { return !m_Maximize; }
  void SetMinimize(bool v)
  { this->SetMaximize(!v); }
  void MinimizeOn()
  { this->MaximizeOff(); }
  void MinimizeOff()
  { this->MaximizeOn(); }
  
  /** Advance one step following the gradient direction. */
  virtual void AdvanceOneStep( void );

  /** Start optimization. */
  void    StartOptimization( void );

  /** Resume previously stopped optimization with current parameters
   * \sa StopOptimization. */
  void    ResumeOptimization( void );

  /** Stop optimization.
   * \sa ResumeOptimization */
  void    StopOptimization( void );

  /** Set the learning rate. */
  itkSetMacro( LearningRate, double );

  /** Get the learning rate. */
  itkGetConstReferenceMacro( LearningRate, double);

  /** Set the number of iterations. */
  itkSetMacro( NumberOfIterations, unsigned long );

  /** Get the number of iterations. */
  itkGetConstReferenceMacro( NumberOfIterations, unsigned long );

  /** Get the current iteration number. */
  itkGetConstMacro( CurrentIteration, unsigned int );

  /** Get the current value. */
  itkGetConstReferenceMacro( Value, double );

  /** Get Stop condition. */
  itkGetConstReferenceMacro( StopCondition, StopConditionType );

  /** Get Stop condition. */
  itkGetConstReferenceMacro( Gradient, DerivativeType );

  itkSetMacro(GradientMagnitudeTolerance, double);
  itkGetConstMacro(GradientMagnitudeTolerance, double);

protected:
  GradientDescentOptimizer();
  virtual ~GradientDescentOptimizer() {};
  void PrintSelf(std::ostream& os, itk::Indent indent) const;


  // made protected so subclass can access
  DerivativeType                m_Gradient; 
  bool                          m_Maximize;
  double                        m_LearningRate;

private:
  GradientDescentOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  bool                          m_Stop;
  double                        m_Value;
  StopConditionType             m_StopCondition;
  unsigned long                 m_NumberOfIterations;
  unsigned long                 m_CurrentIteration;

  double m_GradientMagnitudeTolerance;


};


#endif



