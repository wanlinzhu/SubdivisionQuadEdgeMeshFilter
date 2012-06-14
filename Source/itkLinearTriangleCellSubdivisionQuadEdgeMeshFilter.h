/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef __itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter_h
#define __itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter_h

#include "itkCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
/**
 * \class LinearTriangleCellSubdivisionQuadEdgeMeshFilter
 *
 * \brief FIXME     Add documentation here
 * \ingroup ITK-QuadEdgeMeshFiltering
 */
template< typename TInputMesh, typename TOutputMesh >
class LinearTriangleCellSubdivisionQuadEdgeMeshFilter:
  public CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef LinearTriangleCellSubdivisionQuadEdgeMeshFilter                      Self;
  typedef CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >         Superclass;
  typedef SmartPointer< Self >                                                 Pointer;
  typedef SmartPointer< const Self >                                           ConstPointer;

  typedef typename Superclass::InputMeshType    InputMeshType;
  typedef typename Superclass::InputMeshPointer InputMeshPointer;

  typedef typename Superclass::OutputMeshType                OutputMeshType;
  typedef typename Superclass::OutputMeshPointer             OutputMeshPointer;
  typedef typename Superclass::OutputPointsContainerPointer  OutputPointsContainerPointer;
  typedef typename Superclass::OutputPointsContainerIterator OutputPointsContainerIterator;
  typedef typename Superclass::OutputPointType               OutputPointType;
  typedef typename Superclass::OutputVectorType              OutputVectorType;
  typedef typename Superclass::OutputCoordType               OutputCoordType;
  typedef typename Superclass::OutputPointIdentifier         OutputPointIdentifier;
  typedef typename Superclass::OutputCellIdentifier          OutputCellIdentifier;
  typedef typename Superclass::OutputCellType                OutputCellType;
  typedef typename Superclass::OutputQEType                  OutputQEType;
  typedef typename Superclass::OutputMeshTraits              OutputMeshTraits;
  typedef typename Superclass::OutputPointIdIterator         OutputPointIdIterator;

  typedef typename Superclass::EdgePointIdentifierContainer              EdgePointIdentifierContainer;
  typedef typename Superclass::EdgePointIdentifierContainerPointer       EdgePointIdentifierContainerPointer;
  typedef typename Superclass::EdgePointIdentifierContainerConstIterator EdgePointIdentifierContainerConstIterator;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(LinearTriangleCellSubdivisionQuadEdgeMeshFilter,
               CellSubdivisionQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

protected:
  LinearTriangleCellSubdivisionQuadEdgeMeshFilter() {}
  ~LinearTriangleCellSubdivisionQuadEdgeMeshFilter() {}

  virtual void CellSubdivision(OutputCellType *cell, OutputMeshType *output);

private:
  LinearTriangleCellSubdivisionQuadEdgeMeshFilter(const Self &);
  void operator=(const Self &);
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter.hxx"
#endif

#endif
