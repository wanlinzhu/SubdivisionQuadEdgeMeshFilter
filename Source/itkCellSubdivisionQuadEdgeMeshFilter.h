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

#ifndef __itkCellSubdivisionQuadEdgeMeshFilter_h
#define __itkCellSubdivisionQuadEdgeMeshFilter_h

#include "itkQuadEdgeMeshToQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"
#include "itkTriangleHelper.h"
#include "itkMapContainer.h"

namespace itk
{
/**
 * \class CellSubdivisionQuadEdgeMeshFilter
 *
 * \brief FIXME
 * \ingroup ITK-QuadEdgeMeshFiltering
 */
template< typename TInputMesh, typename TOutputMesh >
class CellSubdivisionQuadEdgeMeshFilter:
  public QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef CellSubdivisionQuadEdgeMeshFilter                           Self;
  typedef QuadEdgeMeshToQuadEdgeMeshFilter< TInputMesh, TOutputMesh > Superclass;
  typedef SmartPointer< Self >                                        Pointer;
  typedef SmartPointer< const Self >                                  ConstPointer;

  typedef TInputMesh                      InputMeshType;
  typedef typename InputMeshType::Pointer InputMeshPointer;

  typedef TOutputMesh                                      OutputMeshType;
  typedef typename OutputMeshType::Pointer                 OutputMeshPointer;
  typedef typename OutputMeshType::PointsContainerPointer  OutputPointsContainerPointer;
  typedef typename OutputMeshType::PointsContainerIterator OutputPointsContainerIterator;
  typedef typename OutputMeshType::CellsContainerPointer   OutputCellsContainerPointer;
  typedef typename OutputMeshType::CellsContainerIterator  OutputCellsContainerIterator;
  typedef typename OutputMeshType::PointType               OutputPointType;
  typedef typename OutputPointType::CoordRepType           OutputCoordType;
  typedef typename OutputMeshType::PointIdentifier         OutputPointIdentifier;
  typedef typename OutputMeshType::CellIdentifier          OutputCellIdentifier;
  typedef typename OutputMeshType::CellType                OutputCellType;
  typedef typename OutputMeshType::QEType                  OutputQEType;
  typedef typename OutputMeshType::MeshTraits              OutputMeshTraits;
  typedef typename OutputMeshType::PointIdIterator         OutputPointIdIterator;

  typedef itk::MapContainer< OutputQEType *, OutputPointIdentifier > EdgePointIdentifierContainer;
  typedef typename EdgePointIdentifierContainer::Pointer             EdgePointIdentifierContainerPointer;
  typedef typename EdgePointIdentifierContainer::Iterator            EdgePointIdentifierContainerIterator;
  typedef typename EdgePointIdentifierContainer::ConstIterator       EdgePointIdentifierContainerConstIterator;

  typedef std::list< OutputCellIdentifier >                     OutputCellIdentifierListType;
  typedef typename OutputCellIdentifierListType::const_iterator OutputCellIdentifierListConstIterator;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(CellSubdivisionQuadEdgeMeshFilter, QuadEdgeMeshToQuadEdgeMeshFilter);

  itkSetMacro(ResolutionLevels, unsigned int);
  itkGetConstMacro(ResolutionLevels, unsigned int);

  itkSetMacro(Uniform, bool);
  itkGetConstMacro(Uniform, bool);
  itkBooleanMacro(Uniform);

  void AddSubdividedCellId(OutputCellIdentifier cellId){m_CellsToBeSubdivided.push_back(cellId);}

protected:
  CellSubdivisionQuadEdgeMeshFilter();

  virtual ~CellSubdivisionQuadEdgeMeshFilter() {}

  virtual void BeforeCellsSubdivision(OutputMeshType *output);

  virtual void CellSubdivision(OutputCellType *cell, OutputMeshType *output) = 0;

  virtual void FixNeighborCells( OutputMeshType* output );

  virtual void AfterCellsSubdivision(OutputMeshType *output);

  virtual void GenerateData();

  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  CellSubdivisionQuadEdgeMeshFilter(const Self &); // purposely not implemented
  void operator=(const Self &);                // purposely not implemented

protected:
  EdgePointIdentifierContainerPointer m_EdgesPointIdentifier;
  OutputCellIdentifierListType        m_CellsToBeSubdivided;
  unsigned int                        m_ResolutionLevels;
  bool                                m_Uniform;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCellSubdivisionQuadEdgeMeshFilter.hxx"
#endif

#endif
