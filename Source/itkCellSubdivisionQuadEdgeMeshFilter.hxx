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

#ifndef __itkCellSubdivisionQuadEdgeMeshFilter_hxx
#define __itkCellSubdivisionQuadEdgeMeshFilter_hxx

#include "itkCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
template< typename TInputMesh, typename TOutputMesh >
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CellSubdivisionQuadEdgeMeshFilter()
{
  m_EdgesPointIdentifier = EdgePointIdentifierContainer::New();
  m_Resolution = 1;
}

template< typename TInputMesh, typename TOutputMesh >
void 
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::GenerateData()
{
  this->CopyInputMeshToOutputMesh();

  while ( m_Resolution )
  {
  OutputMeshPointer result = OutputMeshType::New();
		
  BeforeCellsSubdivision(result);
  this->m_EdgesPointIdentifier->Initialize();
  OutputCellsContainerPointer cells = this->GetOutput()->GetCells();
  for ( OutputCellsContainerIterator cellIt = cells->Begin(); cellIt != cells->End(); ++cellIt )
  {
  CellSubdivision(cellIt->Value(),result);
  }
  AfterCellsSubdivision(result);
		
  --m_Resolution;
  }
}

template< typename TInputMesh, typename TOutputMesh >
void 
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::BeforeCellsSubdivision(OutputMeshType *output)
{
  if ( !output )
  {
  itkExceptionMacro(<< "output Mesh is not defined");
  }
	
  OutputPointsContainerPointer points = this->GetOutput()->GetPoints();
  output->GetPoints()->Reserve( this->GetOutput()->GetNumberOfPoints() );
  for ( OutputPointsContainerIterator ptIt = points->Begin(); ptIt != points->End(); ++ptIt )
  {
  OutputPointType opt;
  opt.CastFrom( ptIt->Value() );
  output->SetPoint(ptIt->Index(), opt);
  }
}

template< typename TInputMesh, typename TOutputMesh >
void 
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::AfterCellsSubdivision(OutputMeshType *output)
{
  this->GetOutput()->Graft(output);
}

template< typename TInputMesh, typename TOutputMesh >
void
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  std::cout << indent << "Subdivision Resolution: " << m_Resolution << std::endl;
}
}
#endif
