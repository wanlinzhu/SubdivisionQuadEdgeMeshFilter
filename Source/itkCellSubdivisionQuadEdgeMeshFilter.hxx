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
  m_ResolutionLevels = 1;
  m_Uniform = true;
}

template< typename TInputMesh, typename TOutputMesh >
void
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::GenerateData()
{
  this->CopyInputMeshToOutputMeshGeometry();

  while ( m_ResolutionLevels != 0 )
    {
    OutputMeshPointer result = OutputMeshType::New();

    BeforeCellsSubdivision(result);

    this->m_EdgesPointIdentifier->Initialize();

    OutputCellsContainerPointer cells = this->GetOutput()->GetCells();

    if( this->m_Uniform )
      {
      for ( OutputCellsContainerIterator cellIt = cells->Begin(); cellIt != cells->End(); ++cellIt )
        {
        CellSubdivision(cellIt->Value(),result);
        }
      }
    else
      {
      OutputCellIdentifierListConstIterator it  = this->m_CellsToBeSubdivided.begin();
      OutputCellIdentifierListConstIterator end = this->m_CellsToBeSubdivided.end();
      while( it != end )
        {
        OutputCellType* cell = cells->GetElement( *it );
        this->CellSubdivision( cell, result );
        ++it;
        }
      }

    AfterCellsSubdivision(result);

    --m_ResolutionLevels;
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
::FixNeighborCells(OutputMeshType *output)
{

  OutputCellIdentifierListType newCellsToBeSubdivided;

  for ( OutputCellsContainerIterator cellIt = output->GetCells()->Begin();
        cellIt != output->GetCells()->End();
        ++cellIt)
    {
    newCellsToBeSubdivided.push_back(cellIt->Index());
    }

  OutputCellsContainerPointer cells = this->GetOutput()->GetCells();

  for ( OutputCellsContainerIterator cellIt = cells->Begin();
        cellIt != cells->End();
        ++cellIt )
    {
    OutputCellType* cell = cellIt->Value();

    if ( cell->GetType() != OutputCellType::POLYGON_CELL || cell->GetNumberOfPoints() != 3 )
      {
      continue;
      }

    OutputPointIdentifier oldPointIdArray[3];
    OutputPointIdentifier newPointIdArray[3];

    OutputPointIdIterator it = cell->PointIdsBegin();
    unsigned int          n = 0;

    while ( it != cell->PointIdsEnd() )
      {
      oldPointIdArray[n] = *it;
      ++it;
      ++n;
      }

    unsigned int splitEdges[3];
    n = 0;
    OutputQEType *edge;
    bool copyCellFlag = true;
    for ( unsigned int ii = 0; ii < 3; ++ii )
      {
      unsigned int jj = ( ii + 1 ) % 3;

     // OutputQEType *edge = this->GetOutput()->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);
      edge = this->GetOutput()->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);

      if ( this->m_EdgesPointIdentifier->IndexExists(edge) )
        {
        copyCellFlag = false;
        break;
        }
      else if(this->m_EdgesPointIdentifier->IndexExists(edge->GetSym()) )
        {
        newPointIdArray[ii] = this->m_EdgesPointIdentifier->GetElement(edge->GetSym());
        splitEdges[n] = ii ;
        ++n;
        }
      }

    if( copyCellFlag && n == 0 )
      {
      // this face has no subdivided face as neighbor, copy it
      output->AddFaceTriangle(oldPointIdArray[0], oldPointIdArray[1], oldPointIdArray[2]);
      }
    else if( n == 1 )
      {
      unsigned int ii = splitEdges[0];
      unsigned int jj = ( ii + 1 ) % 3;
      unsigned int kk = ( ii + 2 ) % 3;

      output->AddFaceTriangle( newPointIdArray[ii], oldPointIdArray[jj], oldPointIdArray[kk] );
      output->AddFaceTriangle( newPointIdArray[ii], oldPointIdArray[kk], oldPointIdArray[ii] );
      }
    else if( n == 2 )
      {
      unsigned int ii = splitEdges[0];
      unsigned int jj = splitEdges[1];

      if( ii == 0 && jj == 1 )
        {
        // ii = 0, jj = 1
        output->AddFaceTriangle( oldPointIdArray[2], oldPointIdArray[0], newPointIdArray[0] );
        output->AddFaceTriangle( oldPointIdArray[2], newPointIdArray[0], newPointIdArray[1] );
        output->AddFaceTriangle( newPointIdArray[0], oldPointIdArray[1], newPointIdArray[1] );
        }
      else if( ii == 0 && jj == 2 )
        {
        // ii = 0, jj = 2
        output->AddFaceTriangle( oldPointIdArray[1], oldPointIdArray[2], newPointIdArray[0] );
        output->AddFaceTriangle( oldPointIdArray[2], newPointIdArray[2], newPointIdArray[0] );
        output->AddFaceTriangle( newPointIdArray[2], oldPointIdArray[0], newPointIdArray[0] );
        }
      else if( ii == 1 && jj == 2 )
        {
        // ii = 1, jj = 2
        output->AddFaceTriangle( oldPointIdArray[0], oldPointIdArray[1], newPointIdArray[1] );
        output->AddFaceTriangle( oldPointIdArray[0], newPointIdArray[1], newPointIdArray[2] );
        output->AddFaceTriangle( newPointIdArray[1], oldPointIdArray[2], newPointIdArray[2] );
        }
      }
    else if( n == 3 )
      {
      // this face was not supposed to be subdivided but all neighbors are
      output->AddFaceTriangle(oldPointIdArray[0], newPointIdArray[0], newPointIdArray[2]);
      output->AddFaceTriangle(newPointIdArray[0], oldPointIdArray[1], newPointIdArray[1]);
      output->AddFaceTriangle(newPointIdArray[1], oldPointIdArray[2], newPointIdArray[2]);
      output->AddFaceTriangle(newPointIdArray[0], newPointIdArray[1], newPointIdArray[2]);
      }
    }

  this->m_CellsToBeSubdivided.swap(newCellsToBeSubdivided);
}

template< typename TInputMesh, typename TOutputMesh >
void
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::AfterCellsSubdivision(OutputMeshType *output)
{
  if( !this->m_Uniform )
    {
    this->FixNeighborCells( output );
    }

  this->GetOutput()->Graft(output);
}

template< typename TInputMesh, typename TOutputMesh >
void
CellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
  std::cout << indent << "Subdivision Resolution Levels: " << m_ResolutionLevels << std::endl;
}
}
#endif
