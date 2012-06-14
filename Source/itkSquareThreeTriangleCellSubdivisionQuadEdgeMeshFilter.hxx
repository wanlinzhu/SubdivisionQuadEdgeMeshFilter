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
#ifndef __itkSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter_hxx
#define __itkSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter_hxx

#include "itkSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
template< typename TInputMesh, typename TOutputMesh >
void
SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CellSubdivision( OutputCellType *cell, OutputMeshType *output )
{
  if ( cell->GetType() != OutputCellType::POLYGON_CELL || cell->GetNumberOfPoints() != 3 )
    {
    itkExceptionMacro(<<" The input cell is not a triangle cell");
    }

  OutputPointIdentifier oldPointIdArray[3];
  OutputPointIdentifier newPointId;

  OutputPointType outPoint;

  outPoint.Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);

  OutputPointIdIterator pter = cell->PointIdsBegin();
  OutputPointIdentifier nn = 0;

  while ( pter != cell->PointIdsEnd() )
    {
    oldPointIdArray[nn] = *pter;
    outPoint += this->GetOutput()->GetPoints()->ElementAt(*pter).GetVectorFromOrigin();
    ++pter;
    ++nn;
    }

  typename OutputPointType::ValueType den = 1. / static_cast< typename OutputPointType::ValueType >( nn );

  for ( unsigned int kk = 0; kk < OutputPointType::PointDimension; ++kk )
    {
    outPoint[kk] *= den;
    }

  newPointId = output->GetNumberOfPoints();
  output->SetPoint(newPointId, outPoint);

  for ( unsigned int ii = 0; ii < 3; ++ii )
    {
    unsigned int jj = ( ii + 1 ) % 3;
    output->AddFaceTriangle(oldPointIdArray[ii], oldPointIdArray[jj], newPointId);

    OutputQEType *edge = output->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);

    if ( edge &&
         !this->m_EdgesPointIdentifier->IndexExists(edge) &&
         !this->m_EdgesPointIdentifier->IndexExists( edge->GetSym() ) )
      {
      this->m_EdgesPointIdentifier->InsertElement(edge, newPointId);
      }
    }
}

template< typename TInputMesh, typename TOutputMesh >
void
SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::SwapEdges( OutputMeshType *output )
{
  if ( !output )
    {
    itkExceptionMacro(<< "output Mesh is not defined");
    }

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0).GetPointer() );
  CopyMeshToMeshPoints( output, this->GetOutput() );

  OutputPointIdentifier pointIdArray[2][2];
  for ( EdgePointIdentifierContainerIterator et = this->m_EdgesPointIdentifier->Begin();
        et != this->m_EdgesPointIdentifier->End();
        ++et )
    {
    pointIdArray[0][0] = et->Index()->GetOrigin();
    pointIdArray[0][1] = et->Index()->GetDestination();

    if ( et->Index()->IsAtBorder() )
      {
      if ( et->Index()->IsLeftSet() )
        {
        pointIdArray[1][0] = et->Index()->GetOnext()->GetDestination();

        this->GetOutput()->AddFaceTriangle(pointIdArray[0][0], pointIdArray[0][1],  pointIdArray[1][0]);
        }
      else if ( et->Index()->IsRightSet() )
        {
        pointIdArray[1][0] = et->Index()->GetOprev()->GetDestination();

        this->GetOutput()->AddFaceTriangle(pointIdArray[0][1], pointIdArray[0][0],  pointIdArray[1][0]);
        }
      }
    else
      {
      pointIdArray[1][0] = et->Index()->GetOnext()->GetDestination();
      pointIdArray[1][1] = et->Index()->GetOprev()->GetDestination();

      this->GetOutput()->AddFaceTriangle(pointIdArray[1][0], pointIdArray[1][1],  pointIdArray[0][1]);
      this->GetOutput()->AddFaceTriangle(pointIdArray[1][1], pointIdArray[1][0],  pointIdArray[0][0]);
      }
    }
}

template< typename TInputMesh, typename TOutputMesh >
void
SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::AdaptiveSwapEdges( OutputMeshType *output )
{
  if ( !output )
    {
    itkExceptionMacro(<< "output Mesh is not defined");
    }

  OutputMeshPointer mesh = this->GetOutput();
  OutputCellsContainerPointer cells = mesh->GetCells();

  this->ProcessObject::SetNthOutput( 0, this->MakeOutput(0).GetPointer() );
  CopyMeshToMeshPoints( output, this->GetOutput() );

  typename Superclass::OutputCellIdentifierListType newCellsToBeSubdivided;

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

    OutputPointIdIterator it = cell->PointIdsBegin();
    unsigned int          n = 0;

    while ( it != cell->PointIdsEnd() )
      {
      oldPointIdArray[n] = *it;
      ++it;
      ++n;
      }

    n = 0;
    unsigned int nn = 0;
    OutputQEType  * splitEdges[3];
    OutputQEType  * newEdges[3];

    for ( unsigned int ii = 0; ii < 3; ++ii )
      {
      unsigned int jj = ( ii + 1 ) % 3;

      OutputQEType *edge = output->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);

      if ( edge )
        {
        newEdges[nn] = edge;
        ++nn;
        }

      if ( edge && this->m_EdgesPointIdentifier->IndexExists(edge) )
        {
        splitEdges[n] = edge;
        ++n;
        }
      }

    OutputPointIdentifier pointIdArray[2][2];
    if ( n == 0 )
      {
      if( nn == 0 || !newEdges[0]->IsLeftSet() )
        {
        //copy the input cell
        this->GetOutput()->AddFaceTriangle(oldPointIdArray[0], oldPointIdArray[1],  oldPointIdArray[2]);
        }
      }
    else
      {
      for(unsigned int ii = 0; ii < n; ++ii)
        {

        pointIdArray[0][0] = splitEdges[ii]->GetOrigin();
        pointIdArray[0][1] = splitEdges[ii]->GetDestination();

        if ( splitEdges[ii]->IsAtBorder() )
          {
          if ( splitEdges[ii]->IsLeftSet() )
            {
            pointIdArray[1][0] = splitEdges[ii]->GetOnext()->GetDestination();

            newCellsToBeSubdivided.push_back(this->GetOutput()->AddFaceTriangle(pointIdArray[0][0], pointIdArray[0][1],  pointIdArray[1][0])->GetLeft());
            }
          else if ( splitEdges[ii]->IsRightSet() )
            {
            pointIdArray[1][0] = splitEdges[ii]->GetOprev()->GetDestination();

            newCellsToBeSubdivided.push_back(this->GetOutput()->AddFaceTriangle(pointIdArray[0][1], pointIdArray[0][0],  pointIdArray[1][0])->GetLeft());
            }
          }
        else
          {
          if ( splitEdges[ii]->IsRightSet() )
            {
            pointIdArray[1][0] = splitEdges[ii]->GetOnext()->GetDestination();
            pointIdArray[1][1] = splitEdges[ii]->GetOprev()->GetDestination();

            newCellsToBeSubdivided.push_back(this->GetOutput()->AddFaceTriangle(pointIdArray[1][0], pointIdArray[1][1],  pointIdArray[0][1])->GetLeft());
            newCellsToBeSubdivided.push_back(this->GetOutput()->AddFaceTriangle(pointIdArray[1][1], pointIdArray[1][0],  pointIdArray[0][0])->GetLeft());
            }
          else
            {
            pointIdArray[1][0] = splitEdges[ii]->GetOnext()->GetDestination();
            newCellsToBeSubdivided.push_back(this->GetOutput()->AddFaceTriangle(pointIdArray[0][0], pointIdArray[0][1],  pointIdArray[1][0])->GetLeft());
            }
          }
        }
      }
    }
  this->m_CellsToBeSubdivided.swap(newCellsToBeSubdivided);
}

template< typename TInputMesh, typename TOutputMesh >
void
SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::AfterCellsSubdivision( OutputMeshType *output )
{
  if ( this->m_Uniform )
    {
    this->SwapEdges(output);
    }
  else
    {
    this->AdaptiveSwapEdges(output);
    }
}
}
#endif
