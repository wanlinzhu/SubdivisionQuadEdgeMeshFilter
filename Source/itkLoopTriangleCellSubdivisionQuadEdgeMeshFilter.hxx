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
#ifndef __itkLoopTriangleCellSubdivisionQuadEdgeMeshFilter_hxx
#define __itkLoopTriangleCellSubdivisionQuadEdgeMeshFilter_hxx

#include "itkLoopTriangleCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
template< typename TInputMesh, typename TOutputMesh >
void
LoopTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CellSubdivision( OutputCellType *cell, OutputMeshType *output)
{
  if ( cell->GetType() != OutputCellType::POLYGON_CELL || cell->GetNumberOfPoints() != 3 )
    {
    itkExceptionMacro(<<" The input cell is not a triangle cell");
    }

  OutputPointIdentifier oldPointIdArray[3];
  OutputPointIdentifier newPointIdArray[3];

  OutputCoordType pointWeight[4] = {0.375, 0.375, 0.125, 0.125};
  OutputPointType pointArray[4];

  OutputPointIdIterator it = cell->PointIdsBegin();
  OutputPointIdentifier numberOfPoints = output->GetNumberOfPoints();
  OutputPointIdentifier n = 0;

  while ( it != cell->PointIdsEnd() )
    {
    oldPointIdArray[n++] = *it;
    ++it;
    }

  for ( unsigned int ii = 0; ii < 3; ++ii )
    {
    unsigned int jj = ( ii + 1 ) % 3;

    OutputQEType *edge = this->GetOutput()->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);

    if ( this->m_EdgesPointIdentifier->IndexExists(edge) )
      {
      newPointIdArray[ii] = this->m_EdgesPointIdentifier->GetElement(edge);
      }
    else if ( this->m_EdgesPointIdentifier->IndexExists(edge->GetSym()) )
      {
      newPointIdArray[ii] = this->m_EdgesPointIdentifier->GetElement(edge->GetSym());
      this->m_EdgesPointIdentifier->InsertElement(edge, newPointIdArray[ii]);
      }
    else if ( edge->IsInternal() )
      {
      OutputPointType outPoint;
      outPoint.Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);

      this->GetOutput()->GetPoint(oldPointIdArray[ii], &pointArray[0]);
      this->GetOutput()->GetPoint(oldPointIdArray[jj], &pointArray[1]);

      if ( edge->GetLnext() )
        {
        this->GetOutput()->GetPoint(edge->GetLnext()->GetDestination(), &pointArray[2]);
        }
      else
        {
        pointArray[2].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      if ( edge->GetRprev() )
        {
        this->GetOutput()->GetPoint(edge->GetRprev()->GetDestination(), &pointArray[3]);
        }
      else
        {
        pointArray[3].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      for ( unsigned int kk = 0; kk < 3; kk++ )
        {
        for ( unsigned int mm = 0; mm < 4; mm++ )
          {
          outPoint[kk] += pointWeight[mm] * pointArray[mm][kk];
          }
        }

      newPointIdArray[ii] = numberOfPoints++;
      this->m_EdgesPointIdentifier->InsertElement(edge, newPointIdArray[ii]);
      output->SetPoint(newPointIdArray[ii], outPoint);
      }
    else if ( edge->IsAtBorder() )
      {
      this->GetOutput()->GetPoint(oldPointIdArray[ii], &pointArray[0]);
      this->GetOutput()->GetPoint(oldPointIdArray[jj], &pointArray[1]);

      OutputPointType outPoint;
      outPoint.SetToMidPoint(pointArray[0], pointArray[1]);
      newPointIdArray[ii] = numberOfPoints++;
      this->m_EdgesPointIdentifier->InsertElement(edge, newPointIdArray[ii]);
      output->SetPoint(newPointIdArray[ii], outPoint);
      }
    else
      {
      itkExceptionMacro(<< "Wire edge detected");
      }
    }
  output->AddFaceTriangle(oldPointIdArray[0], newPointIdArray[0], newPointIdArray[2]);
  output->AddFaceTriangle(newPointIdArray[0], oldPointIdArray[1], newPointIdArray[1]);
  output->AddFaceTriangle(newPointIdArray[1], oldPointIdArray[2], newPointIdArray[2]);
  output->AddFaceTriangle(newPointIdArray[0], newPointIdArray[1], newPointIdArray[2]);
}

template< typename TInputMesh, typename TOutputMesh >
void
LoopTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::BeforeCellsSubdivision( OutputMeshType *output )
{
  if ( !output )
    {
    itkExceptionMacro(<< "output Mesh is not defined");
    }

  OutputPointsContainerPointer points = this->GetOutput()->GetPoints();
  output->GetPoints()->Reserve( this->GetOutput()->GetNumberOfPoints() );
  for ( OutputPointsContainerIterator ptIt = points->Begin(); ptIt != points->End(); ++ptIt )
    {
    OutputPointType ipt = ptIt->Value();
    OutputPointType opt;
    opt.Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
    unsigned int nn = 0;

    OutputPointType bpt;
    bpt.Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
    unsigned int nb = 0;

    OutputQEType *edge = this->GetOutput()->FindEdge( ptIt->Index() );
    typename OutputQEType::IteratorGeom q_it = edge->BeginGeomOnext();
    while ( q_it != edge->EndGeomOnext() )
      {
      if ( q_it.Value()->IsAtBorder() )
        {
        bpt += points->ElementAt( q_it.Value()->GetDestination() ).GetVectorFromOrigin();
        ++nb;
        }

      opt += points->ElementAt( q_it.Value()->GetDestination() ).GetVectorFromOrigin();
      ++nn;
      ++q_it;
      }

    if ( nb )
      {
      for ( unsigned int kk = 0; kk < 3; ++kk )
        {
        opt[kk] = 0.75 * ipt[kk] + 0.125 * bpt[kk];
        }
      }
    else
      {
      OutputCoordType beta =
        ( 0.625 - ( 0.375 + 0.25 * vcl_cos(2.0 * vnl_math::pi / nn) ) * ( 0.375 + 0.25 * vcl_cos(2.0 * vnl_math::pi / nn) ) ) / nn;
      for ( unsigned int kk = 0; kk < 3; ++kk )
        {
        opt[kk] = ( 1.0 - nn * beta ) * ipt[kk] + beta * opt[kk];
        }
      }

    output->SetPoint(ptIt->Index(), opt);
    }
}
}
#endif
