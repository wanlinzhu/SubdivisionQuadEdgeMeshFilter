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
#ifndef __itkModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter_hxx
#define __itkModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter_hxx

#include "itkModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
template< typename TInputMesh, typename TOutputMesh >
void
ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CellSubdivision( OutputCellType *cell, OutputMeshType *output )
{
  if ( cell->GetType() != OutputCellType::POLYGON_CELL || cell->GetNumberOfPoints() != 3 )
  {
  return;
  }
	
  OutputPointIdentifier oldPointIdArray[3];
  OutputPointIdentifier newPointIdArray[3];

  OutputCoordType pointWeight[8] = {0.5, 0.5, 0.125, 0.125, -0.0625, -0.0625, -0.0625, -0.0625};
  OutputPointType pointArray[8];

  OutputPointIdIterator it = cell->PointIdsBegin();
  unsigned int          numberOfPoints = output->GetNumberOfPoints();
  unsigned int          n = 0;

  while ( it != cell->PointIdsEnd() )
    {
    oldPointIdArray[n++] = *it;
    ++it;
    }

  for ( unsigned int ii = 0; ii < 3; ++ii )
    {
    int jj = ( ii + 1 ) % 3;

    OutputQEType *edge = this->GetOutput()->FindEdge(oldPointIdArray[ii], oldPointIdArray[jj]);

    if ( this->m_EdgesPointIdentifier->IndexExists(edge) )
      {
      newPointIdArray[ii] = this->m_EdgesPointIdentifier->GetElement(edge);
      }
    else
      {
      OutputPointType outPoint;
      outPoint.Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);

      this->GetOutput()->GetPoint(oldPointIdArray[ii], &pointArray[0]);
      this->GetOutput()->GetPoint(oldPointIdArray[jj], &pointArray[1]);

      if ( edge->GetLnext() )
        {
        this->GetOutput()->GetPoint(edge->GetLnext()->GetDestination(), &pointArray[2]);

        if ( edge->GetLnext()->GetRprev() )
          {
          this->GetOutput()->GetPoint(edge->GetLnext()->GetRprev()->GetDestination(), &pointArray[4]);
          }
        else
          {
          pointArray[4].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
          }
        }
      else
        {
        pointArray[2].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      if ( edge->GetRprev() )
        {
        this->GetOutput()->GetPoint(edge->GetRprev()->GetDestination(), &pointArray[3]);
        if ( edge->GetRprev()->GetLnext() )
          {
          this->GetOutput()->GetPoint(edge->GetRprev()->GetLnext()->GetDestination(), &pointArray[5]);
          }
        else
          {
          pointArray[5].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
          }
        }
      else
        {
        pointArray[3].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      if ( edge->GetLprev() && edge->GetLprev()->GetRprev() )
        {
        this->GetOutput()->GetPoint(edge->GetLprev()->GetRprev()->GetDestination(), &pointArray[6]);
        }
      else
        {
        pointArray[6].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      if ( edge->GetRnext() && edge->GetRnext()->GetLnext() )
        {
        this->GetOutput()->GetPoint(edge->GetRnext()->GetLnext()->GetDestination(), &pointArray[7]);
        }
      else
        {
        pointArray[7].Fill(NumericTraits< typename OutputPointType::ValueType >::Zero);
        }

      for ( unsigned kk = 0; kk < 3; ++kk )
        {
        for ( unsigned int mm = 0; mm < 8; ++mm )
          {
          outPoint[kk] += pointWeight[mm] * pointArray[mm][kk];
          }
        }

      newPointIdArray[ii] = numberOfPoints++;
      this->m_EdgesPointIdentifier->InsertElement(edge, newPointIdArray[ii]);
      this->m_EdgesPointIdentifier->InsertElement(edge->GetSym(), newPointIdArray[ii]);
      output->SetPoint(newPointIdArray[ii], outPoint);
      }
    }

  output->AddFaceTriangle(oldPointIdArray[0], newPointIdArray[0], newPointIdArray[2]);
  output->AddFaceTriangle(newPointIdArray[0], oldPointIdArray[1], newPointIdArray[1]);
  output->AddFaceTriangle(newPointIdArray[1], oldPointIdArray[2], newPointIdArray[2]);
  output->AddFaceTriangle(newPointIdArray[0], newPointIdArray[1], newPointIdArray[2]);
}
}
#endif
