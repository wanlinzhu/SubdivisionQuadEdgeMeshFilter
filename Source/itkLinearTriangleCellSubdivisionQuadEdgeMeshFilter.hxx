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
#ifndef __itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter_hxx
#define __itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter_hxx

#include "itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter.h"

namespace itk
{
template< typename TInputMesh, typename TOutputMesh >
void
LinearTriangleCellSubdivisionQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
::CellSubdivision( OutputCellType *cell, OutputMeshType *output )
{
  if ( cell->GetType() != OutputCellType::POLYGON_CELL || cell->GetNumberOfPoints() != 3 )
    {
    itkExceptionMacro(<<" The input cell is not a triangle cell");
    }

  OutputPointIdentifier oldPointIdArray[3];
  OutputPointIdentifier newPointIdArray[3];
  OutputPointType       pointArray[3];

  OutputPointIdIterator it = cell->PointIdsBegin();
  OutputPointIdentifier n = 0;
  OutputPointIdentifier numberOfPoints = output->GetNumberOfPoints();

  while ( it != cell->PointIdsEnd() )
    {
    oldPointIdArray[n] = *it;
    this->GetOutput()->GetPoint(*it, &pointArray[n]);
    ++it;
    ++n;
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
    else
      {
      OutputPointType outPoint;
      outPoint.SetToMidPoint(pointArray[ii], pointArray[jj]);

      newPointIdArray[ii] = numberOfPoints++;

      this->m_EdgesPointIdentifier->InsertElement(edge, newPointIdArray[ii]);

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
