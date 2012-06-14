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
#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif

#include "itkModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter.h"
#include "itkLinearTriangleCellSubdivisionQuadEdgeMeshFilter.h"
#include "itkLoopTriangleCellSubdivisionQuadEdgeMeshFilter.h"
#include "itkSquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter.h"
#include "itkQuadEdgeMesh.h"
#include "itkQuadEdgeMeshParamMatrixCoefficients.h"
#if ( ITK_VERSION_MAJOR < 4 )
#include "itkQuadEdgeMeshSmoothing.h"
#else
#include "itkSmoothingQuadEdgeMeshFilter.h"
#endif
#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"

int main(int argc, char *argv[])
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputMeshFile  outputMeshFile subdivisionType Resolution non-uniform" << std::endl;
    std::cerr << " 0 : ModifiedButterfly " << std::endl;
    std::cerr << " 1 : Linear " << std::endl;
    std::cerr << " 2 : Loop " << std::endl;
    std::cerr << " 3 : Squarethree " << std::endl;
    return EXIT_FAILURE;
    }

  typedef float MeshPixelType;
  const unsigned int Dimension = 3;

  typedef itk::QuadEdgeMesh< MeshPixelType, Dimension > InputMeshType;
  typedef itk::QuadEdgeMesh< MeshPixelType, Dimension > OutputMeshType;

  typedef itk::CellSubdivisionQuadEdgeMeshFilter< InputMeshType, OutputMeshType >                    CellSubdivisionFilterType;
  typedef itk::ModifiedButterflyTriangleCellSubdivisionQuadEdgeMeshFilter< InputMeshType, OutputMeshType >   ModifiedButterflySubdivisionFilterType;
  typedef itk::LinearTriangleCellSubdivisionQuadEdgeMeshFilter< InputMeshType, OutputMeshType >      LinearSubdivisionFilterType;
  typedef itk::LoopTriangleCellSubdivisionQuadEdgeMeshFilter< InputMeshType, OutputMeshType >        LoopSubdivisionFilterType;
  typedef itk::SquareThreeTriangleCellSubdivisionQuadEdgeMeshFilter< InputMeshType, OutputMeshType > SquareThreeSubdivisionFilterType;

  typedef itk::VTKPolyDataReader< InputMeshType >  ReaderType;
  typedef itk::VTKPolyDataWriter< OutputMeshType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while reading the input file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  CellSubdivisionFilterType::Pointer subdivision;

  if ( argc >= 4 )
    {
    int type = std::atoi(argv[3]);

    switch ( type )
      {
      case 0:
        subdivision = ModifiedButterflySubdivisionFilterType::New().GetPointer();
        break;
      case 1:
        subdivision = LinearSubdivisionFilterType::New().GetPointer();
        break;
      case 2:
        subdivision = LoopSubdivisionFilterType::New().GetPointer();
        break;
      case 3:
        subdivision = SquareThreeSubdivisionFilterType::New().GetPointer();
        break;
      default:
        std::cerr << "Invalid subdivision type : " << type << std::endl;
        return EXIT_FAILURE;
      }
    }
  else
    {
    std::cerr << "You must have subdivision type " << std::endl;
    return EXIT_FAILURE;
    }

  if ( argc >= 5 )
    {
    unsigned int n = std::atoi(argv[4]);
    subdivision->SetResolutionLevels(n);
    }

  if ( argc >= 6 )
    {
    subdivision->UniformOff();
    subdivision->AddSubdividedCellId(0);
    subdivision->AddSubdividedCellId(1);
    subdivision->AddSubdividedCellId(2);
    subdivision->AddSubdividedCellId(3);
    subdivision->AddSubdividedCellId(5);
    subdivision->AddSubdividedCellId(6);
    subdivision->AddSubdividedCellId(9);
    }
  subdivision->SetInput( reader->GetOutput() );
  subdivision->Update();

  bool smoothing = true;
  if ( argc >= 7 )
  {
  smoothing = false;
  }

  if ( smoothing )
    {
#if ( ITK_VERSION_MAJOR < 4 )
    typedef itk::QuadEdgeMeshSmoothing< OutputMeshType, OutputMeshType > OutputMeshSmoothingFilterType;
#else
    typedef itk::SmoothingQuadEdgeMeshFilter< OutputMeshType, OutputMeshType > OutputMeshSmoothingFilterType;
#endif
    typedef itk::MatrixCoefficients< OutputMeshType >                          MatrixCoefficientsType;
    typedef itk::OnesMatrixCoefficients< OutputMeshType >                      OnesMatrixCoefficientsType;
    typedef itk::InverseEuclideanDistanceMatrixCoefficients< OutputMeshType >  InverseEuclideanDistanceMatrixCoefficientsType;
    typedef itk::ConformalMatrixCoefficients< OutputMeshType >                 ConformalMatrixCoefficientsType;
    typedef itk::AuthalicMatrixCoefficients< OutputMeshType >                  AuthalicMatrixCoefficientsType;
    typedef itk::IntrinsicMatrixCoefficients< OutputMeshType >                 IntrinsicMatrixCoefficientsType;
    typedef itk::HarmonicMatrixCoefficients< OutputMeshType >                  HarmonicMatrixCoefficientsType;

    OnesMatrixCoefficientsType             coef;
    OutputMeshSmoothingFilterType::Pointer meshSmoothingFilter = OutputMeshSmoothingFilterType::New();
    meshSmoothingFilter->SetInput( subdivision->GetOutput() );
    meshSmoothingFilter->SetCoefficientsMethod(&coef);
    meshSmoothingFilter->SetDelaunayConforming(1);
    meshSmoothingFilter->SetNumberOfIterations(1);
    meshSmoothingFilter->Update();

    subdivision->GetOutput()->Graft( meshSmoothingFilter->GetOutput() );
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput( subdivision->GetOutput() );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject & exp )
    {
    std::cerr << "Exception thrown while writting the output file " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
