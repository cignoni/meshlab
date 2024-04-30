/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/
#pragma once

#include <deque>
#include <functional>
#include <random>
#include <vcg/simplex/face/pos.h>
#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/update/quality.h>

namespace vcg{
namespace tri{

template <class MeshType>
class RansacPlaneFitting{
  public:
  typedef typename MeshType::VertexType VertexType;
  typedef typename MeshType::ScalarType  ScalarType;
  typedef typename MeshType::FacePointer FacePointer;
  typedef Plane3<ScalarType> Plane3m;
  
  class Stats
  {
  public:
      int improvedCnt;
      int correctGuessCnt;
      
      void reset()
      {
          improvedCnt=0;
          correctGuessCnt=0;
      }
  };
  
  class Parameters
  {
  public:
      
      float  planeDistanceThresholdPerc = 0.1f; // distance threshold for inlier/outlier classification
      float  planeDistanceThresholdAbs = -1.0f; // distance threshold for inlier/outlier classification (neg values means not initialized) 
      float sampleDistanceThresholdPerc = 0.1f; // distance used to be sure that sampling points are far apart
      float sampleDistanceThresholdAbs = -1.0f; // distance used to be sure that sampling points are far apart
      int  maxIterations = 1000; // the number of planes we sample to search the best one. The higher the number the more accurate guess we get.
      float  minInlierFraction = 0.5f; // minimum fraction of inliers to accept a plane
      int quickInlierEstimationSamples = 500; // number of samples used for the quick estimate of the inlier fraction
      bool onlySelected = false; // if true, only selected vertices are used for the plane fitting
      
      Parameters() {}
      
      void Init(MeshType &mesh)
      {
          // let fill the permuation vector with the indices of the selected vertices
          
          
          if(planeDistanceThresholdAbs < 0)
              planeDistanceThresholdAbs=   mesh.bbox.Diag()*  planeDistanceThresholdPerc;
          if(sampleDistanceThresholdAbs < 0)
              sampleDistanceThresholdAbs= mesh.bbox.Diag()*sampleDistanceThresholdPerc;
      }
  };

void Init(MeshType &mesh, Parameters &param)
{
      p=param;
      m=&mesh;
      for(int i=0;i<mesh.VN();i++)
          if(!p.onlySelected || mesh.vert[i].IsS())
              permutationVec.push_back(i);
      
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(permutationVec.begin(),permutationVec.end(),g);
}

private:
MeshType *m;
Parameters p;
std::vector<int> permutationVec;
    
public:
Stats s;
    
static math::MarsenneTwisterRNG &rnd() 
{
  static math::MarsenneTwisterRNG rnd;
  return rnd;
}

void SearchBestFittingPlane(Plane3m &bestFittingPlane)
{
  Plane3m fitPlane;
  float inliersFraction = 0;
  int ransacIterations = 0;
  float bestInliersFraction = 0;
  s.reset();
  
  while (ransacIterations < p.maxIterations)
  {
    GeneratePlane(fitPlane);
    EstimateInlierFraction(fitPlane,  inliersFraction);
    
    if (inliersFraction > bestInliersFraction * 0.8f)
    {
        ComputeInlierFraction(fitPlane, inliersFraction);  // this is the precise estimate of the inlier fraction
        ++s.correctGuessCnt;
    }
    
    if (inliersFraction > bestInliersFraction)
    {
        ++s.improvedCnt;
      bestInliersFraction = inliersFraction;
      bestFittingPlane = fitPlane;
    }
    ransacIterations++;
  }
}

 void GeneratePlane(Plane3m &fittingPlane)
{
  // First I take three random points over my mesh
  // checking that they are far apart and not aligned
 
  int p0,p1,p2;
  p0 = permutationVec[rnd().generate(permutationVec.size())];
  while(true)
  {
    p1 = permutationVec[rnd().generate(permutationVec.size())];
    if (Distance(m->vert[p0].P(),m->vert[p1].P()) > p.sampleDistanceThresholdAbs)  
      break;
  }
 while(true)
  {
    p2 = permutationVec[rnd().generate(permutationVec.size())];
    
    if ( (Distance(m->vert[p0].P(),m->vert[p2].P()) > p.sampleDistanceThresholdAbs) &&
    (Distance(m->vert[p1].P(),m->vert[p2].P()) > p.sampleDistanceThresholdAbs) )
      break;
  }
 // Then I compute the plane passing through them
  fittingPlane.Init(m->vert[p0].P(),m->vert[p1].P(),m->vert[p2].P());
}

// This function returns the fraction of the dataset that matches with the current hypotesis
// In this case the fraction is intended as the fraction of the total surface of the mesh that
// is composed by triangles with all three vertexes closer to the plane than a given trheshold 
 void ComputeInlierFractionAsArea(const Plane3m &fittingPlane, float &inlierFraction)
{
  
}


// This function returns the fraction of the dataset that matches with the current hypotesis
// in this case the fraction is intended as the number of points closer to the plane than the given threshold

 void ComputeInlierFraction(const Plane3m &fittingPlane, float &inlierFraction)
{
  EstimateInlierFraction(fittingPlane, inlierFraction, true);
}

// This function computes an estimate of the inlier fraction of a plane by choosing a random subset of points to check
 void EstimateInlierFraction(const Plane3m &fittingPlane, float &inlierFraction, bool useAllPoints=false)
{
 int pointsOnPlane = 0;
  int checkingPoints = permutationVec.size();
 if(!useAllPoints)
    checkingPoints = std::min(checkingPoints,p.quickInlierEstimationSamples);
 
  for(int i=0;i<checkingPoints;++i)
  {
    int ip = permutationVec[i];
    if (fabs(SignedDistancePlanePoint(fittingPlane,m->vert[ip].P())) < p.planeDistanceThresholdAbs)
        pointsOnPlane++;
  }
  inlierFraction = (float)pointsOnPlane/(float)checkingPoints; 
}

 void DeselectMeshByPlaneDistance(Plane3m &fittingPlane)
{
  for (int i = 0; i < m->vert.size(); i++)
  {
    if (fabs(SignedDistancePlanePoint(fittingPlane,m->vert[i].P())) < p.planeDistanceThresholdAbs)
      m->vert[i].ClearS();
  }
}

 void ColorizeMeshByPlaneDistance(Plane3m &fittingPlane, Color4b c=Color4b::Red)
{
  for (int i = 0; i < m->vert.size(); i++)
  {
    if (fabs(SignedDistancePlanePoint(fittingPlane,m->vert[i].P())) < p.planeDistanceThresholdAbs)
      m->vert[i].C() = c;
  }
}



};
}   // end namespace tri
}   // end namespace vcg
