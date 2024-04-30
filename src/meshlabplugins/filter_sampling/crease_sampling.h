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
namespace vcg
{

template <class MeshType>
void UniformCreaseEdgeSampling(MeshType &m, // the mesh that has to be sampled
                     std::vector<typename MeshType::CoordType> &UniformEdgeSamples, // the vector that will contain the set of points sample on the crease edges (including corners)
                     typename MeshType::ScalarType &radius,  // the Poisson Disk Radius (used if sampleNum==0, setted if sampleNum!=0)
                     float creaseAngleDeg // the angle that defines the crease edges
)
{
// we search for all the edges where the difference between the normals of the two incident faces 
 // is greater than a given threshold

  MeshType edgeMesh; // will contains all the crease edges
  float cosThr = cosf( math::ToRad(creaseAngleDeg) );
  for(auto fi = m.face.begin(); fi != m.face.end(); ++fi)
  {
    auto &f = *fi;
    for(int i = 0; i < 3; i++)
    {
      float dot = f.N() * f.FFp(i)->N();
      if((dot < cosThr) && (f.FFp(i) > &f) )
        tri::Allocator<MeshType>::AddEdge(edgeMesh, f.P0(i), f.P1(i));
    }
  }

  tri::UpdateBounding<MeshType>::Box(edgeMesh); 
  tri::Clean<MeshType>::RemoveDuplicateVertex(edgeMesh);
  tri::Allocator<MeshType>::CompactEveryVector(edgeMesh);
  tri::Clean<MeshType>::SelectNonManifoldVertexOnEdgeMesh(edgeMesh);
  tri::Clean<MeshType>::SplitSelectedVertexOnEdgeMesh(edgeMesh);
  tri::Allocator<MeshType>::CompactEveryVector(edgeMesh);
  
  tri::TrivialSampler<MeshType> ts(UniformEdgeSamples);  
  tri::SurfaceSampling<MeshType, tri::TrivialSampler < MeshType> >::EdgeMeshUniform(edgeMesh, ts, radius); 

  printf("Sampled %i points\n",ts.SampleVec().size());
}

template <class MeshType>
void CreaseAwarePoissonDiskSampling(MeshType &m, // the mesh that has to be sampled
                     std::vector<typename MeshType::CoordType> &poissonSamples, // the vector that will contain the set of points
                     typename MeshType::ScalarType &radius,  // the Poisson Disk Radius (used if sampleNum==0, setted if sampleNum!=0)
                    float creaseAngleDeg, // the angle that defines the crease edges
                     unsigned int randSeed=0)

{
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  typedef tri::MeshSampler<MeshType> MontecarloSampler;
  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;

  int t0=clock();
  std::vector<typename MeshType::CoordType> samples;
  UniformCreaseEdgeSampling(m, samples, radius, creaseAngleDeg);
  MeshType meshCreaseEdgeSamples;
  tri::BuildMeshFromCoordVector(meshCreaseEdgeSamples,samples);
 
  MeshType MontecarloMesh; // it will contains the samples generated by the montecarlo sampling over the input mesh 
  tri::MeshSampler<MeshType> mcSampler(MontecarloMesh);
  tri::SurfaceSampling<MeshType, tri::MeshSampler<MeshType> >::Montecarlo(m, mcSampler, 50000);
  tri::UpdateBounding<MeshType>::Box(MontecarloMesh);

  typename tri::SurfaceSampling<MeshType, tri::MeshSampler < MeshType> >::PoissonDiskParam pp2;
  pp2.preGenMesh = &meshCreaseEdgeSamples;
  pp2.preGenFlag = true;
  pp2.bestSamplePoolSize = 20;
  
  tri::TrivialSampler<MeshType> pdSampler(poissonSamples);   
  tri::SurfaceSampling<MeshType, tri::MeshSampler < MeshType> >::PoissonDiskPruning(pdSampler, MontecarloMesh, radius, pp2);
}


}
