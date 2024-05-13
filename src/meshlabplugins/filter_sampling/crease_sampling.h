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
#include "vcg/complex/algorithms/update/color.h"
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/point_sampling.h>
namespace vcg
{

template <class MeshType>
void UniformCreaseEdgeSampling(MeshType &m, // the mesh that has to be sampled
							   MeshType &UniformEdgeSampleMesh, // the vector that will contain the set of points sample on the crease edges (including corners)
							   typename MeshType::ScalarType &radius,  // the Poisson Disk Radius
							   float creaseAngleDeg, // the angle that defines the crease edges
							   bool includeBoundaryEdges=true
)
{
// we search for all the edges where the difference between the normals of the two incident faces 
 // is greater than a given threshold
	
	MeshType edgeMesh; // will contains all the crease edges
	float cosThr = cosf( math::ToRad(creaseAngleDeg) );
	int boundaryCnt=0;
	int creaseCnt=0;
	for(auto fi = m.face.begin(); fi != m.face.end(); ++fi)
	{
		auto &f = *fi;
		for(int i = 0; i < 3; i++)
		{
			float dot = f.N() * f.FFp(i)->N();
			if( ((dot < cosThr) && (f.FFp(i) > &f) )){
				creaseCnt++;
				tri::Allocator<MeshType>::AddEdge(edgeMesh, f.P0(i), f.P1(i));
				edgeMesh.edge.back().V(0)->N() = f.V0(i)->N();
				edgeMesh.edge.back().V(1)->N() = f.V1(i)->N();
			}
			else
			{			
				if(includeBoundaryEdges && face::IsBorder(f,i)) {
					boundaryCnt++;
					tri::Allocator<MeshType>::AddEdge(edgeMesh, f.P0(i), f.P1(i));
					edgeMesh.edge.back().V(0)->N() = f.V0(i)->N();
					edgeMesh.edge.back().V(1)->N() = f.V1(i)->N();
				}
			}
		}
  }
	
  //printf("Collected %i crease edges and %i boundary edges\n",creaseCnt,boundaryCnt);
  

  tri::UpdateBounding<MeshType>::Box(edgeMesh); 
  tri::Clean<MeshType>::RemoveDuplicateVertex(edgeMesh);
  tri::Allocator<MeshType>::CompactEveryVector(edgeMesh);
  tri::Clean<MeshType>::SelectNonManifoldVertexOnEdgeMesh(edgeMesh);
  int nonManifVertNum = tri::UpdateSelection<MeshType>::VertexCount(edgeMesh);
  printf("Selected %i non manifold vertex number\n",nonManifVertNum);
  tri::UpdateColor<MeshType>::PerVertexConstant(edgeMesh,Color4b::Red,true);
  tri::Clean<MeshType>::SplitSelectedVertexOnEdgeMesh(edgeMesh);
  
  tri::Allocator<MeshType>::CompactEveryVector(edgeMesh);
  tri::Clean<MeshType>::SelectCreaseVertexOnEdgeMeshEE(edgeMesh, math::ToRad(creaseAngleDeg));
  int creaseVertNum = tri::UpdateSelection<MeshType>::VertexCount(edgeMesh);
  printf("Selected %i crease vertex number\n",creaseVertNum);
  tri::Clean<MeshType>::SplitSelectedVertexOnEdgeMesh(edgeMesh);
  tri::Allocator<MeshType>::CompactEveryVector(edgeMesh);
  
  tri::MeshSampler<MeshType> ts(UniformEdgeSampleMesh);  
  tri::SurfaceSampling<MeshType, tri::MeshSampler < MeshType> >::EdgeMeshUniform(edgeMesh, ts, radius); 

  printf("Sampled %i points\n",edgeMesh.vn);
}

template <class MeshType>
void CreaseAwarePoissonDiskSampling(MeshType &m, // the mesh that has to be sampled
					 MeshType &edgeSampleMesh, // the point cloud that will contain the set of points on the edges
                     MeshType &poissonMesh, // the point cloud that will contain the set of points of the poisson disk sampling
                     typename MeshType::ScalarType &radius,  // the Poisson Disk Radius (used if sampleNum==0, setted if sampleNum!=0)
                     float creaseAngleDeg, // the angle that defines the crease edges
                     unsigned int randSeed=0)

{
  typedef tri::TrivialSampler<MeshType> BaseSampler;
  typedef tri::MeshSampler<MeshType> MontecarloSampler;
  typename tri::SurfaceSampling<MeshType, BaseSampler>::PoissonDiskParam pp;
  
  int estimatedSampleNum = tri::SurfaceSampling<MeshType, BaseSampler>::ComputePoissonSampleNum(m, radius);
  
  int t0=clock();
  UniformCreaseEdgeSampling(m, edgeSampleMesh, radius, creaseAngleDeg);
 
  MeshType MontecarloMesh; // it will contains the samples generated by the montecarlo sampling over the input mesh 
  tri::MeshSampler<MeshType> mcSampler(MontecarloMesh);
  tri::SurfaceSampling<MeshType, tri::MeshSampler<MeshType> >::Montecarlo(m, mcSampler, estimatedSampleNum*20);
  tri::UpdateBounding<MeshType>::Box(MontecarloMesh);

  typename tri::SurfaceSampling<MeshType, tri::MeshSampler < MeshType> >::PoissonDiskParam pp2;
  pp2.preGenMesh = &edgeSampleMesh;
  pp2.preGenFlag = true;
  pp2.bestSamplePoolSize = 20;
  
  tri::MeshSampler<MeshType> pdSampler(poissonMesh);   
  tri::SurfaceSampling<MeshType, tri::MeshSampler < MeshType> >::PoissonDiskPruning(pdSampler, MontecarloMesh, radius, pp2);
}


}
