/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005                                                \/)\/    *
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

#include "filter_ransac.h"
#include "ransac_plane_fitting.h"

using namespace vcg;

/**
 * @brief Constructor usually performs only two simple tasks of filling the two lists
 *  - typeList: with all the possible id of the filtering actions
 *  - actionList with the corresponding actions. If you want to add icons to
 *  your filtering actions you can do here by construction the QActions accordingly
 */
FilterRansacPlugin::FilterRansacPlugin() 
{ 
	typeList = {FP_FIT_BEST_PLANE,FP_DECOMPOSE_IN_PLANES};

	for(ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

FilterRansacPlugin::~FilterRansacPlugin()
{
}

QString FilterRansacPlugin::pluginName() const
{
    return "FilterRansac";
}

/**
 * @brief ST() must return the very short string describing each filtering action
 * (this string is used also to define the menu entry)
 * @param filterId: the id of the filter
 * @return the name of the filter
 */
QString FilterRansacPlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_FIT_BEST_PLANE :
		return "Ransac Best Plane Fitting";
    case FP_DECOMPOSE_IN_PLANES :
        return "Ransac Plane Segmentation";
	default :
		assert(0);
		return QString();
	}
}

/**
 * @brief FilterRansacPlugin::pythonFilterName if you want that your filter should have a different
 * name on pymeshlab, use this function to return its python name.
 * @param f
 * @return
 */
QString FilterRansacPlugin::pythonFilterName(ActionIDType f) const
{
	switch(f) {
	case FP_FIT_BEST_PLANE :
		return "computeXXXXXX";
    case FP_DECOMPOSE_IN_PLANES :
        return "computeXXXXXX";
	default :
		assert(0);
		return QString();
	}
}


/**
 * @brief // Info() must return the longer string describing each filtering action
 * (this string is used in the About plugin dialog)
 * @param filterId: the id of the filter
 * @return an info string of the filter
 */
 QString FilterRansacPlugin::filterInfo(ActionIDType filterId) const
{
	switch(filterId) {
    case FP_FIT_BEST_PLANE :
        return "Find the plane that best fit a given mesh and select the vertexes belonging to it";
    case FP_DECOMPOSE_IN_PLANES :
        return "Find the set of planes that fit well a given mesh.";
	default :
		assert(0);
		return "Unknown Filter";
	}
}

 /**
 * @brief The FilterClass describes in which generic class of filters it fits.
 * This choice affect the submenu in which each filter will be placed
 * More than a single class can be chosen.
 * @param a: the action of the filter
 * @return the class od the filter
 */
FilterRansacPlugin::FilterClass FilterRansacPlugin::getClass(const QAction *a) const
{
	switch(ID(a)) {
	case FP_FIT_BEST_PLANE :
        return FilterPlugin::Selection;
    case FP_DECOMPOSE_IN_PLANES :
        return FilterPlugin::FaceColoring;
	default :
		assert(0);
		return FilterPlugin::Generic;
	}
}

/**
 * @brief FilterRansacPlugin::filterArity
 * @return
 */
FilterPlugin::FilterArity FilterRansacPlugin::filterArity(const QAction*) const
{
	return SINGLE_MESH;
}

/**
 * @brief FilterRansacPlugin::getPreConditions
 * @return
 */
int FilterRansacPlugin::getPreConditions(const QAction*a) const
{
    switch(ID(a)) {
    case FP_FIT_BEST_PLANE :
        return MeshModel::MM_NONE;
    case FP_DECOMPOSE_IN_PLANES :
        return MeshModel::MM_FACENUMBER;
    default :
        assert(0);
    }
    return 0;
}

/**
 * @brief FilterRansacPlugin::postCondition
 * @return
 */
int FilterRansacPlugin::postCondition(const QAction*) const
{
    return MeshModel::MM_ALL;
}

/**
 * @brief This function define the needed parameters for each filter. Return true if the filter has some parameters
 * it is called every time, so you can set the default value of parameters according to the mesh
 * For each parameter you need to define,
 * - the name of the parameter,
 * - the default value
 * - the string shown in the dialog
 * - a possibly long string describing the meaning of that parameter (shown as a popup help in the dialog)
 * @param action
 * @param m
 * @param parlst
 */
RichParameterList FilterRansacPlugin::initParameterList(const QAction *action,const MeshModel &m)
{
	RichParameterList parlst;
	switch(ID(action)) {
    case FP_FIT_BEST_PLANE :
        parlst.addParam(RichInt("RansacIteration", 500, "Ransac Iteration", "The number of attempt that the ransac algorithm makes to sample three points and guess a plane."));
        parlst.addParam(RichPercentage("PlaneThr",m.cm.bbox.Diag()/100.0f,0.0f,m.cm.bbox.Diag(),"Plane Threshold","The plane distance threshold used by ransac to decide if a point belong or not to a plane"));
        parlst.addParam(RichInt("InlierEstimationNum",100,"Inlier Estim. Num","The number of vertices that are used to estimate the fraction of points that fits a plane."));
        parlst.addParam(RichBool("ExcludeSel",true,"Exclude Selection","If the "));
		break;
    case FP_DECOMPOSE_IN_PLANES:
        break;
    default :
		assert(0);
	}
	return parlst;
}

/**
 * @brief The Real Core Function doing the actual mesh processing.
 * @param action
 * @param md: an object containing all the meshes and rasters of MeshLab
 * @param par: the set of parameters of each filter
 * @param cb: callback object to tell MeshLab the percentage of execution of the filter
 * @return true if the filter has been applied correctly, false otherwise
 */
std::map<std::string, QVariant> FilterRansacPlugin::applyFilter(const QAction * action, const RichParameterList & parameters, MeshDocument &md, unsigned int& /*postConditionMask*/, CallBackPos *cb)
{
	switch(ID(action)) {
	case FP_FIT_BEST_PLANE :        
    {
        md.mm()->updateDataMask(MeshModel::MM_VERTCOLOR);
       
        CMeshO &m = md.mm()->cm;
        tri::RansacPlaneFitting<CMeshO>::Parameters params;
        params.planeDistanceThresholdAbs = parameters.getAbsPerc("PlaneThr");
        params.maxIterations = parameters.getInt("RansacIteration");
        params.onlySelected=true;
        params.quickInlierEstimationSamples=parameters.getInt("InlierEstimationNum");
        tri::UpdateSelection<CMeshO>::VertexAll(m);
        tri::RansacPlaneFitting<CMeshO> rpf;
        rpf.Init(m,params);
        
        Plane3m bestPlane1;
        Plane3m bestPlane2; 
        float bestInliersFraction; 
        rpf.SearchBestFittingPlane(bestPlane1);
        rpf.ComputeInlierFraction(bestPlane1,bestInliersFraction);
        tri::UpdateColor<CMeshO>::PerVertexConstant(m,Color4b::White);
        rpf.ColorizeMeshByPlaneDistance(bestPlane1,Color4b::Green);
        
        log("fitted plane  vertices of a mesh of %i vertices found %5.3f inliers",m.vn, bestInliersFraction);
        log("done %i iterations and %i improv and %i guess",params.maxIterations,rpf.s.improvedCnt, rpf.s.correctGuessCnt );
        
    }
		break;
    case FP_DECOMPOSE_IN_PLANES:
        break;    
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterRansacPlugin)
