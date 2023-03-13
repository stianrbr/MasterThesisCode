// Parametric modelling of the INO WINDMOOR 12 MW FOWT
// author: Stian Brurås


// This code models the wetted part of the hull

GenieRules.Units.setInputUnit(Angle, "deg");

GenieRules.Tolerances.useTolerantModelling = true;

GenieRules.Compatibility.version = "V8.3-31";
GenieRules.Tolerances.useTolerantModelling = true;
GenieRules.Tolerances.angleTolerance = 2 deg;

GenieRules.Meshing.autoSimplifyTopology = true;
GenieRules.Meshing.eliminateInternalEdges = true;
GenieRules.BeamCreation.DefaultCurveOffset = ReparameterizedBeamCurveOffset();
GenieRules.Transformation.DefaultConnectedCopy = false;

Draft = 20.0m; // Draft of floater
Col_diameter = 13.32m; // Diameter of column
Col_height = 35.0m; // Column height
Pon_width = 11.1m; // Width of pontoon
Pon_height =  2.24m; // Height of pontoon
CC_distance = 77.23m; // Center-Center distance columns

Deck_height = 3.5m;
Deck_width = 3.5m;

fairlead_depth = 14.00 m; // Fairlead_depth

damage_lower_extent = 3m;
damage_upper_extent = 5m;
horizontal_extent_sector = 0.125; // fraction of circumference
penetration = 1.5m;

t1 = 40.0 mm; // Steel thickness, external plates
t2 = 25.0 mm; // Steel thickness, internal plates
Tck1 = Thickness(t1);
Tck2 = Thickness(t2);
Tck1.setDefault();


Steel = MaterialLinear(360e6, 7850 kg/m^3, 2.1e+11 Pa, 0.3, 1.2e-05 delC^-1, 0.03 N*s/m);
Steel.setDefault();


Mesh_length = 0.8m; // Size of elements
Mesh_radius = 300 m;  // Radius for free surface mesh
FS_mesh_length = 5.0 m;  // Mesh size for free surface mesh
mesh_lid = 0.5 m;  // Mesh size for internal lids

Md_def = MeshDensity(Mesh_length);
Md_def.enforceDensity = true;
Lid_density = MeshDensity(mesh_lid);

mesh_square_ratio = 0.5;
Half_pon_height = NumberOfElements(Math.Ceil((Pon_height/2)/Mesh_length));

p_straight = 58.0 m;
mesh_rad  = 6.0 m;


Md_def.setDefault();

// Wet Surfaces
outer_wet = WetSurface();
wetted_surface = Set();


//Meshing rules
GenieRules.Meshing.elementType = mp1stOrder;
GenieRules.Meshing.superElementType = 1;
GenieRules.Meshing.autoSimplifyTopology = false;
GenieRules.Meshing.autoSplitPeriodicGeometry = false;
GenieRules.Meshing.preference(mpPreferRectangularMesh, false);
GenieRules.Meshing.preference(mpAllowTriangularElements, false);
GenieRules.Meshing.preference(mpPreferPointMassAsNodeMass, true);
GenieRules.Meshing.preference(mpUseDrillingElements, false);
GenieRules.Meshing.preference(mpUseEccentricHinges, true);
GenieRules.Meshing.eliminateInternalEdges = false;
GenieRules.Meshing.eliminateInternalVertices = true;
GenieRules.Meshing.preference(mpIncludeUnusedProperties, false);
GenieRules.Meshing.preference(mpUseLongLoadcaseNames, false);
GenieRules.Meshing.preference(mpUseLongSetNames, false);
GenieRules.Meshing.preference(mpUseLongPropertyNames, false);
GenieRules.Meshing.preference(mpMeshDensityRounded, true);
GenieRules.Meshing.scantlings = msGross;
GenieRules.Meshing.ignoreEccentricities = false;
GenieRules.Meshing.useCocentricBeams = false;
//GenieRules.Meshing.faceMeshStrategy = SesamQuadMesher;
GenieRules.Meshing.edgeMeshStrategy = UniformDistributionEdge;
GenieRules.Meshing.activate(mpMaxAngle, mpFail, true);
GenieRules.Meshing.setLimit(mpMaxAngle, mpFail, 179 deg);
GenieRules.Meshing.activate(mpMaxAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxAngle, mpSplit, 165 deg);
GenieRules.Meshing.activate(mpMinAngle, mpFail, false);
GenieRules.Meshing.setLimit(mpMinAngle, mpFail, 1 deg);
GenieRules.Meshing.activate(mpMinAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMinAngle, mpSplit, 15 deg);
GenieRules.Meshing.activate(mpMaxRelativeJacobi, mpFail, false);
GenieRules.Meshing.setLimit(mpMaxRelativeJacobi, mpFail, 10);
GenieRules.Meshing.activate(mpMaxRelativeJacobi, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxRelativeJacobi, mpSplit, 5);
GenieRules.Meshing.activate(mpMinNormalizedJacobi, mpFail, false);
GenieRules.Meshing.setLimit(mpMinNormalizedJacobi, mpFail, 0);
GenieRules.Meshing.activate(mpMinNormalizedJacobi, mpSplit, false);
GenieRules.Meshing.setLimit(mpMinNormalizedJacobi, mpSplit, 0.2);
GenieRules.Meshing.activate(mpMinEdge, false);
GenieRules.Meshing.setLimit(mpMinEdge, 0.1);
GenieRules.Meshing.activate(mpMaxChord, false);
GenieRules.Meshing.setLimit(mpMaxChord, 0.2);
GenieRules.Meshing.activate(mpMaxTwistAngle, mpFail, false);
GenieRules.Meshing.setLimit(mpMaxTwistAngle, mpFail, 30 deg);
GenieRules.Meshing.activate(mpMaxTwistAngle, mpSplit, false);
GenieRules.Meshing.setLimit(mpMaxTwistAngle, mpSplit, 10 deg);
GenieRules.Meshing.IgnoreStructureNotPartOfSubsetWhenIdealizing = true;
GenieRules.Meshing.useUniformizedFaceParameterization = true;

//Tolerances Rules
GenieRules.Tolerances.angleTolerance = 2 deg;
GenieRules.Tolerances.pointTolerance = 0.01 m;
GenieRules.Tolerances.useTolerantModelling = true;


LC1 = DummyHydroLoadCase();
LC1.setFemLoadcase(1);
LC1.designCondition(lcOperating);
LC1.wetSurface = outer_wet;
//Analyses
Meshing_panels = Analysis(true);
Meshing_panels.add(MeshActivity());

//
// Geometrical modelling
//


// Column A points

P_A_Cx = -1/3*CC_distance*Math.cos(30*Math.PI/180);
P_A_Cy = CC_distance/2;
P_A_c0 = Point(P_A_Cx, P_A_Cy, -Draft);  // Center point of column A at keel
P_A_c1 = Point(P_A_Cx, P_A_Cy, 0);  // Center point of column A at MWL
P_A_c2 = Point(P_A_Cx, P_A_Cy, -Draft+Col_height);  // Center point of column A at top
P_A_cp = Point(P_A_Cx, P_A_Cy, -Draft+Pon_height);  // Center point of column A at top

sq_length = mesh_square_ratio/2*Col_diameter;
P_A_F1 = Point(P_A_Cx, P_A_Cy+sq_length, -Draft);
P_A_F2 = Point(P_A_Cx+sq_length, P_A_Cy+sq_length, -Draft);
P_A_F3 = Point(P_A_Cx+sq_length, P_A_Cy, -Draft);


P_A_Cir1 = Point(P_A_Cx+Math.cos(90*Math.PI/180)*Col_diameter/2, P_A_Cy+Math.sin(90*Math.Pi/180)*Col_diameter/2, -Draft);
P_A_Cir2 = Point(P_A_Cx+Math.cos(67.5*Math.PI/180)*Col_diameter/2, P_A_Cy+Math.sin(67.5*Math.Pi/180)*Col_diameter/2, -Draft);
P_A_Cir3 = Point(P_A_Cx+Math.cos(45*Math.PI/180)*Col_diameter/2, P_A_Cy+Math.sin(45*Math.Pi/180)*Col_diameter/2, -Draft);
P_A_Cir4 = Point(P_A_Cx+Math.cos(22.5*Math.PI/180)*Col_diameter/2, P_A_Cy+Math.sin(22.5*Math.Pi/180)*Col_diameter/2, -Draft);
P_A_Cir5 = Point(P_A_Cx+Math.cos(0*Math.PI/180)*Col_diameter/2, P_A_Cy+Math.sin(0*Math.Pi/180)*Col_diameter/2, -Draft);



// Column B points

P_B_Cx = 2/3*CC_distance*Math.cos(30*Math.PI/180);
P_B_Cy = 0m;
P_B_c0 = Point(P_B_Cx, P_B_Cy, -Draft); // Center point of column A at keel

P_B_c1 = Point(P_B_Cx, P_B_Cy, 0);  // Center point of column A at MWL

P_B_c2 = Point(P_B_Cx, P_B_Cy, -Draft+Col_height);  // Center point of column A at top
P_B_cp = Point(P_B_Cx, P_B_Cy, -Draft+Pon_height);  // Center point of column A at top


phi = Math.asin((Pon_width/2)/(Col_diameter/2))*180/Math.PI;  // Angle for pontoon-column intersection

el_elliptic = NumberOfElements(Math.Round(((phi+30)/360*Math.PI*Col_diameter)/Mesh_length));

// Pontoon 1 points

P_A1_cx = P_A_Cx + (Col_diameter/2*Math.cos(-90*Math.PI/180));
P_A1_cy = P_A_Cy + (Col_diameter/2*Math.sin(-90*Math.PI/180));
P_A1_c0 = Point(P_A1_cx, P_A1_cy, -Draft);  // Center point of pontoon 1 @ column A at keel
P_A1_c1 = P_A1_c0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Center point of pontoon 1 @ column A at pontoon top

P_A1_lx = P_A_Cx + (Col_diameter/2*Math.cos((-90-phi)*Math.PI/180));
P_A1_ly = P_A_Cy + (Col_diameter/2*Math.sin((-90-phi)*Math.PI/180));
P_A1_l0 = Point(P_A1_lx, P_A1_ly, -Draft);  // Left point of pontoon 1 @ column A at keel
P_A1_l1 = P_A1_l0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Left point of pontoon 1 @ column A at pontoon top

P_A1_rx = P_A_Cx + (Col_diameter/2*Math.cos((-90+phi)*Math.PI/180));
P_A1_ry = P_A_Cy + (Col_diameter/2*Math.sin((-90+phi)*Math.PI/180));
P_A1_r0 = Point(P_A1_rx, P_A1_ry, -Draft);  // Right point of pontoon 1 @ column A at keel
P_A1_r1 = P_A1_r0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Right point of pontoon 1 @ column A at pontoon top


P_C1_cx = P_A_Cx;
P_C1_cy = 0m;
P_C1_c0 = Point(P_C1_cx, P_C1_cy, -Draft);  // Center point of pontoon 1 @ y-axis at keel
P_C1_c1 = P_C1_c0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Center point of pontoon 1 @ y-axis at pontoon top

P_C1_lx = P_A_Cx-Pon_width/2;
P_C1_ly = 0m;
P_C1_l0 = Point(P_C1_lx, P_C1_ly, -Draft);  // Left point of pontoon 1 @ y-axis at keel
P_C1_l1 = P_C1_l0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Left point of pontoon 1 @ y-axis at pontoon top

P_C1_rx = P_A_Cx+Pon_width/2;
P_C1_ry = 0m;
P_C1_r0 = Point(P_C1_rx, P_C1_ry, -Draft);  // Right point of pontoon 1 @ y-axis at keel
P_C1_r1 = P_C1_r0.CopyTranslate(Vector3d(0, 0, Pon_height)); // Right point of pontoon 1 @ y-axis at pontoon top


// Pontoon 2 points

P_A2_cx = P_A_Cx + (Col_diameter/2*Math.cos(-30*Math.PI/180));
P_A2_cy = P_A_Cy + (Col_diameter/2*Math.sin(-30*Math.PI/180));
P_A2_c0 = Point(P_A2_cx, P_A2_cy, -Draft);  // Center point of pontoon 2 @ column A at keel
P_A2_c1 = P_A2_c0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Center point of pontoon 2 @ column A at pontoon top

P_A2_lx = P_A_Cx + (Col_diameter/2*Math.cos((-30-phi)*Math.PI/180));
P_A2_ly = P_A_Cy + (Col_diameter/2*Math.sin((-30-phi)*Math.PI/180));
P_A2_l0 = Point(P_A2_lx, P_A2_ly, -Draft);  // Left point of pontoon 2 @ column A at keel
P_A2_l1 = P_A2_l0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Left point of pontoon 2 @ column A at pontoon top

P_A2_rx = P_A_Cx + (Col_diameter/2*Math.cos((-30+phi)*Math.PI/180));
P_A2_ry = P_A_Cy + (Col_diameter/2*Math.sin((-30+phi)*Math.PI/180));
P_A2_r0 = Point(P_A2_rx, P_A2_ry, -Draft);  // Right point of pontoon 2 @ column A at keel
P_A2_r1 = P_A2_r0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Right point of pontoon 2 @ column A at pontoon top

P_B2_cx = P_B_Cx + Col_diameter/2*Math.cos(150*Math.PI/180);
P_B2_cy = P_B_Cy + Col_diameter/2*Math.sin(150*Math.PI/180);
P_B2_c0 = Point(P_B2_cx, P_B2_cy, -Draft);  // Center point of pontoon 2 @ column B at keel
P_B2_c1 = P_B2_c0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Center point of pontoon 2 @ column B at pontoon top

P_B2_lx = P_B_Cx + Col_diameter/2*Math.cos((150+phi)*Math.PI/180);
P_B2_ly = P_B_Cy + Col_diameter/2*Math.sin((150+phi)*Math.PI/180);
P_B2_l0 = Point(P_B2_lx, P_B2_ly, -Draft);  // Left point of pontoon 2 @ column B at keel
P_B2_l1 = P_B2_l0.CopyTranslate(Vector3d(0, 0, Pon_height)); // Left point of pontoon 2 @ column B at pontoon top


P_B2_rx = P_B_Cx + Col_diameter/2*Math.cos((150-phi)*Math.PI/180);
P_B2_ry = P_B_Cy + Col_diameter/2*Math.sin((150-phi)*Math.PI/180);
P_B2_r0 = Point(P_B2_rx, P_B2_ry, -Draft);  // Right point of pontoon 2 @ column B at keel
P_B2_r1 = P_B2_r0.CopyTranslate(Vector3d(0, 0, Pon_height));  // Right point of pontoon 2 @ column B at pontoon top

// Finding intersection points between pontoons and axes

m2 = (P_B2_ly-P_A2_ly)/(P_B2_lx-P_A2_lx);  // Gradient of line describing pontoon 2
c2 = ((P_A2_ly-P_B2_ly)/(P_B2_lx-P_A2_lx))*P_A2_lx + P_A2_ly;  // Intersect pf line describing pontoon 2

P_Yintersect2 = Point(-c2/m2,0m, -Draft);  // Intersect between pontoon 2 and y-axis
P_Yintersect2_1 = P_Yintersect2.CopyTranslate(Vector3d(0, 0, Pon_height));
P_BYintersect = Point(P_B_Cx-Col_diameter/2, P_B_Cy, -Draft);  // Intersect between column B and y-axis
P_BYintersect_1 = P_BYintersect.CopyTranslate(Vector3d(0, 0, Pon_height));
P_12intersect = Point(P_C1_rx, P_C1_rx*m2+c2,-Draft); // Intersect between pontoon 1 and 2
P_12intersect_1 = P_12intersect.CopyTranslate(Vector3d(0, 0, Pon_height));
P_A_12intersect = Point(P_A_Cx + (Col_diameter/2*Math.cos(-60*Math.PI/180)), P_A_Cy + (Col_diameter/2*Math.sin(-60*Math.PI/180)), -Draft);
P_A_12intersect_1 = P_A_12intersect.CopyTranslate(Vector3d(0, 0, Pon_height)); // Intersect between pontoon 1 and 2 radially towards column A_A_c2


l_inter = (-c2/m2-(P_B_Cx-Col_diameter/2));
P_eliptic_2Ar = P_A2_r0.copyTranslate(Vector3d(-l_inter*Math.cos(-30*Math.PI/180),-l_inter*Math.sin(-30*Math.PI/180) ,0));
P_eliptic_2Br = P_B2_r0.copyTranslate(Vector3d(-l_inter*Math.cos(150*Math.PI/180),-l_inter*Math.sin(150*Math.PI/180) ,0));
P_eliptic_1Al = P_A1_l0.copyTranslate(Vector3d(0, l_inter,0));

el_inter = NumberOfElements(Math.Round(-l_inter/Mesh_length));

//
// Creating lines for pontoons
//

//Pontoon 1
L_1_l0 = CreateLineTwoPoints(P_A1_l0, P_C1_l0);
L_1_l1 = L_1_l0.CopyTranslate(Vector3d(0,0,Pon_height));
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_1_r0 = CreateLineTwoPoints(P_12intersect, P_C1_r0);
}
else
{
L_1_r0 = CreateLineTwoPoints(P_A1_r0, P_C1_r0);
}
L_1_r1 = L_1_r0.CopyTranslate(Vector3d(0,0,Pon_height));
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_12_intersect0 = CreateLineTwoPoints(P_12intersect, P_A_12intersect);
L_12_intersect1 = L_12_intersect0.CopyTranslate(Vector3d(0,0,Pon_height));
}
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_1A_0 = CreateCircularArcFromThreePoints(P_A_12intersect, P_A1_c0, P_A1_l0);
}
else
{
L_1A_0 = CreateCircularArcFromThreePoints(P_A1_r0, P_A1_c0, P_A1_l0);
}
L_1A_1 = L_1A_0.CopyTranslate(Vector3d(0,0,Pon_height));

L_1Y_0 = CreateLineTwoPoints(P_C1_l0, P_C1_r0);
L_1Y_1 = L_1Y_0.CopyTranslate(Vector3d(0,0,Pon_height));

L_C1_lv = CreateLineTwoPoints(P_C1_l0, P_C1_l1);
L_C1_rv = CreateLineTwoPoints(P_C1_r0, P_C1_r1);
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_12_intersect_v = CreateLineTwoPoints(P_12intersect, P_12intersect_1);
L_A12_intersect_v = CreateLineTwoPoints(P_A_12intersect, P_A_12intersect_1);
}
else
{
L_12_intersect_v = CreateLineTwoPoints(P_A1_r0, P_A1_r1);
}

L_A1_lv = CreateLineTwoPoints(P_A1_l0, P_A1_l1);

//Pontoon 2
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_2_l0 = CreateLineTwoPoints(P_12intersect, P_Yintersect2);
}
else
{
L_2_l0 = CreateLineTwoPoints(P_A2_l0, P_B2_l0);
}

L_2_l1 = L_2_l0.CopyTranslate(Vector3d(0,0,Pon_height));
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_2_lY0 = CreateLineTwoPoints(P_Yintersect2, P_BYintersect);
L_2_lY1 = L_2_lY0.CopyTranslate(Vector3d(0,0,Pon_height));

}
else
{
L_2_low = CreateCircularArcFromThreePoints(P_B2_l0, P_B2_c0, P_B2_r0);
L_2_high = L_2_low.CopyTranslate(Vector3d(0,0,Pon_height));
}

L_2B_0 = CreateCircularArcFromThreePoints(P_BYintersect, P_B2_c0, P_B2_r0);

L_2B_1 = L_2B_0.CopyTranslate(Vector3d(0,0,Pon_height));
L_2_R0 = CreateLineTwoPoints(P_B2_r0, P_A2_r0);
L_2_R1 = L_2_R0.CopyTranslate(Vector3d(0,0,Pon_height));
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_2A_0 = CreateCircularArcFromThreePoints(P_A2_r0, P_A2_c0, P_A_12intersect);
}
else
{
L_2A_0 = CreateCircularArcFromThreePoints(P_A2_r0, P_A2_c0, P_A2_l0);
}

L_2A_1 = L_2A_0.CopyTranslate(Vector3d(0,0,Pon_height));
L_A2_rv = CreateLineTwoPoints(P_A2_r0, P_A2_r1);
L_B2_rv = CreateLineTwoPoints(P_B2_r0, P_B2_r1);

if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_By_intersect_v = CreateLineTwoPoints(P_BYintersect, P_BYintersect_1);
L_Y2_intersect_v = CreateLineTwoPoints(P_Yintersect2, P_Yintersect2_1);
}
else
{
L_Y2_intersect_v = CreateLineTwoPoints(P_B2_l0, P_B2_l1);
L_A2_lv = CreateLineTwoPoints(P_A2_l0, P_A2_l1);
}
// Column A

L_ColA_bottom = CreateCircularArcFromThreePoints(P_A1_l0, Point(P_A_Cx, P_A_Cy+Col_diameter/2,-Draft),P_A2_r0);  // Remaining of guidecurve
if(Pon_width/2<Math.sin(30*Math.PI/180)*Col_diameter/2)
{
L_ColA_between = CreateCircularArcFromThreePoints(P_A1_r0, P_A_12intersect, P_A2_l0);

}
// Column B

L_ColB_Yaxis = CreateLineTwoPoints(P_BYintersect, Point(P_B_Cx+Col_diameter/2, P_B_Cy, -Draft));
L_ColB_bottom = CreateCircularArcFromThreePoints(Point(P_B_Cx+Col_diameter/2, P_B_Cy, -Draft), Point(P_B_Cx, P_B_Cy+Col_diameter/2, -Draft), P_B2_r0);

//
// Creating plates
//

// Pontoon 1

Pl_1_l = CoverCurves(L_1_l1,L_1_l0,L_C1_lv,L_A1_lv);
Pl_1_l.flipNormal();
if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
Pl_1_t = CoverCurves(L_1_l1,L_1Y_1,L_1_r1,L_1A_1,L_12_intersect1);
Pl_1_r = CoverCurves(L_1_r1,L_1_r0,L_C1_rv,L_12_intersect_v);
Pl_1_b = CoverCurves(L_12_intersect0,L_1_r0,L_1_l0,L_1A_0,L_1Y_0);
}
else
{
Pl_1_t = CoverCurves(L_1_l1,L_1Y_1,L_1_r1,L_1A_1);
Pl_1_r = CoverCurves(L_1_r1,L_1_r0,L_C1_rv,L_12_intersect_v);
Pl_1_b = CoverCurves(L_1_r0,L_1_l0,L_1A_0,L_1Y_0);
}

Pl_1_b.flipNormal();

Pl_1_l.front.WetSurface = outer_wet;
Pl_1_t.front.WetSurface = outer_wet;
Pl_1_r.front.WetSurface = outer_wet;
Pl_1_b.front.WetSurface = outer_wet;

wetted_surface.add(Pl_1_l);
wetted_surface.add(Pl_1_t);
wetted_surface.add(Pl_1_r);
wetted_surface.add(Pl_1_b);

// Pontoon 2

if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
Pl_2_l = CoverCurves(L_2_l0,L_12_intersect_v,L_2_l1,L_Y2_intersect_v);
Pl_2_t = CoverCurves(L_2_R1,L_2A_1,L_12_intersect1,L_2_l1,L_2_lY1,L_2B_1);
Pl_2_b = CoverCurves(L_2_l0,L_2B_0,L_2_lY0,L_12_intersect0,L_2_R0,L_2A_0);
}
else
{
Pl_2_l = CoverCurves(L_2_l0,L_2_l1,L_Y2_intersect_v,L_A2_lv);
Pl_2_t = CoverCurves(L_2_l1,L_2_R1,L_2A_1,L_2_high);
Pl_2_b = CoverCurves(L_2_l0,L_2_low,L_2_R0,L_2A_0);
Pl_2_b.flipNormal();
}
Pl_2_l.flipNormal();

Pl_2_r = CoverCurves(L_2_R1,L_2_R0,L_A2_rv,L_B2_rv);
Pl_2_r.flipNormal();



Pl_2_l.front.WetSurface = outer_wet;
Pl_2_t.front.WetSurface = outer_wet;
Pl_2_r.front.WetSurface = outer_wet;
Pl_2_b.front.WetSurface = outer_wet;

wetted_surface.add(Pl_2_l);
wetted_surface.add(Pl_2_t);
wetted_surface.add(Pl_2_r);
wetted_surface.add(Pl_2_b);

// Column A

Pl_colA_0p = CreateShellCircularConeCylinder(P_A_c0, Col_diameter/2, P_A_cp, Col_diameter/2, phi-30, 270-phi);
Pl_colA_p1 = CreateShellCircularConeCylinder(P_A_cp, Col_diameter/2,P_A_c1, Col_diameter/2, 0, 360);




if(Pon_width/2>Math.sin(30*Math.PI/180)*Col_diameter/2)
{
Pl_colA_bottom = CoverCurves(L_2A_0,L_1A_0,L_ColA_bottom);
}
else
{
Pl_colA_bottom = CoverCurves(L_1A_0,L_2A_0,L_ColA_bottom,L_ColA_between);
}

Pl_colA_bottom.flipNormal();

Pl_colA_0p.front.WetSurface = outer_wet;
Pl_colA_p1.front.WetSurface = outer_wet;
Pl_colA_bottom.front.WetSurface = outer_wet;


wetted_surface.add(Pl_colA_0p);
wetted_surface.add(Pl_colA_p1);
wetted_surface.add(Pl_colA_bottom);



// Column B

Pl_colB_0p = CreateShellCircularConeCylinder(P_B_c0, Col_diameter/2, P_B_cp, Col_diameter/2, 0, 150-phi);
Pl_colB_p1 = CreateShellCircularConeCylinder(P_B_cp, Col_diameter/2, P_B_c1, Col_diameter/2, 0, 180);
Pl_colB_bottom = CoverCurves(L_ColB_bottom,L_ColB_Yaxis,L_2B_0);





Pl_colB_0p.front.WetSurface = outer_wet;
Pl_colB_p1.front.WetSurface = outer_wet;
Pl_colB_bottom.front.WetSurface = outer_wet;



wetted_surface.add(Pl_colB_0p);
wetted_surface.add(Pl_colB_p1);
wetted_surface.add(Pl_colB_bottom);


// Feature edges
Curve13 = CreateLineTwoPoints(L_1_r1.end1, L_1_l1.project(L_1_r1.end1));
Curve11 = CreateLineTwoPoints(Curve13.end2, L_1_l0.project(Curve13.end2));
Curve12 = CreateLineTwoPoints(Curve11.end2, L_1_r0.project(Curve11.end2));
Curve17 = CreateLineTwoPoints(L_1_r1.end1, L_2_R1.project(L_1_r1.end1));
Curve14 = CreateLineTwoPoints(Curve17.end2, L_2_R0.project(Curve17.end2));
Curve15 = CreateLineTwoPoints(Curve14.end2, L_1_r0.project(Curve14.end2));
Curve16 = CreateLineTwoPoints(L_2_l1.end2, L_2_R1.project(L_2_l1.end2));
Curve18 = CreateLineTwoPoints(Curve16.end2, L_2_R0.project(Curve16.end2));
Curve20 = CreateLineTwoPoints(Curve18.end2, L_2_l0.project(Curve18.end2));
FEdge1 = FeatureEdge(L_Y2_intersect_v);
FEdge2 = FeatureEdge(Curve13);
FEdge3 = FeatureEdge(Curve11);
FEdge4 = FeatureEdge(Curve12);
FEdge5 = FeatureEdge(Curve17);
FEdge6 = FeatureEdge(Curve14);
FEdge7 = FeatureEdge(Curve15);
FEdge8 = FeatureEdge(Curve16);
FEdge9 = FeatureEdge(Curve18);
FEdge10 = FeatureEdge(Curve20);

DefaultName(typePoint,"Mesh_Point",1,"");
DefaultName(typeGuideCurve,"Mesh_Curve",1,"");
DefaultName(typeFeatureEdge,"Mesh_FEdge",1,"");
Mesh_Curve1 = L_1A_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.25*Mesh_length));
Mesh_Curve2 = L_2B_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.25*Mesh_length));
Mesh_Curve3 = L_2A_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.25*Mesh_length));
Mesh_Curve4 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-0.25*Mesh_length));
Mesh_Curve5 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-0.25*Mesh_length));
Mesh_Curve6 = L_1A_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.75*Mesh_length));
Mesh_Curve7 = L_2B_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.75*Mesh_length));
Mesh_Curve8 = L_2A_0.copyTranslate(Vector3d(0 m,0 m,Draft-0.75*Mesh_length));
Mesh_Curve9 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-0.75*Mesh_length));
Mesh_Curve10 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-0.75*Mesh_length));
Mesh_Curve11 = L_1A_0.copyTranslate(Vector3d(0 m,0 m,Draft-1.5*Mesh_length));
Mesh_Curve12 = L_2B_0.copyTranslate(Vector3d(0 m,0 m,Draft-1.5*Mesh_length));
Mesh_Curve13 = L_2A_0.copyTranslate(Vector3d(0 m,0 m,Draft-1.5*Mesh_length));
Mesh_Curve14 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-1.5*Mesh_length));
Mesh_Curve15 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,Draft-1.5*Mesh_length));
Mesh_Curve16 = L_1A_0.copyTranslate(Vector3d(0 m,0 m,Draft));
Mesh_Curve17 = L_2B_0.copyTranslate(Vector3d(0 m,0 m,Draft));
Mesh_Curve18 = L_2A_0.copyTranslate(Vector3d(0 m,0 m,Draft));
Mesh_Curve19 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,Draft));
Mesh_Curve20 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,Draft));
Mesh_Curve21 = L_1_l0.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve22 = L_1_r0.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve23 = L_2_l0.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve24 = L_2_R0.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve25 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve26 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve27 = L_1_l1.copyTranslate(Vector3d(0 m,0 m,-0.25*Mesh_length));
Mesh_Curve28 = L_1_r1.copyTranslate(Vector3d(0 m,0 m,-0.25*Mesh_length));
Mesh_Curve29 = L_2_l1.copyTranslate(Vector3d(0 m,0 m,-0.25*Mesh_length));
Mesh_Curve30 = L_2_R1.copyTranslate(Vector3d(0 m,0 m,-0.25*Mesh_length));
Mesh_Curve31 = L_2_l1.offset(-0.25*Mesh_length);
Mesh_Curve32 = L_2_R1.offset(-0.25*Mesh_length);
Mesh_Point1 = Point(Mesh_Curve32.intersect(L_2B_1));
Mesh_Point2 = Point(L_1_l1.end2.x+0.25*Mesh_length, L_1_l1.end2.y, L_1_l1.end2.z);
Mesh_Point3 = Point(P_A1_l1.x+0.25*Mesh_length, P_A1_l1.y, P_A1_l1.z);
Mesh_Curve33 = CreateLineTwoPoints(Mesh_Point3, Mesh_Point2);
Mesh_Point4 = Point(Mesh_Curve33.intersect(L_1A_1));
Mesh_Point5 = Point(P_C1_r1.x-0.25*Mesh_length, P_C1_r1.y, P_C1_r1.z);
Mesh_Point6 = Point(P_12intersect_1.x-0.25*Mesh_length, P_12intersect_1.y+1m, P_12intersect_1.z);
Mesh_Curve34 = CreateLineTwoPoints(Mesh_Point5, Mesh_Point6);
Mesh_Curve35 = CreateLineTwoPoints(Mesh_Curve34.intersect(L_12_intersect1), Mesh_Curve31.end1);
Mesh_Curve36 = L_1A_1.offset(0.25*Mesh_length);
Mesh_Curve37 = L_2A_1.offset(0.25*Mesh_length);
Mesh_Curve38 = L_1A_1.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve39 = L_2A_1.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve40 = L_2B_1.copyTranslate(Vector3d(0 m,0 m,0.25*Mesh_length));
Mesh_Curve41 = CreateLineTwoPoints(L_2_R1.end1, Mesh_Curve40.end2);
Mesh_Curve42 = Mesh_Curve34.divideAt(Mesh_Curve34.parameter(Mesh_Curve35.curvePoint(1)));
Delete(Mesh_Curve42);
Mesh_Curve42 = Mesh_Curve33.divideAt(Mesh_Curve33.parameter(Mesh_Point4));
Delete(Mesh_Curve33);
Mesh_Curve43 = Mesh_Curve32.divideAt(Mesh_Curve32.parameter(Mesh_Point1));
Delete(Mesh_Curve32);
Mesh_Curve44 = CreateLineTwoPoints(Mesh_Curve38.end2, P_A1_l1);
Mesh_Point7 = Point(L_1_l1.intersect(Mesh_Curve36));
Mesh_Curve45 = Mesh_Curve36.divideAt(Mesh_Curve36.parameter(Mesh_Point7));
Delete(Mesh_Curve45);
Mesh_Curve45 = CreateLineTwoPoints(Mesh_Curve39.end1, L_2_R1.end2);
Mesh_Point8 = Point(L_2_R1.intersect(Mesh_Curve37));
Mesh_Curve46 = Mesh_Curve37.divideAt(Mesh_Curve37.parameter(Mesh_Point8));
Delete(Mesh_Curve37);
tmpArrayOfCurves = JoinMultipleCurves(Array(Mesh_Curve31, Mesh_Curve35), true);
Rename(tmpArrayOfCurves[0], "Mesh_Curve47");
Mesh_Point9 = Point(L_2A_1.intersect(Mesh_Curve43));
Mesh_Curve48 = Mesh_Curve43.divideAt(Mesh_Curve43.parameter(Mesh_Point9));
Delete(Mesh_Curve48);
Mesh_Curve48 = Mesh_Curve43.copyTranslate(Vector3d(0 m,0 m,-Pon_height));
Mesh_Curve49 = Mesh_Curve34.copyTranslate(Vector3d(0 m,0 m,-Pon_height));
Mesh_Curve50 = Mesh_Curve47.copyTranslate(Vector3d(0 m,0 m,-Pon_height));
Mesh_Curve51 = Mesh_Curve42.copyTranslate(Vector3d(0 m,0 m,-Pon_height));
Mesh_Curve52 = CreateLineTwoPoints(L_A1_lv.end1, P_A_c0);
Mesh_Curve53 = CreateLineTwoPoints(L_ColA_bottom.end2, Mesh_Curve52.end2);
Mesh_Curve54 = CreateLineTwoPoints(P_A_Cir1, Mesh_Curve53.end2);
Mesh_Curve55 = Mesh_Curve53.divideAt(0.25*Mesh_length/(Col_diameter/2));
Mesh_Curve56 = Mesh_Curve54.divideAt(0.25*Mesh_length/(Col_diameter/2));
Mesh_Curve57 = Mesh_Curve52.divideAt(0.25*Mesh_length/(Col_diameter/2));
Mesh_Curve58 = CreateCircularArcFromThreePoints(Mesh_Curve57.end1, Mesh_Curve56.end1, Mesh_Curve53.end2);
Mesh_Curve59 = CreateLineTwoPoints(L_B2_rv.end1, P_B_c0);
Mesh_Curve60 = CreateLineTwoPoints(L_ColB_bottom.end1, Mesh_Curve59.end2);
Mesh_Curve61 = Mesh_Curve59.divideAt(0.25*Mesh_length/(Col_diameter/2));
Mesh_Curve62 = Mesh_Curve60.divideAt(0.25*Mesh_length/(Col_diameter/2));
Mesh_Curve63 = Mesh_Curve59.copyRotate(P_B_c0,Vector3d(0, 0, 1),-45);
Mesh_Curve64 = CreateCircularArcFromThreePoints(Mesh_Curve61.end1, Mesh_Curve63.end2, Mesh_Curve60.end2);
Mesh_Curve65 = CreateLineTwoPoints(Mesh_Curve59.end2, Mesh_Curve48.end1);
Mesh_Curve66 = CreateLineTwoPoints(Mesh_Curve57.end1, Mesh_Curve51.end1);
Mesh_Curve67 = CreateLineTwoPoints(Mesh_Curve53.end2, P_A_Cir4);
Mesh_Curve68 = L_2B_1.offset(0.25*Mesh_length);
Mesh_Point10 = Point(L_2_R1.intersect(Mesh_Curve68));
Mesh_Curve69 = Mesh_Curve68.divideAt(Mesh_Curve68.parameter(Mesh_Point10));
Delete(Mesh_Curve69);
Mesh_FEdge1 = FeatureEdge(L_2_R1);
Mesh_FEdge2 = FeatureEdge(Mesh_Curve1);
Mesh_FEdge3 = FeatureEdge(Mesh_Curve2);
Mesh_FEdge4 = FeatureEdge(Mesh_Curve3);
Mesh_FEdge5 = FeatureEdge(Mesh_Curve4);
Mesh_FEdge6 = FeatureEdge(Mesh_Curve5);
Mesh_FEdge7 = FeatureEdge(Mesh_Curve6);
Mesh_FEdge8 = FeatureEdge(Mesh_Curve7);
Mesh_FEdge9 = FeatureEdge(Mesh_Curve8);
Mesh_FEdge10 = FeatureEdge(Mesh_Curve9);
Mesh_FEdge11 = FeatureEdge(Mesh_Curve10);
Mesh_FEdge12 = FeatureEdge(Mesh_Curve11);
Mesh_FEdge13 = FeatureEdge(Mesh_Curve12);
Mesh_FEdge14 = FeatureEdge(Mesh_Curve13);
Mesh_FEdge15 = FeatureEdge(Mesh_Curve14);
Mesh_FEdge16 = FeatureEdge(Mesh_Curve15);
Mesh_FEdge17 = FeatureEdge(Mesh_Curve21);
Mesh_FEdge18 = FeatureEdge(Mesh_Curve22);
Mesh_FEdge19 = FeatureEdge(Mesh_Curve23);
Mesh_FEdge20 = FeatureEdge(Mesh_Curve24);
Mesh_FEdge21 = FeatureEdge(Mesh_Curve25);
Mesh_FEdge22 = FeatureEdge(Mesh_Curve26);
Mesh_FEdge23 = FeatureEdge(Mesh_Curve27);
Mesh_FEdge24 = FeatureEdge(Mesh_Curve28);
Mesh_FEdge25 = FeatureEdge(Mesh_Curve29);
Mesh_FEdge26 = FeatureEdge(Mesh_Curve30);
Mesh_FEdge27 = FeatureEdge(Mesh_Curve49);
Mesh_FEdge28 = FeatureEdge(Mesh_Curve44);
Mesh_FEdge29 = FeatureEdge(Mesh_Curve43);
Mesh_FEdge30 = FeatureEdge(Mesh_Curve34);
Mesh_FEdge31 = FeatureEdge(Mesh_Curve48);
Mesh_FEdge32 = FeatureEdge(Mesh_Curve47);
Mesh_FEdge33 = FeatureEdge(Mesh_Curve36);
Mesh_FEdge34 = FeatureEdge(Mesh_Curve38);
Mesh_FEdge35 = FeatureEdge(Mesh_Curve39);
Mesh_FEdge36 = FeatureEdge(Mesh_Curve40);
Mesh_FEdge37 = FeatureEdge(Mesh_Curve41);
Mesh_FEdge38 = FeatureEdge(Mesh_Curve42);
Mesh_FEdge39 = FeatureEdge(Mesh_Curve45);
Mesh_FEdge40 = FeatureEdge(Mesh_Curve46);
Mesh_FEdge41 = FeatureEdge(Mesh_Curve50);
Mesh_FEdge42 = FeatureEdge(Mesh_Curve51);
Mesh_FEdge43 = FeatureEdge(Mesh_Curve58);
Mesh_FEdge44 = FeatureEdge(Mesh_Curve64);
Mesh_FEdge45 = FeatureEdge(Mesh_Curve65);
Mesh_FEdge46 = FeatureEdge(Mesh_Curve66);
Mesh_FEdge47 = FeatureEdge(Mesh_Curve67);
Mesh_FEdge48 = FeatureEdge(Mesh_Curve68);


LC1.generateAppliedLoads();
GenieRules.Meshing.faceMeshStrategy = PatchSurfQuadMesher;
Md_def.setDefault();


LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

GenieRules.Meshing.superElementType = 1;
Meshing_panels.execute();

ExportMeshFem().DoExport("T1.FEM");


/*
Modeling of the structure above MWL, still applying symmetry
 */



Pl_colA_12 = CreateShellCircularConeCylinder(P_A_c1, Col_diameter/2, P_A_c2, Col_diameter/2, 0, 360);


Pl_colB_12 = CreateShellCircularConeCylinder(P_B_c1, Col_diameter/2, P_B_c2, Col_diameter/2, 0, 180);

phi2 = Math.asin((Deck_width/2)/(Col_diameter/2))*180/Math.PI;

P_A1_cx_top = P_A_Cx + (Col_diameter/2*Math.cos(-90*Math.PI/180));
P_A1_cy_top = P_A_Cy + (Col_diameter/2*Math.sin(-90*Math.PI/180));
P_A1_c0_top = Point(P_A1_cx, P_A1_cy, -Draft+Col_height-Deck_height);
P_A1_c1_top = P_A1_c0_top.CopyTranslate(Vector3d(0, 0, Deck_height));

P_A1_lx_top = P_A_Cx + (Col_diameter/2*Math.cos((-90-phi2)*Math.PI/180));
P_A1_ly_top = P_A_Cy + (Col_diameter/2*Math.sin((-90-phi2)*Math.PI/180));
P_A1_l0_top = Point(P_A1_lx_top, P_A1_ly_top, -Draft+Col_height-Deck_height);
P_A1_l1_top = P_A1_l0_top.CopyTranslate(Vector3d(0, 0, Deck_height));

P_A1_rx_top = P_A_Cx + (Col_diameter/2*Math.cos((-90+phi2)*Math.PI/180));
P_A1_ry_top = P_A_Cy + (Col_diameter/2*Math.sin((-90+phi2)*Math.PI/180));
P_A1_r0_top = Point(P_A1_rx_top, P_A1_ry_top, -Draft+Col_height-Deck_height);
P_A1_r1_top = P_A1_r0_top.CopyTranslate(Vector3d(0, 0, Deck_height));


P_C1_cx_top = P_A_Cx;
P_C1_cy_top = 0m;
P_C1_c0_top = Point(P_C1_cx_top, P_C1_cy_top, -Draft+Col_height-Deck_height);
P_C1_c1_top = P_C1_c0_top.CopyTranslate(Vector3d(0, 0, Deck_height));

P_C1_lx_top = P_A_Cx-Deck_width/2;
P_C1_ly_top = 0m;
P_C1_l0_top = Point(P_C1_lx_top, P_C1_ly_top, -Draft+Col_height-Deck_height);
P_C1_l1_top = P_C1_l0_top.CopyTranslate(Vector3d(0, 0, Deck_height));

P_C1_rx_top = P_A_Cx+Deck_width/2;
P_C1_ry_top = 0m;
P_C1_r0_top = Point(P_C1_rx_top, P_C1_ry_top, -Draft+Col_height-Deck_width);
P_C1_r1_top = P_C1_r0_top.CopyTranslate(Vector3d(0, 0, Deck_height));

Point32 = P_A1_r1_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);
Point33 = P_A1_r0_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);
Point34 = P_A1_l1_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);
Point35 = P_A1_l0_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);
Point36 = P_A1_c1_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);
Point37 = P_A1_c0_top.copyRotate(P_A_c2,Vector3d(0, 0, 1),60);

Point38 = Point33.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
Point39 = Point32.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
Point40 = Point36.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
Point41 = Point37.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
Point42 = Point35.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
Point43 = Point34.copyRotate(P_A_c2,Vector3d(0, 0, 1),180);
autoMSet = Set();
autoMSet.add(Point38);
autoMSet.add(Point39);
autoMSet.add(Point40);
autoMSet.add(Point41);
autoMSet.add(Point42);
autoMSet.add(Point43);
autoMSet.moveTranslate(Vector3d(P_A_c2, P_B_c2),geUNCONNECTED);
Delete(autoMSet);

Curve15 = CreateCircularArcFromThreePoints(P_A1_l1_top, P_A1_c1_top, P_A1_r1_top);
Curve16 = CreateCircularArcFromThreePoints(P_A1_l0_top, P_A1_c0_top, P_A1_r0_top);
Curve17 = CreateCircularArcFromThreePoints(Point34, Point36, Point32);
Curve18 = CreateCircularArcFromThreePoints(Point35, Point37, Point33);
Curve19 = CreateCircularArcFromThreePoints(Point39, Point40, Point43);
Curve20 = CreateCircularArcFromThreePoints(Point38, Point41, Point42);
Curve21 = CreateLineTwoPoints(Curve15.end1, Curve16.end1);
Curve22 = CreateLineTwoPoints(Curve15.end2, Curve16.end2);
Curve23 = CreateLineTwoPoints(Curve17.end1, Curve18.end1);
Curve24 = CreateLineTwoPoints(Curve17.end2, Curve18.end2);
Curve25 = CreateLineTwoPoints(P_C1_l1_top, P_C1_l0_top);
Curve26 = CreateLineTwoPoints(P_C1_r0_top, P_C1_r1_top);
Curve27 = CreateLineTwoPoints(Curve26.end2, Curve25.end1);
Curve28 = CreateLineTwoPoints(Curve26.end1, Curve25.end2);
Curve29 = CreateLineTwoPoints(Curve20.end1, Curve19.end1);
Curve30 = CreateLineTwoPoints(Curve20.end2, Curve19.end2);
Curve31 = CreateLineTwoPoints(Curve19.end2, Curve17.end2);
Curve32 = CreateLineTwoPoints(Curve19.end1, Curve17.end1);
Curve33 = CreateLineTwoPoints(Curve18.end1, Curve20.end1);
Curve34 = CreateLineTwoPoints(Curve18.end2, Curve20.end2);
Curve35 = CreateLineTwoPoints(Curve25.end1, Curve21.end1);
Curve36 = CreateLineTwoPoints(Curve15.end2, Curve26.end2);
Curve37 = CreateLineTwoPoints(Curve16.end2, Curve26.end1);
Curve38 = CreateLineTwoPoints(Curve25.end2, Curve16.end1);
Pl21 = CoverCurves(Curve21,Curve25,Curve35,Curve38);
Pl21.flipNormal();
Pl22 = CoverCurves(Curve15,Curve27,Curve35,Curve36);
Pl22.flipNormal();
Pl23 = CoverCurves(Curve22,Curve26,Curve36,Curve37);
Pl23.flipNormal();
Pl24 = CoverCurves(Curve16,Curve28,Curve37,Curve38);
Pl25 = CoverCurves(Curve17,Curve19,Curve31,Curve32);
Pl25.flipNormal();
Pl26 = CoverCurves(Curve24,Curve30,Curve31,Curve34);
Pl27 = CoverCurves(Curve18,Curve20,Curve33,Curve34);
Pl28 = CoverCurves(Curve23,Curve29,Curve32,Curve33);
Pl27.flipNormal();
Pl26.flipNormal();

Curve39 = L_1A_0.copyTranslate(Vector3d(0 m,0 m,Col_height));
Curve40 = L_2B_0.copyTranslate(Vector3d(0 m,0 m,Col_height));
Curve41 = L_2A_0.copyTranslate(Vector3d(0 m,0 m,Col_height));
Curve42 = L_ColA_bottom.copyTranslate(Vector3d(0 m,0 m,Col_height));
Curve43 = L_ColB_Yaxis.copyTranslate(Vector3d(0 m,0 m,Col_height));
Curve44 = L_ColB_bottom.copyTranslate(Vector3d(0 m,0 m,Col_height));


Pl29 = CoverCurves(Curve40,Curve43,Curve44);
Pl30 = CoverCurves(Curve39,Curve41,Curve42);
Pl30.flipNormal();

wetted_surface.add(Pl_colA_12);

wetted_surface.add(Pl_colB_12);

wetted_surface.front.wetSurface = outer_wet;


GenieRules.Meshing.superElementType = 2;
Meshing_panels.execute();
ExportMeshFem().DoExport("T2.FEM");


/*
Internal structure
 */

damage_lower_extent = 3m;
damage_upper_extent = 5m;
horizontal_extent_sector = 0.125; // fraction of circumference
penetration = 1.5m;

Tck2.setDefault();

Pl_Y_intersect = CoverCurves(L_2_lY0,L_2_lY1,L_By_intersect_v,L_Y2_intersect_v);
Pl_12_intersect = CoverCurves(L_12_intersect0,L_12_intersect1,L_12_intersect_v,L_A12_intersect_v);
L_1_top_c = CreateLineTwoPoints(P_C1_c1, P_A1_c1);
L_1_colA_v = CreateLineTwoPoints(P_A1_c1, P_A1_c0);
L_1_bot_c = CreateLineTwoPoints(L_1_colA_v.end2, P_C1_c0);
L_1_Y_v = CreateLineTwoPoints(L_1_top_c.end1, L_1_bot_c.end2);
Pl_1_internal = CoverCurves(L_1_top_c,L_1_colA_v,L_1_bot_c,L_1_Y_v);
Pl_1Y = CoverCurves(L_1Y_0,L_1Y_1,L_C1_lv,L_C1_rv);
Pl_2mid = Pl_1Y.copyRotate(P_A_c0, Vector3d(0, 0, 1), 60);
Pl261 = CoverCurves(L_1A_0,L_1A_1,L_A12_intersect_v,L_A1_lv);
Pl271 = CoverCurves(L_A12_intersect_v,L_2A_0,L_2A_1,L_A2_rv);
Pl281 = CoverCurves(L_2B_0,L_2B_1,L_B2_rv,L_By_intersect_v);
L_2_ColA_v = CreateLineTwoPoints(P_A2_c0, P_A2_c1);
L_2_top_c = CreateLineTwoPoints(L_2_ColA_v.end2, P_B2_c1);
L_2_ColB_v = CreateLineTwoPoints(L_2_top_c.end2, P_B2_c0);
L_2_bot_c = CreateLineTwoPoints(L_2_ColB_v.end2, L_2_ColA_v.end1);
Pl_2_internal = CoverCurves(L_2_ColA_v,L_2_top_c,L_2_ColB_v,L_2_bot_c);
Point26 = P_A_c1.copyTranslate(Vector3d(0 m,0 m, -damage_lower_extent));
Point27 = P_A_c1.copyTranslate(Vector3d(0 m,0 m, damage_upper_extent));
Point28 = Point26.copyTranslate(Vector3d(0 m,Col_diameter/2,0 m));
Point29 = Point27.copyTranslate(Vector3d(0 m,Col_diameter/2,0 m));
Point30 = Point26.copyTranslate(Vector3d(0 m,Col_diameter/2-penetration,0 m));
Point31 = Point27.copyTranslate(Vector3d(0 m,Col_diameter/2-penetration,0 m));
Pl462 = Plate(Point28,Point30,Point31,Point29);
Point32 = Point28.copyRotate(P_A_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point33 = Point29.copyRotate(P_A_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point34 = Point30.copyRotate(P_A_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point35 = Point31.copyRotate(P_A_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl462, "Pl472");
ModelTransformer(MyModelTransformerMap).copyRotate(P_A_c1, Vector3d(0, 0, 1), 360*horizontal_extent_sector, Math.Round(1/horizontal_extent_sector-1));
Curve15 = CreateCircleFromCenterAndRadius(Point26, Point30, Point34);
Curve16 = CreateCircleFromCenterAndRadius(Point26, Point28, Point32);
Curve19 = CreateLineTwoPoints(Curve15.end1, Point31);
Pl541 = SweepCurve(Curve15, Curve19);
Pl551 = CoverCurves(Curve15,Curve16);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl551, "Pl561");
ModelTransformer(MyModelTransformerMap).copyTranslate(Vector3d(0 m,0 m, damage_lower_extent+damage_upper_extent));
Point36 = Point26.copyTranslate(Vector3d(P_A_c1, P_B_c1));
Point37 = Point27.copyTranslate(Vector3d(P_A_c1, P_B_c1));
Point38 = Point36.copyTranslate(Vector3d(Col_diameter/2,0m, 0 m));
Point39 = Point37.copyTranslate(Vector3d(Col_diameter/2, 0m, 0 m));
Point40 = Point36.copyTranslate(Vector3d(Col_diameter/2-penetration, 0m, 0 m));
Point41 = Point37.copyTranslate(Vector3d(Col_diameter/2-penetration, 0m, 0 m));
Point42 = Point38.copyRotate(P_B_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point43 = Point39.copyRotate(P_B_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point44 = Point40.copyRotate(P_B_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point45 = Point41.copyRotate(P_B_c1,Vector3d(0, 0, 1),360*horizontal_extent_sector);
Point46 = Point36.copyTranslate(Vector3d(-Col_diameter/2,0m, 0 m));
Point47 = Point37.copyTranslate(Vector3d(-Col_diameter/2, 0m, 0 m));
Point48 = Point36.copyTranslate(Vector3d(-Col_diameter/2+penetration, 0m, 0 m));
Point49 = Point37.copyTranslate(Vector3d(-Col_diameter/2+penetration, 0m, 0 m));
Pl571 = Plate(Point39,Point38,Point40,Point41);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl571, "Pl581");
ModelTransformer(MyModelTransformerMap).copyRotate(P_B_c1, Vector3d(0, 0, 1), 360*horizontal_extent_sector, Math.Round(0.5/horizontal_extent_sector));
Curve18 = CreateCircularArcFromThreePoints(Point40, Point44, Point48);
Curve17 = CreateCircularArcFromThreePoints(Point38, Point42, Point46);
Curve20 = CreateLineTwoPoints(Curve18.end1, Point41);
Pl621 = SweepCurve(Curve18, Curve20);
Curve21 = CreateLineTwoPoints(Curve17.end1, Curve20.end1);
Curve22 = CreateLineTwoPoints(Curve18.end2, Curve17.end2);
Pl631 = CoverCurves(Curve18,Curve17,Curve21,Curve22);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl631, "Pl641");
ModelTransformer(MyModelTransformerMap).copyTranslate(Vector3d(0 m,0 m,damage_lower_extent+damage_upper_extent));
Pl4600 = CoverCurves(Curve15);
Curve45 = CreateLineTwoPoints(Point(39.42876129 m,0 m,-3 m), Point(49.74876129 m,0 m,-3 m));
Pl4700 = CoverCurves(Curve18,Curve45);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl4600, "Pl4800");
MyModelTransformerMap.Add(Pl4700, "Pl4900");
ModelTransformer(MyModelTransformerMap).copyTranslate(Vector3d(0 m,0 m,damage_lower_extent+damage_upper_extent));

/*
Mirror model
 */

MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(FEdge1, "Mesh_FEdge49");
MyModelTransformerMap.Add(FEdge2, "Mesh_FEdge50");
MyModelTransformerMap.Add(FEdge3, "Mesh_FEdge51");
MyModelTransformerMap.Add(FEdge4, "Mesh_FEdge52");
MyModelTransformerMap.Add(FEdge5, "Mesh_FEdge53");
MyModelTransformerMap.Add(FEdge6, "Mesh_FEdge54");
MyModelTransformerMap.Add(FEdge7, "Mesh_FEdge55");
MyModelTransformerMap.Add(FEdge8, "Mesh_FEdge56");
MyModelTransformerMap.Add(FEdge9, "Mesh_FEdge57");
MyModelTransformerMap.Add(FEdge10, "Mesh_FEdge58");
MyModelTransformerMap.Add(Mesh_FEdge1, "Mesh_FEdge59");
MyModelTransformerMap.Add(Mesh_FEdge2, "Mesh_FEdge60");
MyModelTransformerMap.Add(Mesh_FEdge3, "Mesh_FEdge61");
MyModelTransformerMap.Add(Mesh_FEdge4, "Mesh_FEdge62");
MyModelTransformerMap.Add(Mesh_FEdge5, "Mesh_FEdge63");
MyModelTransformerMap.Add(Mesh_FEdge6, "Mesh_FEdge64");
MyModelTransformerMap.Add(Mesh_FEdge7, "Mesh_FEdge65");
MyModelTransformerMap.Add(Mesh_FEdge8, "Mesh_FEdge66");
MyModelTransformerMap.Add(Mesh_FEdge9, "Mesh_FEdge67");
MyModelTransformerMap.Add(Mesh_FEdge10, "Mesh_FEdge68");
MyModelTransformerMap.Add(Mesh_FEdge11, "Mesh_FEdge69");
MyModelTransformerMap.Add(Mesh_FEdge12, "Mesh_FEdge70");
MyModelTransformerMap.Add(Mesh_FEdge13, "Mesh_FEdge71");
MyModelTransformerMap.Add(Mesh_FEdge14, "Mesh_FEdge72");
MyModelTransformerMap.Add(Mesh_FEdge15, "Mesh_FEdge73");
MyModelTransformerMap.Add(Mesh_FEdge16, "Mesh_FEdge74");
MyModelTransformerMap.Add(Mesh_FEdge17, "Mesh_FEdge75");
MyModelTransformerMap.Add(Mesh_FEdge18, "Mesh_FEdge76");
MyModelTransformerMap.Add(Mesh_FEdge19, "Mesh_FEdge77");
MyModelTransformerMap.Add(Mesh_FEdge20, "Mesh_FEdge78");
MyModelTransformerMap.Add(Mesh_FEdge21, "Mesh_FEdge79");
MyModelTransformerMap.Add(Mesh_FEdge22, "Mesh_FEdge80");
MyModelTransformerMap.Add(Mesh_FEdge23, "Mesh_FEdge81");
MyModelTransformerMap.Add(Mesh_FEdge24, "Mesh_FEdge82");
MyModelTransformerMap.Add(Mesh_FEdge25, "Mesh_FEdge83");
MyModelTransformerMap.Add(Mesh_FEdge26, "Mesh_FEdge84");
MyModelTransformerMap.Add(Mesh_FEdge27, "Mesh_FEdge85");
MyModelTransformerMap.Add(Mesh_FEdge28, "Mesh_FEdge86");
MyModelTransformerMap.Add(Mesh_FEdge29, "Mesh_FEdge87");
MyModelTransformerMap.Add(Mesh_FEdge30, "Mesh_FEdge88");
MyModelTransformerMap.Add(Mesh_FEdge31, "Mesh_FEdge89");
MyModelTransformerMap.Add(Mesh_FEdge32, "Mesh_FEdge90");
MyModelTransformerMap.Add(Mesh_FEdge33, "Mesh_FEdge91");
MyModelTransformerMap.Add(Mesh_FEdge34, "Mesh_FEdge92");
MyModelTransformerMap.Add(Mesh_FEdge35, "Mesh_FEdge93");
MyModelTransformerMap.Add(Mesh_FEdge36, "Mesh_FEdge94");
MyModelTransformerMap.Add(Mesh_FEdge37, "Mesh_FEdge95");
MyModelTransformerMap.Add(Mesh_FEdge38, "Mesh_FEdge96");
MyModelTransformerMap.Add(Mesh_FEdge39, "Mesh_FEdge97");
MyModelTransformerMap.Add(Mesh_FEdge40, "Mesh_FEdge98");
MyModelTransformerMap.Add(Mesh_FEdge41, "Mesh_FEdge99");
MyModelTransformerMap.Add(Mesh_FEdge42, "Mesh_FEdge100");
MyModelTransformerMap.Add(Mesh_FEdge43, "Mesh_FEdge101");
MyModelTransformerMap.Add(Mesh_FEdge44, "Mesh_FEdge102");
MyModelTransformerMap.Add(Mesh_FEdge45, "Mesh_FEdge103");
MyModelTransformerMap.Add(Mesh_FEdge46, "Mesh_FEdge104");
MyModelTransformerMap.Add(Mesh_FEdge47, "Mesh_FEdge105");
MyModelTransformerMap.Add(Mesh_FEdge48, "Mesh_FEdge106");
MyModelTransformerMap.Add(Pl21, "Pl44");
MyModelTransformerMap.Add(Pl22, "Pl45");
MyModelTransformerMap.Add(Pl23, "Pl46");
MyModelTransformerMap.Add(Pl24, "Pl47");
MyModelTransformerMap.Add(Pl25, "Pl48");
MyModelTransformerMap.Add(Pl26, "Pl49");
MyModelTransformerMap.Add(Pl27, "Pl50");
MyModelTransformerMap.Add(Pl28, "Pl51");
MyModelTransformerMap.Add(Pl29, "Pl52");
MyModelTransformerMap.Add(Pl30, "Pl53");
MyModelTransformerMap.Add(Pl261, "Pl54");
MyModelTransformerMap.Add(Pl271, "Pl55");
MyModelTransformerMap.Add(Pl281, "Pl56");
MyModelTransformerMap.Add(Pl462, "Pl57");
MyModelTransformerMap.Add(Pl472, "Pl58");
MyModelTransformerMap.Add(Pl473, "Pl59");
MyModelTransformerMap.Add(Pl474, "Pl60");
MyModelTransformerMap.Add(Pl475, "Pl61");
MyModelTransformerMap.Add(Pl476, "Pl62");
MyModelTransformerMap.Add(Pl477, "Pl63");
MyModelTransformerMap.Add(Pl478, "Pl64");
MyModelTransformerMap.Add(Pl541, "Pl65");
MyModelTransformerMap.Add(Pl551, "Pl66");
MyModelTransformerMap.Add(Pl561, "Pl67");
MyModelTransformerMap.Add(Pl571, "Pl68");
MyModelTransformerMap.Add(Pl581, "Pl69");
MyModelTransformerMap.Add(Pl582, "Pl70");
MyModelTransformerMap.Add(Pl583, "Pl71");
MyModelTransformerMap.Add(Pl584, "Pl72");
MyModelTransformerMap.Add(Pl621, "Pl73");
MyModelTransformerMap.Add(Pl631, "Pl74");
MyModelTransformerMap.Add(Pl641, "Pl75");
MyModelTransformerMap.Add(Pl4600, "Pl76");
MyModelTransformerMap.Add(Pl4700, "Pl77");
MyModelTransformerMap.Add(Pl4800, "Pl78");
MyModelTransformerMap.Add(Pl4900, "Pl79");
MyModelTransformerMap.Add(Pl_1Y, "Pl80");
MyModelTransformerMap.Add(Pl_1_b, "Pl81");
MyModelTransformerMap.Add(Pl_1_internal, "Pl82");
MyModelTransformerMap.Add(Pl_1_l, "Pl83");
MyModelTransformerMap.Add(Pl_1_r, "Pl84");
MyModelTransformerMap.Add(Pl_1_t, "Pl85");
MyModelTransformerMap.Add(Pl_2mid, "Pl86");
MyModelTransformerMap.Add(Pl_2_b, "Pl87");
MyModelTransformerMap.Add(Pl_2_internal, "Pl88");
MyModelTransformerMap.Add(Pl_2_l, "Pl89");
MyModelTransformerMap.Add(Pl_2_r, "Pl90");
MyModelTransformerMap.Add(Pl_2_t, "Pl91");
MyModelTransformerMap.Add(Pl_12_intersect, "Pl92");
MyModelTransformerMap.Add(Pl_colA_0p, "Pl93");
MyModelTransformerMap.Add(Pl_colA_12, "Pl94");
MyModelTransformerMap.Add(Pl_colA_bottom, "Pl95");
MyModelTransformerMap.Add(Pl_colA_p1, "Pl96");
MyModelTransformerMap.Add(Pl_colB_0p, "Pl97");
MyModelTransformerMap.Add(Pl_colB_12, "Pl98");
MyModelTransformerMap.Add(Pl_colB_bottom, "Pl99");
MyModelTransformerMap.Add(Pl_colB_p1, "Pl100");
MyModelTransformerMap.Add(Pl_Y_intersect, "Pl101");
ModelTransformer(MyModelTransformerMap).copyMirror(P_C1_c0, Vector3d(0, 1, 0));
Mesh_Point53 = P_A_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point54 = P_A_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point55 = P_A_c2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point56 = P_A_cp.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point57 = P_A_F1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point58 = P_A_F2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point59 = P_A_F3.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point60 = P_A_Cir1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point61 = P_A_Cir2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point62 = P_A_Cir3.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point63 = P_A_Cir4.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point64 = P_A_Cir5.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point65 = P_B_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point66 = P_B_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point67 = P_B_c2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point68 = P_B_cp.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point69 = P_A1_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point70 = P_A1_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point71 = P_A1_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point72 = P_A1_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point73 = P_A1_r0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point74 = P_A1_r1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point75 = P_C1_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point76 = P_C1_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point77 = P_C1_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point78 = P_C1_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point79 = P_C1_r0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point80 = P_C1_r1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point81 = P_A2_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point82 = P_A2_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point83 = P_A2_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point84 = P_A2_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point85 = P_A2_r0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point86 = P_A2_r1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point87 = P_B2_c0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point88 = P_B2_c1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point89 = P_B2_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point90 = P_B2_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point91 = P_B2_r0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point92 = P_B2_r1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point93 = P_Yintersect2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point94 = P_Yintersect2_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point95 = P_BYintersect.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point96 = P_BYintersect_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point97 = P_12intersect.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point98 = P_12intersect_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point99 = P_A_12intersect.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point100 = P_A_12intersect_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point101 = P_eliptic_2Ar.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point102 = P_eliptic_2Br.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point103 = P_eliptic_1Al.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point104 = Point37.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point105 = Point36.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point106 = Point38.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point107 = Point39.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point108 = Point40.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point109 = Point41.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point110 = Point42.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point111 = Point43.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point112 = Point26.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point113 = Point27.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point114 = Point28.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point115 = Point29.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point116 = Point30.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point117 = Point31.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point118 = Point32.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point119 = Point33.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point120 = Point34.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point121 = Point35.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point122 = Mesh_Point1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point123 = P_A1_r1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point124 = P_C1_c0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point125 = P_C1_c1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point126 = Mesh_Point2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point127 = Mesh_Point3.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point128 = P_C1_l0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point129 = Mesh_Point4.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point130 = Mesh_Point5.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point131 = Mesh_Point6.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point132 = P_C1_l1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point133 = P_C1_r0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point134 = P_C1_r1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve75 = L_1_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve76 = L_1_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve77 = L_1_r0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve78 = L_1_r1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve79 = L_12_intersect0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve80 = L_12_intersect1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve81 = L_1A_0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve82 = L_1A_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve83 = L_1Y_0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve84 = L_1Y_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve85 = L_C1_lv.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve86 = L_C1_rv.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve87 = L_12_intersect_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve88 = L_A12_intersect_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve89 = L_A1_lv.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve90 = L_2_l0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve91 = L_2_l1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve92 = L_2_lY0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve93 = L_2_lY1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve94 = L_2B_0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve95 = L_2B_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve96 = L_2_R0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve97 = L_2_R1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve98 = L_2A_0.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve99 = L_2A_1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve100 = L_A2_rv.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve101 = L_B2_rv.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve102 = L_By_intersect_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve103 = L_Y2_intersect_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve104 = L_ColA_bottom.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve105 = L_ColB_Yaxis.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve106 = L_ColB_bottom.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve107 = Curve13.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve108 = Curve11.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve109 = Curve12.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve110 = Curve17.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve111 = Curve14.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve112 = Curve19.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve113 = Curve20.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve114 = Curve18.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve115 = Curve22.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve116 = Mesh_Curve1.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve117 = Mesh_Curve2.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve118 = Mesh_Curve3.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve119 = Mesh_Curve4.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve120 = Mesh_Curve5.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve121 = Mesh_Curve6.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve122 = Mesh_Curve7.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve123 = Mesh_Curve8.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve124 = Mesh_Curve9.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve125 = Mesh_Curve10.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve126 = Mesh_Curve11.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve127 = Mesh_Curve12.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve128 = Mesh_Curve13.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve129 = Mesh_Curve14.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve130 = Mesh_Curve15.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve131 = Mesh_Curve16.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve132 = Mesh_Curve17.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve133 = Mesh_Curve18.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve134 = Mesh_Curve19.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve135 = Mesh_Curve20.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve136 = Mesh_Curve21.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve137 = Mesh_Curve22.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve138 = Mesh_Curve23.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve139 = Mesh_Curve24.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve140 = Mesh_Curve25.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve141 = Mesh_Curve26.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve142 = Mesh_Curve27.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve143 = Mesh_Curve28.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve144 = Mesh_Curve29.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve145 = Mesh_Curve30.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve146 = Mesh_Curve49.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve147 = Mesh_Curve44.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve148 = Mesh_Curve43.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve149 = Mesh_Curve34.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve150 = Mesh_Curve48.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve151 = Mesh_Curve36.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve152 = Mesh_Curve47.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve153 = Mesh_Curve38.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve154 = Mesh_Curve39.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve155 = Mesh_Curve40.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve156 = Mesh_Curve41.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve157 = Mesh_Curve42.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve158 = Mesh_Curve45.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve159 = Mesh_Curve46.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve160 = Mesh_Curve50.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve161 = Mesh_Curve51.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve162 = Mesh_Curve52.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve163 = Mesh_Curve53.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve164 = Mesh_Curve54.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve165 = Mesh_Curve55.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve166 = Mesh_Curve56.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve167 = Mesh_Curve57.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve168 = Mesh_Curve58.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve169 = Mesh_Curve59.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve170 = Mesh_Curve60.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve171 = Mesh_Curve61.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve172 = Mesh_Curve62.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve173 = Mesh_Curve63.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve174 = Mesh_Curve64.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve175 = Mesh_Curve65.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve176 = Mesh_Curve66.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve177 = Mesh_Curve67.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve178 = Mesh_Curve68.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve179 = Curve16.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve180 = Curve21.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve181 = Curve45.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve182 = Curve23.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve183 = Curve24.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve184 = Curve25.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point135 = Mesh_Point7.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point136 = Mesh_Point8.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point137 = Mesh_Point9.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point138 = Point44.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point139 = Point45.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point140 = Point46.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point141 = Point47.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point142 = Point48.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point143 = Point49.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point144 = P_A1_r0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point145 = P_A1_l1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point146 = P_A1_l0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point147 = P_A1_c1_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point148 = P_A1_c0_top.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Point149 = Mesh_Point10.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve185 = Curve26.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve186 = Curve27.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve187 = Curve28.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve188 = Curve29.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve189 = Curve30.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve190 = Curve31.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve191 = Curve32.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve192 = Curve33.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve193 = Curve34.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve194 = Curve35.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve195 = Curve36.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve196 = Curve37.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve197 = Curve38.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve198 = Curve39.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve199 = Curve40.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve200 = Curve41.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve201 = Curve42.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve202 = Curve43.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve203 = Curve44.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve204 = L_1_top_c.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve205 = L_1_colA_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve206 = L_1_bot_c.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve207 = L_1_Y_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve208 = L_2_ColA_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve209 = L_2_top_c.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve210 = L_2_ColB_v.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve211 = L_2_bot_c.copyMirror(P_C1_c0,Vector3d(0, 1, 0));
Mesh_Curve212 = Curve15.copyMirror(P_C1_c0,Vector3d(0, 1, 0));

DefaultName(typePoint,"Point",20,"");
DefaultName(typeGuideCurve,"Curve",11,"");
DefaultName(typeFeatureEdge,"FEdge",11,"");

wetted_surface.add(Pl100);
wetted_surface.add(Pl81);
wetted_surface.add(Pl83);
wetted_surface.add(Pl85);
wetted_surface.add(Pl87);
wetted_surface.add(Pl90);
wetted_surface.add(Pl91);
wetted_surface.add(Pl93);
wetted_surface.add(Pl94);
wetted_surface.add(Pl95);
wetted_surface.add(Pl96);
wetted_surface.add(Pl97);
wetted_surface.add(Pl98);
wetted_surface.add(Pl99);
wetted_surface.add(Pl84);
wetted_surface.add(Pl89);

wetted_surface.front.wetSurface = outer_wet;

LC1.generateAppliedLoads();
GenieRules.Meshing.faceMeshStrategy = PatchSurfQuadMesher;

Md_def.setDefault();


LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);
GenieRules.Meshing.superElementType = 4;
Meshing_panels.execute();
ExportMeshFem().DoExport("T4.FEM");

/*
Export mesh of only full wetted surface
 */

Meshing_panels.step(1).subset = wetted_surface;
SimplifyTopology();
SplitPeriodicGeometry();
GenieRules.Meshing.superElementType = 3;
Meshing_panels.execute();
ExportMeshFem().DoExport("T3.FEM");


/*
Generate mesh of free surface for 2nd order analysis
 */

FS_density = MeshDensity(FS_mesh_length);
FS_density.growthRate = 1.04;

DefaultName(typeFlatPlate,"Pl_fs",1,"");
DefaultName(typePoint,"Point_fs",1,"");
DefaultName(typeGuideCurve,"Curve_fs",4,"");
DefaultName(typeFeatureEdge,"FEdge_fs",1,"");

Point_fs1 = Point(Mesh_radius, 0, 0);
Point_fs2 = Point(-Mesh_radius, 0, 0);
Point_fs3 = Point(0, Mesh_radius, 0);
Curve_fs4 = CreateCircularArcFromThreePoints(Point_fs2, Point_fs3, Point_fs1);
Curve_fs5 = CreateLineTwoPoints(Point_fs2, Curve_fs4.end2);

Pl_fs1 = CoverCurves(Curve_fs4,Curve_fs5);

Validate(Pl_fs1.primitivePartCount == 15);
Pl_fs1.explode(IndexedNameMask(2));
Validate(Pl_fs2);
Validate(Pl_fs3);
Validate(Pl_fs4);
Validate(Pl_fs5);
Validate(Pl_fs6);
Validate(Pl_fs7);
Validate(Pl_fs8);
Validate(Pl_fs9);
Validate(Pl_fs10);
Validate(Pl_fs11);
Validate(Pl_fs12);
Validate(Pl_fs13);
Validate(Pl_fs14);
Validate(Pl_fs15);
Validate(Pl_fs16);
Delete(Pl_fs6);
Delete(Pl_fs14);
Delete(Pl_fs12);
Delete(Pl_fs13);
Delete(Pl_fs4);
Delete(Pl_fs9);
Delete(Pl_fs10);
Delete(Pl_fs8);
Delete(Pl_fs7);
Delete(Pl_fs15);
Delete(Pl_fs2);
Delete(Pl_fs5);
Delete(Pl_fs3);
Delete(Pl_fs16);


Curve_fs6 = Curve_fs5.divideAt(0.5);
Curve_fs7 = CreateLineTwoPoints(Curve_fs5.end2, Point_fs3);
Curve_fs8 = Curve_fs7.divideAt(0.5);
Curve_fs9 = Curve_fs5.divideAt(0.5);
Curve_fs10 = Curve_fs6.divideAt(0.5);
Curve_fs11 = CreateCircularArcFromThreePoints(Curve_fs5.end2, Curve_fs7.end2, Curve_fs10.end1);
FEdge_fs1 = FeatureEdge(Curve_fs7);
FEdge_fs2 = FeatureEdge(Curve_fs8);
FEdge_fs3 = FeatureEdge(Curve_fs11);

Curve_fs12 = Curve_fs5.divideAt(0.5);
Curve_fs13 = Curve_fs8.divideAt(0.5);
Curve_fs14 = Curve_fs10.divideAt(0.5);
Curve_fs15 = Curve_fs7.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),45);
Curve_fs16 = Curve_fs8.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),45);
Curve_fs17 = Curve_fs13.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),45);
Curve_fs18 = Curve_fs7.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),-45);
Curve_fs19 = Curve_fs8.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),-45);
Curve_fs20 = Curve_fs13.copyRotate(Curve_fs6.curvePoint(1),Vector3d(0, 0, 1),-45);
Curve_fs21 = CreateCircularArcFromThreePoints(Curve_fs12.end1, Curve_fs13.end1, Curve_fs10.end2);
FEdge_fs4 = FeatureEdge(Curve_fs15);
FEdge_fs5 = FeatureEdge(Curve_fs16);
FEdge_fs6 = FeatureEdge(Curve_fs17);
FEdge_fs7 = FeatureEdge(Curve_fs18);
FEdge_fs8 = FeatureEdge(Curve_fs19);
FEdge_fs9 = FeatureEdge(Curve_fs20);
FEdge_fs10 = FeatureEdge(Curve_fs21);

Pl_fs11.front.wetSurface = outer_wet;
Free_surface = Set();
Free_surface.add(Pl_fs11);

LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

Meshing_panels.step(1).subset = Free_surface;

Pl_fs11.meshDensity = FS_density;
GenieRules.Meshing.superElementType = 73;
Meshing_panels.execute();
ExportMeshFem().DoExport("T73.FEM");

/*
Generate mesh of internal lids for removing irregular frequencies
 */

DefaultName(typeFlatPlate,"Pl_int_lid",1,"");
DefaultName(typePoint,"Point_int_lid",1,"");
DefaultName(typeGuideCurve,"Curve_int_lid",1,"");

Curve_int_lid1 = Curve39.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));
Curve_int_lid2 = Curve40.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));
Curve_int_lid3 = Curve41.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));
Curve_int_lid4 = Curve42.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));
Curve_int_lid5 = Curve43.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));
Curve_int_lid6 = Curve44.copyTranslate(Vector3d(0 m,0 m,Draft-Col_height));

Pl_int_lid1 = CoverCurves(Curve_int_lid2,Curve_int_lid5,Curve_int_lid6);
Pl_int_lid2 = CoverCurves(Curve_int_lid3,Curve_int_lid4,Curve_int_lid1);
Pl_int_lid1.flipNormal();
Internal_lid = Set();
Internal_lid.add(Pl_int_lid1);
Internal_lid.add(Pl_int_lid2);
Internal_lid.meshDensity = Lid_density;
Internal_lid.front.wetSurface = outer_wet;

Meshing_panels.step(1).subset = Internal_lid;
GenieRules.Meshing.superElementType = 23;

Meshing_panels.execute();
ExportMeshFem().DoExport("T23.FEM");

/*
Generate compartment model
 */

Surf_lid = Permeable(true);
Pl_int_lid2.permeability = Surf_lid;
Pl_int_lid1.permeability = Surf_lid;

cm = CompartmentManager();

LC2 = LoadCase();
LC2.compartment(cm.compartment(Point(-3.924695148 m,28.00925713 m,13.25 m))).globalIntensity = DummyHydroPressure();
LC3 = LoadCase();
LC3.compartment(cm.compartment(Point(-4.724306289 m,31.67520684 m,-18.13333333 m))).globalIntensity = DummyHydroPressure();
LC4 = LoadCase();
LC4.compartment(cm.compartment(Point(-6.012776705 m,26.01051759 m,-18.13333333 m))).globalIntensity = DummyHydroPressure();
LC5 = LoadCase();
LC5.compartment(cm.compartment(Point(-12.18511773 m,-34.11349708 m,-18.13333333 m))).globalIntensity = DummyHydroPressure();
LC6 = LoadCase();
LC6.compartment(cm.compartment(Point(-16.05152461 m,-29.17853138 m,-18.13333333 m))).globalIntensity = DummyHydroPressure();
LC7 = LoadCase();
LC7.compartment(cm.compartment(Point(-16.10746603 m,-33.18964185 m,13.25 m))).globalIntensity = DummyHydroPressure();
LC8 = LoadCase();
LC8.compartment(cm.compartment(Point(-16.48989889 m,39.41649063 m,3 m))).globalIntensity = DummyHydroPressure();
LC9 = LoadCase();
LC9.compartment(cm.compartment(Point(-16.60742947 m,-37.24371515 m,3 m))).globalIntensity = DummyHydroPressure();
LC10 = LoadCase();
LC10.compartment(cm.compartment(Point(-17.4590317 m,1.702100158 m,-18.88 m))).globalIntensity = DummyHydroPressure();
LC11 = LoadCase();
LC11.compartment(cm.compartment(Point(-17.17645912 m,-40.59475789 m,3 m))).globalIntensity = DummyHydroPressure();
LC12 = LoadCase();
LC12.compartment(cm.compartment(Point(-17.17645912 m,36.63524211 m,3 m))).globalIntensity = DummyHydroPressure();
LC13 = LoadCase();
LC13.compartment(cm.compartment(Point(-18.04786303 m,-28.90408675 m,-18.88 m))).globalIntensity = DummyHydroPressure();
LC14 = LoadCase();
LC14.compartment(cm.compartment(Point(-20.689412 m,-31.19432378 m,13.25 m))).globalIntensity = DummyHydroPressure();
LC15 = LoadCase();
LC15.compartment(cm.compartment(Point(-20.60608412 m,44.22607218 m,3 m))).globalIntensity = DummyHydroPressure();
LC16 = LoadCase();
LC16.compartment(cm.compartment(Point(-20.83936815 m,-32.95159337 m,3 m))).globalIntensity = DummyHydroPressure();
LC17 = LoadCase();
LC17.compartment(cm.compartment(Point(-21.15791312 m,37.22727717 m,8.25 m))).globalIntensity = DummyHydroPressure();
LC18 = LoadCase();
LC18.compartment(cm.compartment(Point(-21.79036661 m,-39.15898115 m,8.25 m))).globalIntensity = DummyHydroPressure();
LC19 = LoadCase();
LC19.compartment(cm.compartment(Point(-21.79265764 m,-39.16485352 m,1 m))).globalIntensity = DummyHydroPressure();
LC20 = LoadCase();
LC20.compartment(cm.compartment(Point(-21.79265764 m,38.06514648 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC21 = LoadCase();
LC21.compartment(cm.compartment(Point(-21.86900672 m,-39.06786844 m,-10.38 m))).globalIntensity = DummyHydroPressure();
LC22 = LoadCase();
LC22.compartment(cm.compartment(Point(-21.94276225 m,38.23808953 m,-10.38 m))).globalIntensity = DummyHydroPressure();
LC23 = LoadCase();
LC23.compartment(cm.compartment(Point(-22.08650764 m,-43.94989766 m,3 m))).globalIntensity = DummyHydroPressure();
LC24 = LoadCase();
LC24.compartment(cm.compartment(Point(-22.10895116 m,33.29601412 m,3 m))).globalIntensity = DummyHydroPressure();
LC25 = LoadCase();
LC25.compartment(cm.compartment(Point(-23.00694789 m,1.803839581 m,-18.88 m))).globalIntensity = DummyHydroPressure();
LC26 = LoadCase();
LC26.compartment(cm.compartment(Point(-23.32454036 m,32.70676362 m,-2 m))).globalIntensity = DummyHydroPressure();
LC27 = LoadCase();
LC27.compartment(cm.compartment(Point(-23.72427661 m,-27.47563638 m,-18.88 m))).globalIntensity = DummyHydroPressure();
LC28 = LoadCase();
LC28.compartment(cm.compartment(Point(-24.27413853 m,-33.49707848 m,3 m))).globalIntensity = DummyHydroPressure();
LC29 = LoadCase();
LC29.compartment(cm.compartment(Point(-24.27413853 m,43.73292152 m,3 m))).globalIntensity = DummyHydroPressure();
LC30 = LoadCase();
LC30.compartment(cm.compartment(Point(-26.08347489 m,-43.09318711 m,3 m))).globalIntensity = DummyHydroPressure();
LC31 = LoadCase();
LC31.compartment(cm.compartment(Point(-26.77229381 m,34.82554875 m,-2 m))).globalIntensity = DummyHydroPressure();
LC32 = LoadCase();
LC32.compartment(cm.compartment(Point(-26.86221545 m,-42.28726309 m,3 m))).globalIntensity = DummyHydroPressure();
LC33 = LoadCase();
LC33.compartment(cm.compartment(Point(-27.62675202 m,38.81925117 m,-2 m))).globalIntensity = DummyHydroPressure();
LC34 = LoadCase();
LC34.compartment(cm.compartment(Point(-27.69026694 m,-38.31217305 m,3 m))).globalIntensity = DummyHydroPressure();
LC35 = LoadCase();
LC35.compartment(cm.compartment(Point(14.51558152 m,-12.06135168 m,-19.62666667 m))).globalIntensity = DummyHydroPressure();
LC36 = LoadCase();
LC36.compartment(cm.compartment(Point(17.47848924 m,18.85641686 m,-18.13333333 m))).globalIntensity = DummyHydroPressure();
LC37 = LoadCase();
LC37.compartment(cm.compartment(Point(19.55601888 m,-15.73551455 m,-19.62666667 m))).globalIntensity = DummyHydroPressure();
LC38 = LoadCase();
LC38.compartment(cm.compartment(Point(24.63528531 m,8.315850735 m,-19.62666667 m))).globalIntensity = DummyHydroPressure();
LC39 = LoadCase();
LC39.compartment(cm.compartment(Point(39.31434448 m,0.1268765731 m,3 m))).globalIntensity = DummyHydroPressure();
LC40 = LoadCase();
LC40.compartment(cm.compartment(Point(40.02092648 m,-3.672263095 m,3 m))).globalIntensity = DummyHydroPressure();
LC41 = LoadCase();
LC41.compartment(cm.compartment(Point(40.79966704 m,-4.478187108 m,3 m))).globalIntensity = DummyHydroPressure();
LC42 = LoadCase();
LC42.compartment(cm.compartment(Point(42.6090034 m,5.117921525 m,3 m))).globalIntensity = DummyHydroPressure();
LC43 = LoadCase();
LC43.compartment(cm.compartment(Point(44.0918491 m,0.4718293734 m,-10.38 m))).globalIntensity = DummyHydroPressure();
LC44 = LoadCase();
LC44.compartment(cm.compartment(Point(44.79663429 m,-5.334897664 m,3 m))).globalIntensity = DummyHydroPressure();
LC45 = LoadCase();
LC45.compartment(cm.compartment(Point(45.0704127 m,-0.5175745489 m,8.25 m))).globalIntensity = DummyHydroPressure();
LC46 = LoadCase();
LC46.compartment(cm.compartment(Point(45.0904843 m,-0.5498535161 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC47 = LoadCase();
LC47.compartment(cm.compartment(Point(46.08908243 m,5.650691695 m,3 m))).globalIntensity = DummyHydroPressure();
LC48 = LoadCase();
LC48.compartment(cm.compartment(Point(49.70668281 m,-1.979757888 m,3 m))).globalIntensity = DummyHydroPressure();
LC49 = LoadCase();
LC49.compartment(cm.compartment(Point(50.25216792 m,1.455012495 m,-2 m))).globalIntensity = DummyHydroPressure();


Meshing_panels.step(1).subset = Free_surface;
Meshing_panels.step(1).excludeIncludeMeshSubsetOption = anExcludeMeshSubset;
SimplifyTopology();
SplitPeriodicGeometry();

GenieRules.Meshing.superElementType = 5;

Meshing_panels.execute();
ExportMeshFem().DoExport("T5.FEM");

/*
Create Morison model
 */

Col = PipeSection(Col_diameter, t1);
Pon = BoxSection(Pon_height, Pon_width, t1, t1, t1);
Col.setDefault();

Bm1 = StraightBeam(Mesh_Curve57.end2, P_A_c2);
Bm2 = StraightBeam(Mesh_Curve61.end2, P_B_c2);
Bm3 = StraightBeam(Mesh_Curve167.end2, Mesh_Point55);
Point_int_lid1 = Mesh_Point87.copyTranslate(Vector3d(0 m,0 m,Pon_height/2));
Point_int_lid2 = Mesh_Point81.copyTranslate(Vector3d(0 m,0 m,Pon_height/2));
Pon.setDefault();
Bm4 = StraightBeam(Point_int_lid2, Point_int_lid1);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Bm4, "Bm5");
ModelTransformer(MyModelTransformerMap).copyRotate(Mesh_Point53, Vector3d(0, 0, 1), 60);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Bm5, "Bm6");
ModelTransformer(MyModelTransformerMap).copyRotate(P_A_c0, Vector3d(0, 0, 1), 60);
Morison = Set();
Morison.add(Bm1);
Morison.add(Bm2);
Morison.add(Bm3);
Morison.add(Bm4);
Morison.add(Bm5);
Morison.add(Bm6);


Delete(el_elliptic);
Delete(el_inter);
Delete(Half_pon_height);
Delete(Md_def);
Delete(Lid_density);
Delete(FS_density);

GenieRules.Meshing.autoSimplifyTopology = false;
GenieRules.Meshing.autoSplitPeriodicGeometry = false;
GenieRules.Meshing.superElementType = 6;

Meshing_panels.step(1).subset = Morison;
Meshing_panels.step(1).excludeIncludeMeshSubsetOption = anIncludeMeshSubset;
Meshing_panels.execute();
ExportMeshFem().DoExport("T6.FEM");