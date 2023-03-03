// Parametric modelling of the UMaine VolturnUS-S reference platform
// author: Stian Brurås
GenieRules.Units.setInputUnit(Angle, "deg");
GenieRules.Tolerances.useTolerantModelling = true;
GenieRules.Compatibility.version = "V8.3-31";
GenieRules.Tolerances.useTolerantModelling = true;
GenieRules.Tolerances.angleTolerance = 2 deg;
GenieRules.Meshing.autoSimplifyTopology = true;
GenieRules.Meshing.eliminateInternalEdges = true;
GenieRules.BeamCreation.DefaultCurveOffset = ReparameterizedBeamCurveOffset();
GenieRules.Transformation.DefaultConnectedCopy = false;


/*
Generate and mesh mean wetted surface using XZ as symmetry plane
 */



mesh_size = 0.5 m;  // Mesh size
Mesh_radius = 150 m;  // Radius for free surface mesh
FS_mesh_length = 5 m;  // Mesh size for free surface mesh
mesh_lid = 0.5 m;  // Mesh size for internal lids

d1 = 10.0 m; // Diameter of center column
d2 = 12.5 m; // Diameter of radial columns
d_brace = 0.91 m; // Diameter of braces

w_pontoon = 12.5 m; // Horizontal width of pontoon
h_pontoon = 7.00 m;  // Vertical height of pontoon

draft = 20.00 m; // Draft of floater in operating condition
freeboard = 15.00 m; // Freeboard of floater in operating condition
fairlead_depth = 14.00 m; // Fairlead_depth

cc_distance = 51.75 m;  // Center-center distance radial to center column

deck_collition_low = 3 m;
deck_collition_high = 5 m;
damage_extent = 1.5 m;

t1 = 40.0 mm; // Steel thickness, external plates
t2 = 25.0 mm; // Steel thickness, internal plates


Steel = MaterialLinear(200e6, 7850 kg/m^3, 2.1e+11 Pa, 0.3, 1.2e-05 delC^-1, 0.03 N*s/m);
Steel.setDefault();

Tck1 = Thickness(t1);
Tck2 = Thickness(t2);
Tck1.setDefault();

circumference2 = Math.Pi*d2;
r2 = d2/2;

Md_def = MeshDensity(mesh_size);

// Wet Surfaces
outer_wet = WetSurface();
wetted_surface = Set();

//Meshing rules
GenieRules.Meshing.elementType = mp1stOrder;
GenieRules.Meshing.superElementType = 1;
GenieRules.Meshing.autoSimplifyTopology = true;
GenieRules.Meshing.autoSplitPeriodicGeometry = true;
GenieRules.Meshing.preference(mpPreferRectangularMesh, true);
GenieRules.Meshing.preference(mpAllowTriangularElements, true);
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
GenieRules.Meshing.edgeMeshStrategy = LinearDistributionEdge;
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

//Tolerances Rules
GenieRules.Tolerances.angleTolerance = 2 deg;
GenieRules.Tolerances.pointTolerance = 0.01 m;
GenieRules.Tolerances.useTolerantModelling = true;

// Generate dummy hydro load case for HydroD
LC1 = DummyHydroLoadCase();
LC1.setFemLoadcase(1);
LC1.designCondition(lcOperating);
LC1.wetSurface = outer_wet;

Meshing_panels = Analysis(true);
Meshing_panels.add(MeshActivity());

Point0 = Point(0m, 0m, -draft);
Point1 = Point(cc_distance, 0m, -draft);
Point2 = Point(cc_distance, d2/2, -draft);
Point3 = Point(cc_distance+d2/2, 0m, -draft);
Point4 = Point(cc_distance-d2/2, 0m, -draft);
Curve1 = CreateCircularArcFromThreePoints(Point3, Point2, Point4);
Curve2 = CreateLineTwoPoints(Point4, Point1);
Curve3 = CreateLineTwoPoints(Point1, Point3);
x_trans = cc_distance*(1+Math.cos(60*Math.PI/180));
y_trans = cc_distance*Math.sin(60*Math.PI/180);
Point5 = Point1.copyTranslate(Vector3d(-x_trans, y_trans, 0));
Point6 = Point2.copyTranslate(Vector3d(-x_trans, y_trans, 0));
Point7 = Point3.copyTranslate(Vector3d(-x_trans, y_trans, 0));
Point8 = Point4.copyTranslate(Vector3d(-x_trans, y_trans, 0));
Curve4 = CreateCircleFromThreePoints(Point7, Point6, Point8);
temp = w_pontoon/2*Math.tan(120*Math.PI/180);
Point9 = Point(-temp, w_pontoon/2, -draft);

Curve5 = CreateLineTwoPoints(Point2, Point9);
Curve6 = CreateLineTwoPoints(Point0, Point5);
Curve7 = CreateLineTwoPoints(Curve1.end1, Curve6.end1);
temp = d2/2*Math.cos(Math.asin(w_pontoon/d2));

r = w_pontoon/(2*Math.sin(60*Math.PI/180));
Point9 = Point(r*Math.cos(60*Math.PI/180), w_pontoon/2, -draft);
Delete(Curve5);
Curve8 = CreateLineTwoPoints(Point2, Point9);
x = cc_distance*Math.cos(60*Math.PI/180)+d2/2*Math.cos(150*Math.PI/180);
y= cc_distance*Math.sin(60*Math.PI/180)+d2/2*Math.sin(150*Math.PI/180);

Point10 = Point(-x,y,-draft);
Curve9 = CreateLineTwoPoints(Curve8.end2, Point10);
x = cc_distance*Math.cos(60*Math.PI/180)+d2/2*Math.cos(-30*Math.PI/180);
y= cc_distance*Math.sin(60*Math.PI/180)+d2/2*Math.sin(-30*Math.PI/180);

Point11 = Point(-x,y,-draft);
x = w_pontoon/(2*Math.tan(60*Math.PI/180))+Math.sqrt(r*r - (w_pontoon/2)*(w_pontoon/2));

Point12 = Point(-x, 0, -draft);
Curve10 = CreateLineTwoPoints(Point11, Point12);
Curve5 = CreateLineTwoPoints(Point0, Point12);
Point13 = Point0.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point14 = Point1.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point15 = Point2.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point16 = Point3.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point17 = Point4.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point18 = Point5.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point19 = Point6.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point20 = Point7.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point21 = Point8.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point22 = Point9.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point23 = Point10.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point24 = Point12.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Point25 = Point11.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve11 = Curve1.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve12 = Curve2.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve13 = Curve3.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve14 = Curve4.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve15 = Curve8.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve16 = Curve6.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve17 = Curve7.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve18 = Curve9.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve19 = Curve10.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve20 = Curve5.copyTranslate(Vector3d(0m, 0m, h_pontoon));
Curve21 = CreateLineTwoPoints(Point4, Point17);
Curve22 = CreateLineTwoPoints(Point2, Point15);
Curve23 = CreateLineTwoPoints(Point3, Point16);
Curve24 = CreateLineTwoPoints(Point9, Point22);
Curve25 = CreateLineTwoPoints(Point0, Point13);
Curve26 = CreateLineTwoPoints(Point12, Point24);
Curve27 = CreateLineTwoPoints(Point10, Point23);
Curve28 = CreateLineTwoPoints(Point11, Point25);
Curve29 = CreateLineTwoPoints(Point5, Point18);
Point26 = Point(-d1/2, 0m, -draft+h_pontoon);
Point27 = Point(d1/2, 0m, -draft+h_pontoon);
Point28 = Point(0m, d1/2, -draft+h_pontoon);
Curve30 = CreateCircularArcFromThreePoints(Point26, Point28, Point27);
Curve31 = Curve1.divideAt(0.5);
Curve32 = Curve11.divideAt(0.5);

Pl1 = CoverCurves(Curve1,Curve11,Curve22,Curve23);
Pl1.flipNormal();
Pl2 = CoverCurves(Curve1,Curve2,Curve3,Curve31);
Pl2.flipNormal();
Curve33 = CreateLineTwoPoints(Point0, Curve9.end1);
Curve34 = CreateLineTwoPoints(Point13, Curve15.end2);

Pl3 = CoverCurves(Curve1,Curve8,Curve7,Curve33);
Pl4 = CoverCurves(Curve8,Curve15,Curve22,Curve24);
Pl5 = CoverCurves(Curve9,Curve18,Curve24,Curve27);
Pl4.flipNormal();
Pl5.flipNormal();
Curve35 = Curve4.divideAt(0.08333333381);
Curve36 = Curve35.divideAt(0.5454545465);
tmpArrayOfCurves = JoinMultipleCurves(Array(Curve4, Curve36), true);
Rename(tmpArrayOfCurves[0], "Curve36");
Pl6 = SweepCurve(Curve35, Curve28);
Curve37 = Curve36.divideAt(0.5000000014);
Curve38 = Curve6.divideAt(0.8792270531);
Pl7 = CoverCurves(Curve6,Curve9,Curve33,Curve37);
Pl7.flipNormal();
Pl8 = CoverCurves(Curve6,Curve10,Curve5,Curve36);
Pl8.flipNormal();
Pl9 = CoverCurves(Curve35,Curve37,Curve36);
Pl9.flipNormal();
Pl10 = CoverCurves(Curve10,Curve19,Curve26,Curve28);
Pl10.flipNormal();
Curve39 = Curve17.divideAt(0.2155172414);
Delete(Curve17);
Pl11 = CoverCurves(Curve15,Curve32,Curve34,Curve39);
Pl11.flipNormal();
Curve40 = Curve14.divideAt(0.08333333381);
Curve41 = Curve40.divideAt(0.5454545465);
tmpArrayOfCurves = JoinMultipleCurves(Array(Curve14, Curve41), true);
Rename(tmpArrayOfCurves[0], "Curve41");
Curve42 = Curve41.divideAt(0.5000000014);
Curve43 = Curve16.divideAt(0.8792270531);
Pl12 = CoverCurves(Curve16,Curve18,Curve34,Curve42);
Pl13 = CoverCurves(Curve16,Curve19,Curve20,Curve41);

Curve44 = Curve39.divideAt(0.8901098901);
Curve45 = Curve20.divideAt(0.692820323);
Point29 = Point21.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point30 = Point23.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point31 = Point14.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point32 = Point13.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point33 = Point20.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point34 = Point19.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point35 = Point17.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point36 = Point16.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point37 = Point15.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point38 = Point25.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point39 = Point26.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point40 = Point27.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point41 = Point28.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve46 = Curve11.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve47 = Curve12.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve48 = Curve13.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve49 = Curve40.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve50 = Curve20.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve51 = Curve30.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve52 = Curve32.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve53 = Curve42.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve54 = Curve41.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve55 = Curve44.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));

Curve56 = CreateLineTwoPoints(Curve13.end2, Curve48.end2);
Curve57 = CreateLineTwoPoints(Curve30.end2, Curve51.end2);
Curve58 = CreateLineTwoPoints(Curve18.end2, Curve49.end1);

Delete(Pl11);
Delete(Pl12);
Delete(Pl13);

Curve59 = Curve34.divideAt(0.692820323);
Curve60 = Curve30.divideAt(0.6666666676);
Curve61 = Curve16.divideAt(0.1098901099);
Curve62 = Curve30.divideAt(0.5000000111);
Pl14 = CoverCurves(Curve15,Curve32,Curve39,Curve59,Curve60);
Pl14.flipNormal();
Pl15 = CoverCurves(Curve18,Curve42,Curve59,Curve61,Curve62);
Pl16 = CoverCurves(Curve19,Curve30,Curve41,Curve45,Curve61);
Pl17 = SweepCurve(Curve11, Curve56);
Pl18 = SweepCurve(Curve32, Curve56);
Pl19 = SweepCurve(Curve30, Curve57);
Pl20 = SweepCurve(Curve60, Curve57);
Pl21 = SweepCurve(Curve62, Curve57);
Pl19.flipNormal();
Pl20.flipNormal();
Pl21.flipNormal();
Pl22 = SweepCurve(Curve40, Curve58);
Pl23 = SweepCurve(Curve42, Curve58);
Pl24 = SweepCurve(Curve41, Curve58);
Pl22.flipNormal();


wetted_surface.add(Pl1);
wetted_surface.add(Pl10);
wetted_surface.add(Pl14);
wetted_surface.add(Pl15);
wetted_surface.add(Pl16);
wetted_surface.add(Pl17);
wetted_surface.add(Pl18);
wetted_surface.add(Pl19);
wetted_surface.add(Pl2);
wetted_surface.add(Pl20);
wetted_surface.add(Pl21);
wetted_surface.add(Pl22);
wetted_surface.add(Pl23);
wetted_surface.add(Pl24);
wetted_surface.add(Pl3);
wetted_surface.add(Pl4);
wetted_surface.add(Pl5);
wetted_surface.add(Pl6);
wetted_surface.add(Pl7);
wetted_surface.add(Pl8);
wetted_surface.add(Pl9);
wetted_surface.setActive();
wetted_surface.front.wetSurface = outer_wet;


FEdge1 = FeatureEdge(Curve30.end2, Curve19.project(Curve30.end2));
FEdge2 = FeatureEdge(Curve19.project(Curve30.end2), Curve10.project(Curve19.project(Curve30.end2)));
FEdge3 = FeatureEdge(Curve10.project(Curve19.project(Curve30.end2)), Curve6.project(Curve10.project(Curve19.project(Curve30.end2))));
FEdge5 = FeatureEdge(Point(-2.5 m,4.330127019 m,-20 m), Curve9.project(Point(-2.5 m,4.330127019 m,-20 m)));
FEdge6 = FeatureEdge(Curve9.project(Point(-2.5 m,4.330127019 m,-20 m)), Curve18.project(Curve9.project(Point(-2.5 m,4.330127019 m,-20 m))));
FEdge7 = FeatureEdge(Curve18.project(Curve9.project(Point(-2.5 m,4.330127019 m,-20 m))), Point(-2.49999988 m,4.330127088 m,-13 m));
FEdge8 = FeatureEdge(Curve60.end2, Curve15.project(Curve60.end2));
FEdge9 = FeatureEdge(Curve15.project(Curve60.end2), Curve8.project(Curve15.project(Curve60.end2)));
FEdge10 = FeatureEdge(Curve8.project(Curve15.project(Curve60.end2)), Curve7.project(Curve8.project(Curve15.project(Curve60.end2))));
FEdge11 = FeatureEdge(Point(45.5 m,0 m,-13 m), Curve15.project(Point(45.5 m,0 m,-13 m)));
FEdge12 = FeatureEdge(Curve15.project(Point(45.5 m,0 m,-13 m)), Curve8.project(Curve15.project(Point(45.5 m,0 m,-13 m))));
FEdge13 = FeatureEdge(Curve8.project(Curve15.project(Point(45.5 m,0 m,-13 m))), Curve7.project(Curve8.project(Curve15.project(Point(45.5 m,0 m,-13 m)))));
FEdge14 = FeatureEdge(Point(-22.75 m,39.40415587 m,-13 m), Curve18.project(Point(-22.75 m,39.40415587 m,-13 m)));
FEdge15 = FeatureEdge(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m)), Curve9.project(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m))));
FEdge16 = FeatureEdge(Curve9.project(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m))), Curve6.project(Curve9.project(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m)))));
FEdge17 = FeatureEdge(Curve6.project(Curve9.project(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m)))), Curve10.project(Curve6.project(Curve9.project(Curve18.project(Point(-22.75 m,39.40415587 m,-13 m))))));
FEdge19 = FeatureEdge(Point(-28.16265877 m,36.27915587 m,-20 m), Curve19.project(Point(-28.16265877 m,36.27915587 m,-20 m)));
FEdge20 = FeatureEdge(Curve19.project(Point(-28.16265877 m,36.27915587 m,-20 m)), Point(-22.75 m,39.40415587 m,-13 m));




GenieRules.Meshing.faceMeshStrategy = PatchSurfQuadMesher;
Md_def.setDefault();
LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

Meshing_panels.execute();
ExportMeshFem().DoExport("T1.FEM");


/*
Modeling of the structure above MWL, still applying symmetry
 */

Curve81 = Curve46.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve82 = Curve47.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve83 = Curve48.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve84 = Curve49.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve85 = Curve50.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve86 = Curve51.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve87 = Curve52.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve88 = Curve53.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve89 = Curve54.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve90 = Curve55.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve91 = CreateLineTwoPoints(Curve46.end1, Curve81.project(Curve46.end1));
Pl25 = SweepCurve(Curve46, Curve91);
Pl26 = SweepCurve(Curve49, Curve91);
Pl27 = SweepCurve(Curve51, Curve91);
Pl28 = SweepCurve(Curve52, Curve91);
Pl29 = SweepCurve(Curve53, Curve91);
Pl30 = SweepCurve(Curve54, Curve91);
Pl25.flipNormal();
Pl28.flipNormal();
Pl29.flipNormal();
Pl30.flipNormal();
wetted_surface.add(Pl25);
wetted_surface.add(Pl26);
wetted_surface.add(Pl27);
wetted_surface.add(Pl28);
wetted_surface.add(Pl29);
wetted_surface.add(Pl30);
wetted_surface.front.wetSurface = outer_wet;

LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);
GenieRules.Meshing.superElementType = 2;
Meshing_panels.execute();
ExportMeshFem().DoExport("T2.FEM");

/*
Generate internal plates
 */

Tck2.setDefault();

Pl31 = CoverCurves(Curve81,Curve82,Curve83,Curve87);
Pl32 = CoverCurves(Curve85,Curve86,Curve90);
Pl32.flipNormal();
Pl33 = CoverCurves(Curve84,Curve88,Curve89);
Pl34 = CoverCurves(Curve21,Curve22,Curve31,Curve32);
Pl34.flipNormal();
Pl35 = CoverCurves(Curve24,Curve25,Curve33,Curve34,Curve59);
Pl35.flipNormal();
Pl36 = CoverCurves(Curve20,Curve30,Curve44,Curve60,Curve62);
Pl36.flipNormal();
Pl37 = CoverCurves(Curve5,Curve20,Curve25,Curve26,Curve45);
Pl38 = SweepCurve(Curve37, Curve27);
Pl39 = SweepCurve(Curve36, Curve27);

Curve92 = Curve46.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve93 = Curve47.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve94 = Curve48.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve95 = Curve49.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve96 = Curve50.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve97 = Curve51.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve98 = Curve52.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve99 = Curve53.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve100 = Curve54.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve101 = Curve55.copyTranslate(Vector3d(0 m,0 m, deck_collition_high));
Curve102 = Curve46.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve103 = Curve47.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve104 = Curve48.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve105 = Curve49.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve106 = Curve50.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve107 = Curve51.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve108 = Curve52.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve109 = Curve53.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve110 = Curve54.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));
Curve111 = Curve55.copyTranslate(Vector3d(0 m,0 m, -deck_collition_low));

Pl40 = CoverCurves(Curve105,Curve109,Curve110);
Pl41 = CoverCurves(Curve106,Curve107,Curve111);
Pl42 = CoverCurves(Curve102,Curve103,Curve104,Curve108);
Pl43 = CoverCurves(Curve92,Curve93,Curve94,Curve98);
Pl44 = CoverCurves(Curve96,Curve97,Curve101);
Pl45 = CoverCurves(Curve95,Curve99,Curve100);


Curve112 = Curve104.divideAt(0.76);
Curve113 = Curve112.copyRotate(Curve103.curvePoint(2),Vector3d(0, 0, 1),90*1);
Curve114 = Curve112.copyRotate(Curve103.curvePoint(2),Vector3d(0, 0, 1),90*2);
Curve115 = CreateCircularArcFromThreePoints(Curve104.end2, Curve113.end1, Curve114.end1);
Curve116 = CreateLineTwoPoints(Curve112.end1, Curve94.project(Curve112.end1));


Pl46 = SweepCurve(Curve112, Curve116);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl46, "Pl47");
ModelTransformer(MyModelTransformerMap).copyRotate(Curve104.curvePoint(1), Vector3d(0, 0, 1), 45, 4);
Pl51 = SweepCurve(Curve115, Curve116);

Curve117 = Curve106.divideAt(0.7);
Curve118 = Curve117.copyRotate(Curve111.curvePoint(2),Vector3d(0, 0, 1),-90*1);
Curve119 = Curve117.copyRotate(Curve111.curvePoint(2),Vector3d(0, 0, 1),-90*2);
Curve120 = CreateCircularArcFromThreePoints(Curve117.end1, Curve118.project(Curve117.end1), Curve119.project(Curve118.project(Curve117.end1)));
Curve121 = CreateLineTwoPoints(Curve117.end1, Curve96.project(Curve117.end1));

Pl52 = SweepCurve(Curve117, Curve121);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl52, "Pl53");
ModelTransformer(MyModelTransformerMap).copyRotate(Curve111.curvePoint(2), Vector3d(0, 0, 1), -45, 4);

Curve122 = CreateLineTwoPoints(Point(-22.75 m,39.40415587 m,-3 m), Point(-29 m,50.22947342 m,-3 m));
Curve123 = Curve122.divideAt(0.12);
Curve124 = Curve122.copyRotate(Curve43.curvePoint(2),Vector3d(0, 0, 1),90*1);
Curve125 = Curve122.copyRotate(Curve43.curvePoint(2),Vector3d(0, 0, 1),90*2);
Curve126 = CreateCircleFromThreePoints(Curve122.end2, Curve124.project(Curve122.end2), Curve125.project(Curve124.project(Curve122.end2)));
Curve127 = CreateLineTwoPoints(Point(-22.75 m,39.40415587 m,-3 m), Curve99.project(Point(-22.75 m,39.40415587 m,-3 m)));
Pl58 = SweepCurve(Curve122, Curve127);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl58, "Pl59");
ModelTransformer(MyModelTransformerMap).copyRotate(Curve43.curvePoint(2), Vector3d(0, 0, 1), 45, 7);

Pl66 = SweepCurve(Curve126, Curve127);
Curve128 = Curve39.divideAt(0.5);
Curve129 = Curve61.divideAt(0.5);
Curve131 = CreateLineTwoPoints(Curve129.end1, Curve19.project(Curve129.end1));
Curve130 = CreateLineTwoPoints(Curve61.end2, Curve18.project(Curve61.end2));
Curve132 = CreateLineTwoPoints(Curve128.end1, Curve15.project(Curve128.end1));
Pl67 = SweepCurve(Curve131, Curve25);
Pl68 = SweepCurve(Curve130, Curve25);
Pl69 = SweepCurve(Curve132, Curve25);
Curve133 = Curve7.divideAt(0.2155172414);
Pl70 = SweepCurve(Curve133, Curve25);
Pl71 = SweepCurve(Curve6, Curve25);
Pl72 = SweepCurve(Curve120, Curve121);


Validate(Pl40.primitivePartCount == 9);
Pl40.explode(IndexedNameMask(73));
Validate(Pl73);
Validate(Pl74);
Validate(Pl75);
Validate(Pl76);
Validate(Pl77);
Validate(Pl78);
Validate(Pl79);
Validate(Pl80);
Validate(Pl81);
Validate(Pl41.primitivePartCount == 5);
Pl41.explode(IndexedNameMask(82));
Validate(Pl82);
Validate(Pl83);
Validate(Pl84);
Validate(Pl85);
Validate(Pl86);
Validate(Pl42.primitivePartCount == 5);
Pl42.explode(IndexedNameMask(87));
Validate(Pl87);
Validate(Pl88);
Validate(Pl89);
Validate(Pl90);
Validate(Pl91);
Validate(Pl43.primitivePartCount == 5);
Pl43.explode(IndexedNameMask(92));
Validate(Pl92);
Validate(Pl93);
Validate(Pl94);
Validate(Pl95);
Validate(Pl96);
Validate(Pl44.primitivePartCount == 5);
Pl44.explode(IndexedNameMask(97));
Validate(Pl97);
Validate(Pl98);
Validate(Pl99);
Validate(Pl100);
Validate(Pl101);
Validate(Pl45.primitivePartCount == 9);
Pl45.explode(IndexedNameMask(102));
Validate(Pl102);
Validate(Pl103);
Validate(Pl104);
Validate(Pl105);
Validate(Pl106);
Validate(Pl107);
Validate(Pl108);
Validate(Pl109);
Validate(Pl110);
/*
Mirror the structure and export full mesh
 */
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(FEdge1, "FEdge21");
MyModelTransformerMap.Add(FEdge2, "FEdge22");
MyModelTransformerMap.Add(FEdge3, "FEdge23");
MyModelTransformerMap.Add(FEdge5, "FEdge24");
MyModelTransformerMap.Add(FEdge6, "FEdge25");
MyModelTransformerMap.Add(FEdge7, "FEdge26");
MyModelTransformerMap.Add(FEdge8, "FEdge27");
MyModelTransformerMap.Add(FEdge9, "FEdge28");
MyModelTransformerMap.Add(FEdge10, "FEdge29");
MyModelTransformerMap.Add(FEdge11, "FEdge30");
MyModelTransformerMap.Add(FEdge12, "FEdge31");
MyModelTransformerMap.Add(FEdge13, "FEdge32");
MyModelTransformerMap.Add(FEdge14, "FEdge33");
MyModelTransformerMap.Add(FEdge15, "FEdge34");
MyModelTransformerMap.Add(FEdge16, "FEdge35");
MyModelTransformerMap.Add(FEdge17, "FEdge36");
MyModelTransformerMap.Add(FEdge19, "FEdge37");
MyModelTransformerMap.Add(FEdge20, "FEdge38");
MyModelTransformerMap.Add(Pl1, "Pl111");
MyModelTransformerMap.Add(Pl2, "Pl112");
MyModelTransformerMap.Add(Pl3, "Pl113");
MyModelTransformerMap.Add(Pl4, "Pl114");
MyModelTransformerMap.Add(Pl5, "Pl115");
MyModelTransformerMap.Add(Pl6, "Pl116");
MyModelTransformerMap.Add(Pl7, "Pl117");
MyModelTransformerMap.Add(Pl8, "Pl118");
MyModelTransformerMap.Add(Pl9, "Pl119");
MyModelTransformerMap.Add(Pl10, "Pl120");
MyModelTransformerMap.Add(Pl14, "Pl121");
MyModelTransformerMap.Add(Pl15, "Pl122");
MyModelTransformerMap.Add(Pl16, "Pl123");
MyModelTransformerMap.Add(Pl17, "Pl124");
MyModelTransformerMap.Add(Pl18, "Pl125");
MyModelTransformerMap.Add(Pl19, "Pl126");
MyModelTransformerMap.Add(Pl20, "Pl127");
MyModelTransformerMap.Add(Pl21, "Pl128");
MyModelTransformerMap.Add(Pl22, "Pl129");
MyModelTransformerMap.Add(Pl23, "Pl130");
MyModelTransformerMap.Add(Pl24, "Pl131");
MyModelTransformerMap.Add(Pl25, "Pl132");
MyModelTransformerMap.Add(Pl26, "Pl133");
MyModelTransformerMap.Add(Pl27, "Pl134");
MyModelTransformerMap.Add(Pl28, "Pl135");
MyModelTransformerMap.Add(Pl29, "Pl136");
MyModelTransformerMap.Add(Pl30, "Pl137");
MyModelTransformerMap.Add(Pl31, "Pl138");
MyModelTransformerMap.Add(Pl32, "Pl139");
MyModelTransformerMap.Add(Pl33, "Pl140");
MyModelTransformerMap.Add(Pl34, "Pl141");
MyModelTransformerMap.Add(Pl35, "Pl142");
MyModelTransformerMap.Add(Pl36, "Pl143");
MyModelTransformerMap.Add(Pl37, "Pl144");
MyModelTransformerMap.Add(Pl38, "Pl145");
MyModelTransformerMap.Add(Pl39, "Pl146");
MyModelTransformerMap.Add(Pl46, "Pl147");
MyModelTransformerMap.Add(Pl47, "Pl148");
MyModelTransformerMap.Add(Pl48, "Pl149");
MyModelTransformerMap.Add(Pl49, "Pl150");
MyModelTransformerMap.Add(Pl50, "Pl151");
MyModelTransformerMap.Add(Pl51, "Pl152");
MyModelTransformerMap.Add(Pl52, "Pl153");
MyModelTransformerMap.Add(Pl53, "Pl154");
MyModelTransformerMap.Add(Pl54, "Pl155");
MyModelTransformerMap.Add(Pl55, "Pl156");
MyModelTransformerMap.Add(Pl56, "Pl157");
MyModelTransformerMap.Add(Pl58, "Pl158");
MyModelTransformerMap.Add(Pl59, "Pl159");
MyModelTransformerMap.Add(Pl60, "Pl160");
MyModelTransformerMap.Add(Pl61, "Pl161");
MyModelTransformerMap.Add(Pl62, "Pl162");
MyModelTransformerMap.Add(Pl63, "Pl163");
MyModelTransformerMap.Add(Pl64, "Pl164");
MyModelTransformerMap.Add(Pl65, "Pl165");
MyModelTransformerMap.Add(Pl66, "Pl166");
MyModelTransformerMap.Add(Pl67, "Pl167");
MyModelTransformerMap.Add(Pl68, "Pl168");
MyModelTransformerMap.Add(Pl69, "Pl169");
MyModelTransformerMap.Add(Pl70, "Pl170");
MyModelTransformerMap.Add(Pl71, "Pl171");
MyModelTransformerMap.Add(Pl72, "Pl172");
MyModelTransformerMap.Add(Pl73, "Pl173");
MyModelTransformerMap.Add(Pl74, "Pl174");
MyModelTransformerMap.Add(Pl75, "Pl175");
MyModelTransformerMap.Add(Pl76, "Pl176");
MyModelTransformerMap.Add(Pl77, "Pl177");
MyModelTransformerMap.Add(Pl78, "Pl178");
MyModelTransformerMap.Add(Pl79, "Pl179");
MyModelTransformerMap.Add(Pl80, "Pl180");
MyModelTransformerMap.Add(Pl81, "Pl181");
MyModelTransformerMap.Add(Pl82, "Pl182");
MyModelTransformerMap.Add(Pl83, "Pl183");
MyModelTransformerMap.Add(Pl84, "Pl184");
MyModelTransformerMap.Add(Pl85, "Pl185");
MyModelTransformerMap.Add(Pl86, "Pl186");
MyModelTransformerMap.Add(Pl87, "Pl187");
MyModelTransformerMap.Add(Pl88, "Pl188");
MyModelTransformerMap.Add(Pl89, "Pl189");
MyModelTransformerMap.Add(Pl90, "Pl190");
MyModelTransformerMap.Add(Pl91, "Pl191");
MyModelTransformerMap.Add(Pl92, "Pl192");
MyModelTransformerMap.Add(Pl93, "Pl193");
MyModelTransformerMap.Add(Pl94, "Pl194");
MyModelTransformerMap.Add(Pl95, "Pl195");
MyModelTransformerMap.Add(Pl96, "Pl196");
MyModelTransformerMap.Add(Pl97, "Pl197");
MyModelTransformerMap.Add(Pl98, "Pl198");
MyModelTransformerMap.Add(Pl99, "Pl199");
MyModelTransformerMap.Add(Pl100, "Pl200");
MyModelTransformerMap.Add(Pl101, "Pl201");
MyModelTransformerMap.Add(Pl102, "Pl202");
MyModelTransformerMap.Add(Pl103, "Pl203");
MyModelTransformerMap.Add(Pl104, "Pl204");
MyModelTransformerMap.Add(Pl105, "Pl205");
MyModelTransformerMap.Add(Pl106, "Pl206");
MyModelTransformerMap.Add(Pl107, "Pl207");
MyModelTransformerMap.Add(Pl108, "Pl208");
MyModelTransformerMap.Add(Pl109, "Pl209");
MyModelTransformerMap.Add(Pl110, "Pl210");
ModelTransformer(MyModelTransformerMap).copyMirror(Point0, Vector3d(0, 1, 0));
Point42 = Point0.copyMirror(Point0,Vector3d(0, 1, 0));
Point43 = Point1.copyMirror(Point0,Vector3d(0, 1, 0));
Point44 = Point2.copyMirror(Point0,Vector3d(0, 1, 0));
Point45 = Point3.copyMirror(Point0,Vector3d(0, 1, 0));
Point46 = Point4.copyMirror(Point0,Vector3d(0, 1, 0));
Point47 = Point5.copyMirror(Point0,Vector3d(0, 1, 0));
Point48 = Point6.copyMirror(Point0,Vector3d(0, 1, 0));
Point49 = Point7.copyMirror(Point0,Vector3d(0, 1, 0));
Point50 = Point8.copyMirror(Point0,Vector3d(0, 1, 0));
Point51 = Point10.copyMirror(Point0,Vector3d(0, 1, 0));
Point52 = Point9.copyMirror(Point0,Vector3d(0, 1, 0));
Point53 = Point11.copyMirror(Point0,Vector3d(0, 1, 0));
Point54 = Point12.copyMirror(Point0,Vector3d(0, 1, 0));
Point55 = Point13.copyMirror(Point0,Vector3d(0, 1, 0));
Point56 = Point14.copyMirror(Point0,Vector3d(0, 1, 0));
Point57 = Point15.copyMirror(Point0,Vector3d(0, 1, 0));
Point58 = Point16.copyMirror(Point0,Vector3d(0, 1, 0));
Point59 = Point17.copyMirror(Point0,Vector3d(0, 1, 0));
Point60 = Point18.copyMirror(Point0,Vector3d(0, 1, 0));
Point61 = Point19.copyMirror(Point0,Vector3d(0, 1, 0));
Point62 = Point20.copyMirror(Point0,Vector3d(0, 1, 0));
Point63 = Point21.copyMirror(Point0,Vector3d(0, 1, 0));
Point64 = Point22.copyMirror(Point0,Vector3d(0, 1, 0));
Point65 = Point23.copyMirror(Point0,Vector3d(0, 1, 0));
Point66 = Point24.copyMirror(Point0,Vector3d(0, 1, 0));
Point67 = Point25.copyMirror(Point0,Vector3d(0, 1, 0));
Point68 = Point26.copyMirror(Point0,Vector3d(0, 1, 0));
Point69 = Point27.copyMirror(Point0,Vector3d(0, 1, 0));
Point70 = Point28.copyMirror(Point0,Vector3d(0, 1, 0));
Point71 = Point29.copyMirror(Point0,Vector3d(0, 1, 0));
Point72 = Point30.copyMirror(Point0,Vector3d(0, 1, 0));
Point73 = Point31.copyMirror(Point0,Vector3d(0, 1, 0));
Point74 = Point32.copyMirror(Point0,Vector3d(0, 1, 0));
Point75 = Point33.copyMirror(Point0,Vector3d(0, 1, 0));
Point76 = Point34.copyMirror(Point0,Vector3d(0, 1, 0));
Point77 = Point35.copyMirror(Point0,Vector3d(0, 1, 0));
Point78 = Point36.copyMirror(Point0,Vector3d(0, 1, 0));
Point79 = Point37.copyMirror(Point0,Vector3d(0, 1, 0));
Point80 = Point38.copyMirror(Point0,Vector3d(0, 1, 0));
Point81 = Point39.copyMirror(Point0,Vector3d(0, 1, 0));
Point82 = Point40.copyMirror(Point0,Vector3d(0, 1, 0));
Point83 = Point41.copyMirror(Point0,Vector3d(0, 1, 0));
Curve134 = Curve1.copyMirror(Point0,Vector3d(0, 1, 0));
Curve135 = Curve2.copyMirror(Point0,Vector3d(0, 1, 0));
Curve136 = Curve3.copyMirror(Point0,Vector3d(0, 1, 0));
Curve137 = Curve38.copyMirror(Point0,Vector3d(0, 1, 0));
Curve138 = Curve8.copyMirror(Point0,Vector3d(0, 1, 0));
Curve139 = Curve6.copyMirror(Point0,Vector3d(0, 1, 0));
Curve140 = Curve7.copyMirror(Point0,Vector3d(0, 1, 0));
Curve141 = Curve9.copyMirror(Point0,Vector3d(0, 1, 0));
Curve142 = Curve10.copyMirror(Point0,Vector3d(0, 1, 0));
Curve143 = Curve5.copyMirror(Point0,Vector3d(0, 1, 0));
Curve144 = Curve11.copyMirror(Point0,Vector3d(0, 1, 0));
Curve145 = Curve12.copyMirror(Point0,Vector3d(0, 1, 0));
Curve146 = Curve13.copyMirror(Point0,Vector3d(0, 1, 0));
Curve147 = Curve43.copyMirror(Point0,Vector3d(0, 1, 0));
Curve148 = Curve15.copyMirror(Point0,Vector3d(0, 1, 0));
Curve149 = Curve16.copyMirror(Point0,Vector3d(0, 1, 0));
Curve150 = Curve40.copyMirror(Point0,Vector3d(0, 1, 0));
Curve151 = Curve18.copyMirror(Point0,Vector3d(0, 1, 0));
Curve152 = Curve19.copyMirror(Point0,Vector3d(0, 1, 0));
Curve153 = Curve20.copyMirror(Point0,Vector3d(0, 1, 0));
Curve154 = Curve21.copyMirror(Point0,Vector3d(0, 1, 0));
Curve155 = Curve22.copyMirror(Point0,Vector3d(0, 1, 0));
Curve156 = Curve23.copyMirror(Point0,Vector3d(0, 1, 0));
Curve157 = Curve24.copyMirror(Point0,Vector3d(0, 1, 0));
Curve158 = Curve25.copyMirror(Point0,Vector3d(0, 1, 0));
Curve159 = Curve26.copyMirror(Point0,Vector3d(0, 1, 0));
Curve160 = Curve27.copyMirror(Point0,Vector3d(0, 1, 0));
Curve161 = Curve28.copyMirror(Point0,Vector3d(0, 1, 0));
Curve162 = Curve29.copyMirror(Point0,Vector3d(0, 1, 0));
Curve163 = Curve30.copyMirror(Point0,Vector3d(0, 1, 0));
Curve164 = Curve31.copyMirror(Point0,Vector3d(0, 1, 0));
Curve165 = Curve32.copyMirror(Point0,Vector3d(0, 1, 0));
Curve166 = Curve33.copyMirror(Point0,Vector3d(0, 1, 0));
Curve167 = Curve34.copyMirror(Point0,Vector3d(0, 1, 0));
Curve168 = Curve35.copyMirror(Point0,Vector3d(0, 1, 0));
Curve169 = Curve37.copyMirror(Point0,Vector3d(0, 1, 0));
Curve170 = Curve36.copyMirror(Point0,Vector3d(0, 1, 0));
Curve171 = Curve39.copyMirror(Point0,Vector3d(0, 1, 0));
Curve172 = Curve42.copyMirror(Point0,Vector3d(0, 1, 0));
Curve173 = Curve41.copyMirror(Point0,Vector3d(0, 1, 0));
Curve174 = Curve44.copyMirror(Point0,Vector3d(0, 1, 0));
Curve175 = Curve45.copyMirror(Point0,Vector3d(0, 1, 0));
Curve176 = Curve46.copyMirror(Point0,Vector3d(0, 1, 0));
Curve177 = Curve47.copyMirror(Point0,Vector3d(0, 1, 0));
Curve178 = Curve48.copyMirror(Point0,Vector3d(0, 1, 0));
Curve179 = Curve49.copyMirror(Point0,Vector3d(0, 1, 0));
Curve180 = Curve50.copyMirror(Point0,Vector3d(0, 1, 0));
Curve181 = Curve51.copyMirror(Point0,Vector3d(0, 1, 0));
Curve182 = Curve52.copyMirror(Point0,Vector3d(0, 1, 0));
Curve183 = Curve53.copyMirror(Point0,Vector3d(0, 1, 0));
Curve184 = Curve54.copyMirror(Point0,Vector3d(0, 1, 0));
Curve185 = Curve55.copyMirror(Point0,Vector3d(0, 1, 0));
Curve186 = Curve56.copyMirror(Point0,Vector3d(0, 1, 0));
Curve187 = Curve57.copyMirror(Point0,Vector3d(0, 1, 0));
Curve188 = Curve58.copyMirror(Point0,Vector3d(0, 1, 0));
Curve189 = Curve59.copyMirror(Point0,Vector3d(0, 1, 0));
Curve190 = Curve60.copyMirror(Point0,Vector3d(0, 1, 0));
Curve191 = Curve61.copyMirror(Point0,Vector3d(0, 1, 0));
Curve192 = Curve62.copyMirror(Point0,Vector3d(0, 1, 0));
Curve193 = Curve81.copyMirror(Point0,Vector3d(0, 1, 0));
Curve194 = Curve82.copyMirror(Point0,Vector3d(0, 1, 0));
Curve195 = Curve83.copyMirror(Point0,Vector3d(0, 1, 0));
Curve196 = Curve84.copyMirror(Point0,Vector3d(0, 1, 0));
Curve197 = Curve85.copyMirror(Point0,Vector3d(0, 1, 0));
Curve198 = Curve86.copyMirror(Point0,Vector3d(0, 1, 0));
Curve199 = Curve87.copyMirror(Point0,Vector3d(0, 1, 0));
Curve200 = Curve88.copyMirror(Point0,Vector3d(0, 1, 0));
Curve201 = Curve89.copyMirror(Point0,Vector3d(0, 1, 0));
Curve202 = Curve90.copyMirror(Point0,Vector3d(0, 1, 0));
Curve203 = Curve91.copyMirror(Point0,Vector3d(0, 1, 0));
Curve204 = Curve92.copyMirror(Point0,Vector3d(0, 1, 0));
Curve205 = Curve93.copyMirror(Point0,Vector3d(0, 1, 0));
Curve206 = Curve94.copyMirror(Point0,Vector3d(0, 1, 0));
Curve207 = Curve95.copyMirror(Point0,Vector3d(0, 1, 0));
Curve208 = Curve96.copyMirror(Point0,Vector3d(0, 1, 0));
Curve209 = Curve97.copyMirror(Point0,Vector3d(0, 1, 0));
Curve210 = Curve98.copyMirror(Point0,Vector3d(0, 1, 0));
Curve211 = Curve99.copyMirror(Point0,Vector3d(0, 1, 0));
Curve212 = Curve100.copyMirror(Point0,Vector3d(0, 1, 0));
Curve213 = Curve101.copyMirror(Point0,Vector3d(0, 1, 0));
Curve214 = Curve102.copyMirror(Point0,Vector3d(0, 1, 0));
Curve215 = Curve103.copyMirror(Point0,Vector3d(0, 1, 0));
Curve216 = Curve104.copyMirror(Point0,Vector3d(0, 1, 0));
Curve217 = Curve105.copyMirror(Point0,Vector3d(0, 1, 0));
Curve218 = Curve106.copyMirror(Point0,Vector3d(0, 1, 0));
Curve219 = Curve107.copyMirror(Point0,Vector3d(0, 1, 0));
Curve220 = Curve108.copyMirror(Point0,Vector3d(0, 1, 0));
Curve221 = Curve109.copyMirror(Point0,Vector3d(0, 1, 0));
Curve222 = Curve110.copyMirror(Point0,Vector3d(0, 1, 0));
Curve223 = Curve111.copyMirror(Point0,Vector3d(0, 1, 0));
Curve224 = Curve112.copyMirror(Point0,Vector3d(0, 1, 0));
Curve225 = Curve113.copyMirror(Point0,Vector3d(0, 1, 0));
Curve226 = Curve114.copyMirror(Point0,Vector3d(0, 1, 0));
Curve227 = Curve115.copyMirror(Point0,Vector3d(0, 1, 0));
Curve228 = Curve116.copyMirror(Point0,Vector3d(0, 1, 0));
Curve229 = Curve117.copyMirror(Point0,Vector3d(0, 1, 0));
Curve230 = Curve118.copyMirror(Point0,Vector3d(0, 1, 0));
Curve231 = Curve119.copyMirror(Point0,Vector3d(0, 1, 0));
Curve232 = Curve120.copyMirror(Point0,Vector3d(0, 1, 0));
Curve233 = Curve121.copyMirror(Point0,Vector3d(0, 1, 0));
Curve234 = Curve122.copyMirror(Point0,Vector3d(0, 1, 0));
Curve235 = Curve123.copyMirror(Point0,Vector3d(0, 1, 0));
Curve236 = Curve124.copyMirror(Point0,Vector3d(0, 1, 0));
Curve237 = Curve125.copyMirror(Point0,Vector3d(0, 1, 0));
Curve238 = Curve126.copyMirror(Point0,Vector3d(0, 1, 0));
Curve239 = Curve127.copyMirror(Point0,Vector3d(0, 1, 0));
Curve240 = Curve128.copyMirror(Point0,Vector3d(0, 1, 0));
Curve241 = Curve129.copyMirror(Point0,Vector3d(0, 1, 0));
Curve242 = Curve131.copyMirror(Point0,Vector3d(0, 1, 0));
Curve243 = Curve130.copyMirror(Point0,Vector3d(0, 1, 0));
Curve244 = Curve132.copyMirror(Point0,Vector3d(0, 1, 0));
Curve245 = Curve133.copyMirror(Point0,Vector3d(0, 1, 0));
Delete(Pl116);
Delete(Pl115);
Delete(Pl145);
Pl151 = SweepCurve(Curve141, Curve160);
Pl153 = SweepCurve(Curve168, Curve160);
Pl157 = SweepCurve(Curve169, Curve160);
Pl153.flipNormal();
Pl151.flipNormal();

Pl32.join(Pl139);
Pl32.simplifyTopology();
Pl31.join(Pl138);
Pl31.simplifyTopology();
Pl2.join(Pl112);
Pl2.simplifyTopology();


wetted_surface.add(Pl111);
wetted_surface.add(Pl113);
wetted_surface.add(Pl114);
wetted_surface.add(Pl117);
wetted_surface.add(Pl118);
wetted_surface.add(Pl119);
wetted_surface.add(Pl120);
wetted_surface.add(Pl121);
wetted_surface.add(Pl122);
wetted_surface.add(Pl123);
wetted_surface.add(Pl124);
wetted_surface.add(Pl125);
wetted_surface.add(Pl126);
wetted_surface.add(Pl127);
wetted_surface.add(Pl128);
wetted_surface.add(Pl129);
wetted_surface.add(Pl130);
wetted_surface.add(Pl131);
wetted_surface.add(Pl132);
wetted_surface.add(Pl133);
wetted_surface.add(Pl134);
wetted_surface.add(Pl135);
wetted_surface.add(Pl136);
wetted_surface.add(Pl137);
wetted_surface.add(Pl151);
wetted_surface.add(Pl153);
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


Validate(Pl_fs1.primitivePartCount == 20);
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
Validate(Pl_fs17);
Validate(Pl_fs18);
Validate(Pl_fs19);
Validate(Pl_fs20);
Validate(Pl_fs21);

Delete(Pl_fs20);
Delete(Pl_fs16);
Delete(Pl_fs7);
Delete(Pl_fs12);
Delete(Pl_fs3);
Delete(Pl_fs2);
Delete(Pl_fs4);
Delete(Pl_fs19);
Delete(Pl_fs13);
Delete(Pl_fs15);
Delete(Pl_fs11);
Delete(Pl_fs17);
Delete(Pl_fs21);
Delete(Pl_fs18);
Delete(Pl_fs8);
Delete(Pl_fs5);
Delete(Pl_fs6);
Delete(Pl_fs9);
Delete(Pl_fs10);

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

Pl_fs14.front.wetSurface = outer_wet;
Free_surface = Set();
Free_surface.add(Pl_fs14);
Meshing_panels.step(1).subset = Free_surface;

LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

Meshing_panels.step(1).subset = Free_surface;

Pl_fs14.meshDensity = FS_density;
GenieRules.Meshing.superElementType = 73;
Meshing_panels.execute();
ExportMeshFem().DoExport("T73.FEM");


/*
Generate mesh of internal lids for removing irregular frequencies
 */

DefaultName(typeFlatPlate,"Pl_int_lid",1,"");
DefaultName(typePoint,"Point_int_lid",1,"");
DefaultName(typeGuideCurve,"Curve_int_lid",1,"");

Lid_density = MeshDensity(mesh_lid);

Pl_int_lid1 = CoverCurves(Curve50,Curve51,Curve55);
Pl_int_lid2 = CoverCurves(Curve49,Curve53,Curve54);
Pl_int_lid3 = CoverCurves(Curve46,Curve47,Curve48,Curve52);
Pl_int_lid3.flipNormal();



Internal_lid = Set();
Internal_lid.add(Pl_int_lid1);
Internal_lid.add(Pl_int_lid2);
Internal_lid.add(Pl_int_lid3);

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
Pl_int_lid3.permeability = Surf_lid;
cm = CompartmentManager();

LC2 = LoadCase();
LC2.compartment(cm.compartment(Point(-0.3729626562 m,0.3403159927 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC3 = LoadCase();
LC3.compartment(cm.compartment(Point(-1.405334509 m,3.711816231 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC4 = LoadCase();
LC4.compartment(cm.compartment(Point(-2.011992465 m,9.734873243 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC5 = LoadCase();
LC5.compartment(cm.compartment(Point(-2.527153525 m,-3.384004496 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC6 = LoadCase();
LC6.compartment(cm.compartment(Point(-3.23127603 m,-2.71058078 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC7 = LoadCase();
LC7.compartment(cm.compartment(Point(-3.771766569 m,1.347247116 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC8 = LoadCase();
LC8.compartment(cm.compartment(Point(-6.65175327 m,-23.00860423 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC9 = LoadCase();
LC9.compartment(cm.compartment(Point(-7.396609698 m,6.56130387 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC10 = LoadCase();
LC10.compartment(cm.compartment(Point(-12.10092976 m,27.20942516 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC11 = LoadCase();
LC11.compartment(cm.compartment(Point(-12.21368271 m,-19.79481517 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC12 = LoadCase();
LC12.compartment(cm.compartment(Point(-14.42159194 m,-33.48039285 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC13 = LoadCase();
LC13.compartment(cm.compartment(Point(-17.63742339 m,24.29891342 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC14 = LoadCase();
LC14.compartment(cm.compartment(Point(-19.95406542 m,43.42386366 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC15 = LoadCase();
LC15.compartment(cm.compartment(Point(-20.44863582 m,-47.7401885 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC16 = LoadCase();
LC16.compartment(cm.compartment(Point(-20.70641795 m,-43.04563902 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC17 = LoadCase();
LC17.compartment(cm.compartment(Point(-21.42163134 m,-26.12306043 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC18 = LoadCase();
LC18.compartment(cm.compartment(Point(-21.60558948 m,41.60984172 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC19 = LoadCase();
LC19.compartment(cm.compartment(Point(-22.54986804 m,-49.20960506 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC20 = LoadCase();
LC20.compartment(cm.compartment(Point(-24.10382437 m,49.9853967 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC21 = LoadCase();
LC21.compartment(cm.compartment(Point(-25.26729287 m,-45.48281939 m,10 m))).globalIntensity = DummyHydroPressure();
LC22 = LoadCase();
LC22.compartment(cm.compartment(Point(-25.41314258 m,-45.32297825 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC23 = LoadCase();
LC23.compartment(cm.compartment(Point(-25.43625797 m,-45.27820032 m,-8 m))).globalIntensity = DummyHydroPressure();
LC24 = LoadCase();
LC24.compartment(cm.compartment(Point(-26.3811636 m,45.27867206 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC25 = LoadCase();
LC25.compartment(cm.compartment(Point(-26.10107061 m,45.01455811 m,-8 m))).globalIntensity = DummyHydroPressure();
LC26 = LoadCase();
LC26.compartment(cm.compartment(Point(-26.54100474 m,45.42452178 m,10 m))).globalIntensity = DummyHydroPressure();
LC27 = LoadCase();
LC27.compartment(cm.compartment(Point(-27.01032115 m,40.01564277 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC28 = LoadCase();
LC28.compartment(cm.compartment(Point(-27.26795098 m,-38.89588007 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC29 = LoadCase();
LC29.compartment(cm.compartment(Point(-28.11098136 m,-49.79286786 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC30 = LoadCase();
LC30.compartment(cm.compartment(Point(-28.79837385 m,50.24317882 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC31 = LoadCase();
LC31.compartment(cm.compartment(Point(-29.08197292 m,-40.54740412 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC32 = LoadCase();
LC32.compartment(cm.compartment(Point(-29.52453926 m,40.74391567 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC33 = LoadCase();
LC33.compartment(cm.compartment(Point(-29.94789897 m,-48.4663539 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC34 = LoadCase();
LC34.compartment(cm.compartment(Point(-30.26779041 m,48.1419466 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC35 = LoadCase();
LC35.compartment(cm.compartment(Point(-30.67617187 m,-45.95213579 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC36 = LoadCase();
LC36.compartment(cm.compartment(Point(-30.85105321 m,42.58083328 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC37 = LoadCase();
LC37.compartment(cm.compartment(Point(0.3515555134 m,-0.389540646 m,-8 m))).globalIntensity = DummyHydroPressure();
LC38 = LoadCase();
LC38.compartment(cm.compartment(Point(0.4861657038 m,-0.5328037946 m,10 m))).globalIntensity = DummyHydroPressure();
LC39 = LoadCase();
LC39.compartment(cm.compartment(Point(0.1018733216 m,4.787375642 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC40 = LoadCase();
LC40.compartment(cm.compartment(Point(1.347247116 m,-3.771766569 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC41 = LoadCase();
LC41.compartment(cm.compartment(Point(3.711816231 m,-1.405334509 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC42 = LoadCase();
LC42.compartment(cm.compartment(Point(4.783166336 m,0.09967199429 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC43 = LoadCase();
LC43.compartment(cm.compartment(Point(5.024478968 m,4.895893591 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC44 = LoadCase();
LC44.compartment(cm.compartment(Point(5.623991435 m,-1.605280239 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC45 = LoadCase();
LC45.compartment(cm.compartment(Point(26.65725182 m,-0.8212507598 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC46 = LoadCase();
LC46.compartment(cm.compartment(Point(26.65725182 m,5.42874924 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC47 = LoadCase();
LC47.compartment(cm.compartment(Point(46.23765638 m,-0.9832154886 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC48 = LoadCase();
LC48.compartment(cm.compartment(Point(46.86535074 m,0.144687735 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC49 = LoadCase();
LC49.compartment(cm.compartment(Point(48.22845033 m,-4.171108838 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC50 = LoadCase();
LC50.compartment(cm.compartment(Point(49.91265694 m,4.770180465 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC51 = LoadCase();
LC51.compartment(cm.compartment(Point(51.94442224 m,-4.923027517 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC52 = LoadCase();
LC52.compartment(cm.compartment(Point(52.21185742 m,-0.5061636049 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC53 = LoadCase();
LC53.compartment(cm.compartment(Point(52.35770713 m,-0.6660047432 m,10 m))).globalIntensity = DummyHydroPressure();
LC54 = LoadCase();
LC54.compartment(cm.compartment(Point(52.50707062 m,5.400329762 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC55 = LoadCase();
LC55.compartment(cm.compartment(Point(52.89853886 m,-1.417510011 m,-8 m))).globalIntensity = DummyHydroPressure();
LC56 = LoadCase();
LC56.compartment(cm.compartment(Point(56.52018047 m,-1.837343063 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC57 = LoadCase();
LC57.compartment(cm.compartment(Point(57.15032976 m,0.7570706222 m,-1.5 m))).globalIntensity = DummyHydroPressure();

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
DefaultName(typePoint,"Point_Morison",1,"");

Column1 = PipeSection(d2, t1);

Column2 = PipeSection(d1, t1);

Pontoon = BoxSection(h_pontoon, w_pontoon, t1, t1, t1);

Column1.setDefault();

Point_Morison1 = Point47.copyTranslate(Vector3d(0 m,0 m,draft-fairlead_depth));
Point_Morison2 = Point47.copyTranslate(Vector3d(0 m,0 m,draft+freeboard));
Bm1 = StraightBeam(Point47, Point_Morison1);
Bm2 = StraightBeam(Point_Morison1, Point_Morison2);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Bm1, "Bm3");
MyModelTransformerMap.Add(Bm2, "Bm4");
ModelTransformer(MyModelTransformerMap).copyRotate(Point42, Vector3d(0, 0, 1), 120, 2);


Point_Morison3 = Point0.copyTranslate(Vector3d(0 m,0 m,draft-fairlead_depth));
Point_Morison4 = Point0.copyTranslate(Vector3d(0 m,0 m,draft+freeboard));
Column2.setDefault();
Bm7 = StraightBeam(Point0, Point_Morison3);
Bm8 = StraightBeam(Point_Morison3, Point_Morison4);


Point_Morison5 = Point(d1/2,0 m,-draft);Point_Morison6 = Point46.copyTranslate(Vector3d(0 m,0 m, h_pontoon/2));
Point_Morison7 = Point_Morison5.copyTranslate(Vector3d(0 m,0 m, h_pontoon/2));
Pontoon.setDefault();
Bm9 = StraightBeam(Point_Morison6, Point_Morison7);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Bm9, "Bm10");
ModelTransformer(MyModelTransformerMap).copyRotate(Point42, Vector3d(0, 0, 1), 120, 2);

Morison = Set();Morison.add(Bm1);
Morison.add(Bm10);
Morison.add(Bm11);
Morison.add(Bm2);
Morison.add(Bm3);
Morison.add(Bm4);
Morison.add(Bm5);
Morison.add(Bm6);
Morison.add(Bm7);
Morison.add(Bm8);
Morison.add(Bm9);

Delete(Md_def);
Delete(FS_density);
Delete(Lid_density);

GenieRules.Meshing.autoSimplifyTopology = false;
GenieRules.Meshing.autoSplitPeriodicGeometry = false;
GenieRules.Meshing.superElementType = 6;

Meshing_panels.step(1).subset = Morison;
Meshing_panels.step(1).excludeIncludeMeshSubsetOption = anIncludeMeshSubset;
Meshing_panels.execute();
ExportMeshFem().DoExport("T6.FEM");