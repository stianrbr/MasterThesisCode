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



mesh_size = 0.50 m;  // Mesh size
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
Point1 = Point(-cc_distance, 0m, -draft);
Point2 = Point(-cc_distance, d2/2, -draft);
Point3 = Point(-cc_distance+d2/2, 0m, -draft);
Point4 = Point(-cc_distance-d2/2, 0m, -draft);
Curve1 = CreateCircularArcFromThreePoints(Point3, Point2, Point4);
Curve2 = CreateLineTwoPoints(Point4, Point1);
Curve3 = CreateLineTwoPoints(Point1, Point3);
x_trans = cc_distance*(1+Math.cos(60*Math.PI/180));
y_trans = cc_distance*Math.sin(60*Math.PI/180);
Point5 = Point1.copyTranslate(Vector3d(x_trans, y_trans, 0));
Point6 = Point2.copyTranslate(Vector3d(x_trans, y_trans, 0));
Point7 = Point3.copyTranslate(Vector3d(x_trans, y_trans, 0));
Point8 = Point4.copyTranslate(Vector3d(x_trans, y_trans, 0));
Curve4 = CreateCircleFromThreePoints(Point7, Point6, Point8);
temp = w_pontoon/2*Math.tan(120*Math.PI/180);
Point9 = Point(-temp, w_pontoon/2, -draft);
Point9 = Point(temp, w_pontoon/2, -draft);
Curve5 = CreateLineTwoPoints(Point2, Point9);
Curve6 = CreateLineTwoPoints(Point0, Point5);
Curve7 = CreateLineTwoPoints(Curve1.end1, Curve6.end1);
temp = d2/2*Math.cos(Math.asin(w_pontoon/d2));
Point9 = Point(-temp, w_pontoon/2, -draft);
r = w_pontoon/(2*Math.sin(120*Math.PI/180));

Point9 = Point(r*Math.cos(120*Math.PI/180), w_pontoon/2, -draft);
Delete(Curve5);
Curve8 = CreateLineTwoPoints(Point2, Point9);
x = cc_distance*Math.cos(60*Math.PI/180)+d2/2*Math.cos(150*Math.PI/180);
y= cc_distance*Math.sin(60*Math.PI/180)+d2/2*Math.sin(150*Math.PI/180);

Point10 = Point(x,y,-draft);
Curve9 = CreateLineTwoPoints(Curve8.end2, Point10);
x = cc_distance*Math.cos(60*Math.PI/180)+d2/2*Math.cos(-30*Math.PI/180);
y= cc_distance*Math.sin(60*Math.PI/180)+d2/2*Math.sin(-30*Math.PI/180);

Point11 = Point(x,y,-draft);
x = w_pontoon/(2*Math.tan(60*Math.PI/180))+Math.sqrt(r*r - (w_pontoon/2)*(w_pontoon/2));

Point12 = Point(x, 0, -draft);
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
Point26 = Point(d1/2, 0m, -draft+h_pontoon);
Point27 = Point(-d1/2, 0m, -draft+h_pontoon);
Point28 = Point(0m, d1/2, -draft+h_pontoon);
Curve30 = CreateCircularArcFromThreePoints(Point26, Point28, Point27);
Curve31 = Curve1.divideAt(0.5);
Curve32 = Curve11.divideAt(0.5);
Pl1 = CoverCurves(Curve21,Curve22,Curve31,Curve32);
Pl1.flipNormal();
Pl2 = CoverCurves(Curve1,Curve2,Curve3,Curve31);
Pl2.flipNormal();
Curve33 = CreateLineTwoPoints(Point0, Curve9.end1);
Curve34 = CreateLineTwoPoints(Point13, Curve15.end2);
Pl3 = CoverCurves(Curve1,Curve8,Curve7,Curve33);
Pl3.flipNormal();
Pl4 = CoverCurves(Curve8,Curve15,Curve22,Curve24);
Pl5 = CoverCurves(Curve9,Curve18,Curve24,Curve27);
Curve35 = Curve4.divideAt(0.4166666691);
Curve36 = Curve35.divideAt(0.8571428582);
tmpArrayOfCurves = JoinMultipleCurves(Array(Curve4, Curve36), true);
Rename(tmpArrayOfCurves[0], "Curve36");
Pl6 = SweepCurve(Curve36, Curve28);
Pl6.flipNormal();
Curve37 = Curve35.divideAt(0.4999999971);
Curve38 = Curve6.divideAt(0.8792270531);
Pl7 = CoverCurves(Curve6,Curve9,Curve33,Curve35);
Pl8 = CoverCurves(Curve6,Curve10,Curve5,Curve37);
Pl9 = CoverCurves(Curve35,Curve37,Curve36);
Pl10 = CoverCurves(Curve10,Curve19,Curve26,Curve28);
Pl11 = CoverCurves(Curve11,Curve15,Curve17,Curve34);
Curve39 = Curve14.divideAt(0.4166666676);
Curve40 = Curve39.divideAt(0.8571428586);
tmpArrayOfCurves = JoinMultipleCurves(Array(Curve14, Curve40), true);
Rename(tmpArrayOfCurves[0], "Curve40");
Curve41 = Curve39.divideAt(0.4999999986);
Curve42 = Curve16.divideAt(0.8792270531);
Pl12 = CoverCurves(Curve16,Curve18,Curve34,Curve39);
Pl12.flipNormal();
Pl13 = CoverCurves(Curve16,Curve19,Curve20,Curve41);
Pl13.flipNormal();
Curve43 = Curve17.divideAt(0.8901098901);
Curve44 = Curve20.divideAt(0.692820323);
Point29 = Point13.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point30 = Point14.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point31 = Point15.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point32 = Point16.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point33 = Point17.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point34 = Point18.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point35 = Point19.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point36 = Point20.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point37 = Point21.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point38 = Point23.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point39 = Point25.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point40 = Point26.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point41 = Point27.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Point42 = Point28.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve45 = Curve11.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve46 = Curve12.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve47 = Curve13.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve48 = Curve20.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve49 = Curve30.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve50 = Curve32.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve51 = Curve39.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve52 = Curve41.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve53 = Curve40.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve54 = Curve43.copyTranslate(Vector3d(0 m,0 m,draft-h_pontoon));
Curve55 = CreateLineTwoPoints(Curve12.end1, Curve46.end1);
Curve56 = CreateLineTwoPoints(Curve17.end2, Curve49.end2);
Curve57 = CreateLineTwoPoints(Curve40.end1, Curve51.end1);
Delete(Pl11);
Delete(Pl12);
Delete(Pl13);
Curve58 = Curve34.divideAt(0.692820323);
Curve59 = Curve30.divideAt(0.6666666676);
Curve60 = Curve16.divideAt(0.1098901099);
Curve61 = Curve30.divideAt(0.5000000111);
Pl14 = CoverCurves(Curve11,Curve15,Curve17,Curve58,Curve59);
Pl15 = CoverCurves(Curve18,Curve39,Curve58,Curve60,Curve61);
Pl15.flipNormal();
Pl16 = CoverCurves(Curve19,Curve30,Curve41,Curve44,Curve60);
Pl16.flipNormal();
Pl17 = SweepCurve(Curve11, Curve55);
Pl18 = SweepCurve(Curve32, Curve55);
Pl19 = SweepCurve(Curve30, Curve56);
Pl20 = SweepCurve(Curve59, Curve56);
Pl21 = SweepCurve(Curve61, Curve56);
Pl22 = SweepCurve(Curve39, Curve57);
Pl23 = SweepCurve(Curve41, Curve57);
Pl24 = SweepCurve(Curve40, Curve57);


Curve65 = CreateLineTwoPoints(Point(22.75 m,39.40415587 m,-13 m), Curve60.end1);
Curve62 = CreateLineTwoPoints(Curve60.end1, Curve18.project(Curve60.end1));
Curve63 = CreateLineTwoPoints(Point(-5 m,0 m,-13 m), Curve15.project(Point(-5 m,0 m,-13 m)));
Curve64 = CreateLineTwoPoints(Curve60.end1, Curve19.project(Curve60.end1));
Curve66 = CreateLineTwoPoints(Curve64.end2, Curve10.project(Curve64.end2));
Curve67 = CreateLineTwoPoints(Curve66.end2, Curve6.project(Curve66.end2));
Curve68 = CreateLineTwoPoints(Curve67.end2, Curve9.project(Curve67.end2));
Curve69 = CreateLineTwoPoints(Curve68.end2, Curve18.project(Curve68.end2));
Curve70 = CreateLineTwoPoints(Curve63.end2, Curve8.project(Curve63.end2));
Curve71 = CreateLineTwoPoints(Curve70.end2, Curve7.project(Curve70.end2));
Curve72 = CreateLineTwoPoints(Point(-45.5 m,0 m,-13 m), Curve15.project(Point(-45.5 m,0 m,-13 m)));
Curve73 = CreateLineTwoPoints(Curve72.end2, Curve8.project(Curve72.end2));
Curve74 = CreateLineTwoPoints(Curve7.end1, Curve8.project(Curve7.end1));
Curve75 = CreateLineTwoPoints(Curve41.end1, Curve18.project(Curve41.end1));
Curve76 = CreateLineTwoPoints(Curve75.intersect(Curve18), Curve9.project(Curve75.intersect(Curve18)));
Curve77 = CreateLineTwoPoints(Curve76.end2, Point(22.75000005 m,39.40415584 m,-20 m));
Curve78 = CreateLineTwoPoints(Curve6.end2, Curve10.project(Curve6.end2));
Curve79 = CreateLineTwoPoints(Curve78.end2, Curve19.project(Curve78.end2));
Curve80 = CreateLineTwoPoints(Curve79.end2, Point(22.75 m,39.40415587 m,-13 m));
FEdge1 = FeatureEdge(Curve24);
FEdge2 = FeatureEdge(Curve33);
FEdge3 = FeatureEdge(Curve58);
FEdge4 = FeatureEdge(Curve62);
FEdge5 = FeatureEdge(Curve63);
FEdge6 = FeatureEdge(Curve64);
FEdge7 = FeatureEdge(Curve66);
FEdge8 = FeatureEdge(Curve67);
FEdge9 = FeatureEdge(Curve68);
FEdge10 = FeatureEdge(Curve69);
FEdge11 = FeatureEdge(Curve70);
FEdge12 = FeatureEdge(Curve71);
FEdge13 = FeatureEdge(Curve72);
FEdge14 = FeatureEdge(Curve73);
FEdge15 = FeatureEdge(Curve74);
FEdge16 = FeatureEdge(Curve75);
FEdge17 = FeatureEdge(Curve76);
FEdge18 = FeatureEdge(Curve77);
FEdge19 = FeatureEdge(Curve78);
FEdge20 = FeatureEdge(Curve79);
FEdge21 = FeatureEdge(Curve80);

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
wetted_surface.front.wetSurface = outer_wet;


GenieRules.Meshing.faceMeshStrategy = PatchSurfQuadMesher;
Md_def.setDefault();
LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

Meshing_panels.execute();
ExportMeshFem().DoExport("T1.FEM");


/*
Modeling of the structure above MWL, still applying symmetry
 */

Curve62 = Curve45.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve63 = Curve46.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve64 = Curve47.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve65 = Curve48.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve66 = Curve49.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve67 = Curve50.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve68 = Curve51.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve69 = Curve52.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve70 = Curve53.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve71 = Curve54.copyTranslate(Vector3d(0 m,0 m,freeboard));
Curve72 = CreateLineTwoPoints(Point(-58 m,0 m,0 m), Point(-58 m,0 m,15 m));
Pl25 = SweepCurve(Curve45, Curve72);
Pl26 = SweepCurve(Curve49, Curve72);
Pl27 = SweepCurve(Curve50, Curve72);
Pl28 = SweepCurve(Curve51, Curve72);
Pl29 = SweepCurve(Curve52, Curve72);
Pl30 = SweepCurve(Curve53, Curve72);
Pl25.flipNormal();
Pl26.flipNormal();
Pl27.flipNormal();
Pl30.flipNormal();
Pl28.flipNormal();
Pl29.flipNormal();
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

Pl31 = CoverCurves(Curve62,Curve63,Curve64,Curve67);
Pl32 = CoverCurves(Curve65,Curve66,Curve71);
Pl33 = CoverCurves(Curve68,Curve69,Curve70);
Pl34 = CoverCurves(Curve1,Curve11,Curve22,Curve23);
Pl34.flipNormal();
Pl35 = CoverCurves(Curve24,Curve25,Curve33,Curve34,Curve58);
Pl36 = CoverCurves(Curve20,Curve30,Curve43,Curve59,Curve61);
Pl37 = CoverCurves(Curve5,Curve20,Curve25,Curve26,Curve44);
Pl38 = SweepCurve(Curve35, Curve27);
Pl39 = SweepCurve(Curve37, Curve28);


Curve91 = Curve45.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve92 = Curve46.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve93 = Curve47.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve94 = Curve48.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve95 = Curve49.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve96 = Curve50.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve97 = Curve51.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve98 = Curve52.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve99 = Curve53.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve100 = Curve54.copyTranslate(Vector3d(0 m,0 m,deck_collition_high));
Curve101 = Curve45.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve102 = Curve46.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve103 = Curve47.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve104 = Curve48.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve105 = Curve49.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve106 = Curve50.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve107 = Curve51.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve108 = Curve52.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve109 = Curve53.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Curve110 = Curve54.copyTranslate(Vector3d(0 m,0 m,-deck_collition_low));
Pl40 = CoverCurves(Curve101,Curve102,Curve103,Curve106);
Pl41 = CoverCurves(Curve91,Curve92,Curve93,Curve96);
Pl42 = CoverCurves(Curve94,Curve95,Curve100);
Pl43 = CoverCurves(Curve104,Curve105,Curve110);
Pl44 = CoverCurves(Curve97,Curve98,Curve99);
Pl45 = CoverCurves(Curve107,Curve108,Curve109);
Curve111 = Curve103.divideAt(0.76);
Curve112 = Curve111.copyRotate(Curve102.curvePoint(2),Vector3d(0, 0, 1),90*1);
Curve113 = Curve111.copyRotate(Curve102.curvePoint(2),Vector3d(0, 0, 1),90*2);
Curve114 = CreateCircularArcFromThreePoints(Curve113.end1, Curve112.project(Curve113.end1), Curve103.end2);
Curve115 = CreateLineTwoPoints(Curve103.end2, Curve93.project(Curve103.end2));
Pl46 = SweepCurve(Curve111, Curve115);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl46, "Pl47");
ModelTransformer(MyModelTransformerMap).copyRotate(Curve102.curvePoint(2), Vector3d(0, 0, 1), 45, 4);
Pl51 = SweepCurve(Curve114, Curve115);
Curve116 = Curve104.divideAt(0.7);
Curve117 = Curve116.copyRotate(Curve104.curvePoint(1),Vector3d(0, 0, 1),90*1);
Curve118 = Curve116.copyRotate(Curve104.curvePoint(1),Vector3d(0, 0, 1),90*2);
Curve119 = CreateCircularArcFromThreePoints(Curve118.end1, Curve117.project(Curve118.end1), Curve104.end2);
Curve120 = CreateLineTwoPoints(Curve116.end1, Curve94.project(Curve116.end1));
Pl52 = SweepCurve(Curve116, Curve120);
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl52, "Pl53");
ModelTransformer(MyModelTransformerMap).copyRotate(Curve104.curvePoint(1), Vector3d(0, 0, 1), 45, 4);
Pl57 = SweepCurve(Curve119, Curve120);
Curve121 = CreateLineTwoPoints(Point(22.75000005 m,39.40415584 m,-3 m), Point(28.99999995 m,50.22947345 m,-3 m));
Curve122 = Curve121.divideAt(0.12);
Curve123 = Curve121.copyRotate(Point34,Vector3d(0, 0, 1),90*1);
Curve124 = Curve121.copyRotate(Point34,Vector3d(0, 0, 1),90*2);
Curve125 = CreateCircleFromThreePoints(Curve121.end2, Curve123.project(Curve121.end2), Curve124.project(Curve123.project(Curve121.end2)));
Curve126 = CreateLineTwoPoints(Curve121.end1, Curve97.project(Curve121.end1));
Pl58 = SweepCurve(Curve121, Curve126);
Curve127 = CreateLineTwoPoints(Curve121.end2, Point(23.50000004 m,40.70319395 m,5 m));
MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl58, "Pl59");
ModelTransformer(MyModelTransformerMap).copyRotate(Point34, Vector3d(0, 0, 1), 45, 7);
Pl66 = SweepCurve(Curve125, Curve127);


Pl134 = SweepCurve(Curve6, Curve25);
Pl135 = SweepCurve(Curve7, Curve25);
Curve243 = Curve17.divideAt(0.5);
Curve244 = Curve60.divideAt(0.5);
Curve252 = CreateLineTwoPoints(Curve243.end1, Curve15.project(Curve243.end1));
Curve253 = CreateLineTwoPoints(Curve17.end2, Curve7.project(Curve17.end2));
Pl67 = SweepCurve(Curve252, Curve253);
Curve254 = CreateLineTwoPoints(Curve244.end1, Curve18.project(Curve244.end1));
Curve255 = CreateLineTwoPoints(Curve60.end2, Curve19.project(Curve60.end2));
Curve256 = CreateLineTwoPoints(Curve254.end1, Curve6.project(Curve254.end1));
Pl68 = SweepCurve(Curve254, Curve256);
Pl69 = SweepCurve(Curve255, Curve256);

/*
Mirror the structure and export full mesh
 */

MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(FEdge1, "FEdge22");
MyModelTransformerMap.Add(FEdge2, "FEdge23");
MyModelTransformerMap.Add(FEdge3, "FEdge24");
MyModelTransformerMap.Add(FEdge4, "FEdge25");
MyModelTransformerMap.Add(FEdge5, "FEdge26");
MyModelTransformerMap.Add(FEdge6, "FEdge27");
MyModelTransformerMap.Add(FEdge7, "FEdge28");
MyModelTransformerMap.Add(FEdge8, "FEdge29");
MyModelTransformerMap.Add(FEdge9, "FEdge30");
MyModelTransformerMap.Add(FEdge10, "FEdge31");
MyModelTransformerMap.Add(FEdge11, "FEdge32");
MyModelTransformerMap.Add(FEdge12, "FEdge33");
MyModelTransformerMap.Add(FEdge13, "FEdge34");
MyModelTransformerMap.Add(FEdge14, "FEdge35");
MyModelTransformerMap.Add(FEdge15, "FEdge36");
MyModelTransformerMap.Add(FEdge16, "FEdge37");
MyModelTransformerMap.Add(FEdge17, "FEdge38");
MyModelTransformerMap.Add(FEdge18, "FEdge39");
MyModelTransformerMap.Add(FEdge19, "FEdge40");
MyModelTransformerMap.Add(FEdge20, "FEdge41");
MyModelTransformerMap.Add(FEdge21, "FEdge42");
MyModelTransformerMap.Add(Pl1, "Pl70");
MyModelTransformerMap.Add(Pl2, "Pl71");
MyModelTransformerMap.Add(Pl3, "Pl72");
MyModelTransformerMap.Add(Pl4, "Pl73");
MyModelTransformerMap.Add(Pl5, "Pl74");
MyModelTransformerMap.Add(Pl6, "Pl75");
MyModelTransformerMap.Add(Pl7, "Pl76");
MyModelTransformerMap.Add(Pl8, "Pl77");
MyModelTransformerMap.Add(Pl9, "Pl78");
MyModelTransformerMap.Add(Pl10, "Pl79");
MyModelTransformerMap.Add(Pl14, "Pl80");
MyModelTransformerMap.Add(Pl15, "Pl81");
MyModelTransformerMap.Add(Pl16, "Pl82");
MyModelTransformerMap.Add(Pl17, "Pl83");
MyModelTransformerMap.Add(Pl18, "Pl84");
MyModelTransformerMap.Add(Pl19, "Pl85");
MyModelTransformerMap.Add(Pl20, "Pl86");
MyModelTransformerMap.Add(Pl21, "Pl87");
MyModelTransformerMap.Add(Pl22, "Pl88");
MyModelTransformerMap.Add(Pl23, "Pl89");
MyModelTransformerMap.Add(Pl24, "Pl90");
MyModelTransformerMap.Add(Pl25, "Pl91");
MyModelTransformerMap.Add(Pl26, "Pl92");
MyModelTransformerMap.Add(Pl27, "Pl93");
MyModelTransformerMap.Add(Pl28, "Pl94");
MyModelTransformerMap.Add(Pl29, "Pl95");
MyModelTransformerMap.Add(Pl30, "Pl96");
MyModelTransformerMap.Add(Pl31, "Pl97");
MyModelTransformerMap.Add(Pl32, "Pl98");
MyModelTransformerMap.Add(Pl33, "Pl99");
MyModelTransformerMap.Add(Pl34, "Pl100");
MyModelTransformerMap.Add(Pl35, "Pl101");
MyModelTransformerMap.Add(Pl36, "Pl102");
MyModelTransformerMap.Add(Pl37, "Pl103");
MyModelTransformerMap.Add(Pl38, "Pl104");
MyModelTransformerMap.Add(Pl39, "Pl105");
MyModelTransformerMap.Add(Pl40, "Pl106");
MyModelTransformerMap.Add(Pl41, "Pl107");
MyModelTransformerMap.Add(Pl42, "Pl108");
MyModelTransformerMap.Add(Pl43, "Pl109");
MyModelTransformerMap.Add(Pl44, "Pl110");
MyModelTransformerMap.Add(Pl45, "Pl111");
MyModelTransformerMap.Add(Pl46, "Pl112");
MyModelTransformerMap.Add(Pl47, "Pl113");
MyModelTransformerMap.Add(Pl48, "Pl114");
MyModelTransformerMap.Add(Pl49, "Pl115");
MyModelTransformerMap.Add(Pl50, "Pl116");
MyModelTransformerMap.Add(Pl51, "Pl117");
MyModelTransformerMap.Add(Pl52, "Pl118");
MyModelTransformerMap.Add(Pl53, "Pl119");
MyModelTransformerMap.Add(Pl54, "Pl120");
MyModelTransformerMap.Add(Pl55, "Pl121");
MyModelTransformerMap.Add(Pl56, "Pl122");
MyModelTransformerMap.Add(Pl57, "Pl123");
MyModelTransformerMap.Add(Pl58, "Pl124");
MyModelTransformerMap.Add(Pl59, "Pl125");
MyModelTransformerMap.Add(Pl60, "Pl126");
MyModelTransformerMap.Add(Pl61, "Pl127");
MyModelTransformerMap.Add(Pl62, "Pl128");
MyModelTransformerMap.Add(Pl63, "Pl129");
MyModelTransformerMap.Add(Pl64, "Pl130");
MyModelTransformerMap.Add(Pl65, "Pl131");
MyModelTransformerMap.Add(Pl66, "Pl132");
MyModelTransformerMap.Add(Pl134, "Pl138");
MyModelTransformerMap.Add(Pl135, "Pl139");
ModelTransformer(MyModelTransformerMap).copyMirror(Point0, Vector3d(0, 1, 0));
Point43 = Point0.copyMirror(Point0,Vector3d(0, 1, 0));
Point44 = Point1.copyMirror(Point0,Vector3d(0, 1, 0));
Point45 = Point2.copyMirror(Point0,Vector3d(0, 1, 0));
Point46 = Point3.copyMirror(Point0,Vector3d(0, 1, 0));
Point47 = Point4.copyMirror(Point0,Vector3d(0, 1, 0));
Point48 = Point5.copyMirror(Point0,Vector3d(0, 1, 0));
Point49 = Point6.copyMirror(Point0,Vector3d(0, 1, 0));
Point50 = Point7.copyMirror(Point0,Vector3d(0, 1, 0));
Point51 = Point8.copyMirror(Point0,Vector3d(0, 1, 0));
Point52 = Point9.copyMirror(Point0,Vector3d(0, 1, 0));
Point53 = Point10.copyMirror(Point0,Vector3d(0, 1, 0));
Point54 = Point11.copyMirror(Point0,Vector3d(0, 1, 0));
Point55 = Point12.copyMirror(Point0,Vector3d(0, 1, 0));
Point56 = Point13.copyMirror(Point0,Vector3d(0, 1, 0));
Point57 = Point14.copyMirror(Point0,Vector3d(0, 1, 0));
Point58 = Point15.copyMirror(Point0,Vector3d(0, 1, 0));
Point59 = Point16.copyMirror(Point0,Vector3d(0, 1, 0));
Point60 = Point17.copyMirror(Point0,Vector3d(0, 1, 0));
Point61 = Point18.copyMirror(Point0,Vector3d(0, 1, 0));
Point62 = Point19.copyMirror(Point0,Vector3d(0, 1, 0));
Point63 = Point20.copyMirror(Point0,Vector3d(0, 1, 0));
Point64 = Point21.copyMirror(Point0,Vector3d(0, 1, 0));
Point65 = Point22.copyMirror(Point0,Vector3d(0, 1, 0));
Point66 = Point23.copyMirror(Point0,Vector3d(0, 1, 0));
Point67 = Point24.copyMirror(Point0,Vector3d(0, 1, 0));
Point68 = Point25.copyMirror(Point0,Vector3d(0, 1, 0));
Point69 = Point26.copyMirror(Point0,Vector3d(0, 1, 0));
Point70 = Point27.copyMirror(Point0,Vector3d(0, 1, 0));
Point71 = Point28.copyMirror(Point0,Vector3d(0, 1, 0));
Point72 = Point29.copyMirror(Point0,Vector3d(0, 1, 0));
Point73 = Point30.copyMirror(Point0,Vector3d(0, 1, 0));
Point74 = Point31.copyMirror(Point0,Vector3d(0, 1, 0));
Point75 = Point32.copyMirror(Point0,Vector3d(0, 1, 0));
Point76 = Point33.copyMirror(Point0,Vector3d(0, 1, 0));
Point77 = Point34.copyMirror(Point0,Vector3d(0, 1, 0));
Point78 = Point35.copyMirror(Point0,Vector3d(0, 1, 0));
Point79 = Point36.copyMirror(Point0,Vector3d(0, 1, 0));
Point80 = Point37.copyMirror(Point0,Vector3d(0, 1, 0));
Point81 = Point38.copyMirror(Point0,Vector3d(0, 1, 0));
Point82 = Point39.copyMirror(Point0,Vector3d(0, 1, 0));
Point83 = Point40.copyMirror(Point0,Vector3d(0, 1, 0));
Point84 = Point41.copyMirror(Point0,Vector3d(0, 1, 0));
Point85 = Point42.copyMirror(Point0,Vector3d(0, 1, 0));
Curve133 = Curve1.copyMirror(Point0,Vector3d(0, 1, 0));
Curve134 = Curve2.copyMirror(Point0,Vector3d(0, 1, 0));
Curve135 = Curve3.copyMirror(Point0,Vector3d(0, 1, 0));
Curve136 = Curve38.copyMirror(Point0,Vector3d(0, 1, 0));
Curve137 = Curve8.copyMirror(Point0,Vector3d(0, 1, 0));
Curve138 = Curve6.copyMirror(Point0,Vector3d(0, 1, 0));
Curve139 = Curve7.copyMirror(Point0,Vector3d(0, 1, 0));
Curve140 = Curve9.copyMirror(Point0,Vector3d(0, 1, 0));
Curve141 = Curve10.copyMirror(Point0,Vector3d(0, 1, 0));
Curve142 = Curve5.copyMirror(Point0,Vector3d(0, 1, 0));
Curve143 = Curve11.copyMirror(Point0,Vector3d(0, 1, 0));
Curve144 = Curve12.copyMirror(Point0,Vector3d(0, 1, 0));
Curve145 = Curve13.copyMirror(Point0,Vector3d(0, 1, 0));
Curve146 = Curve42.copyMirror(Point0,Vector3d(0, 1, 0));
Curve147 = Curve15.copyMirror(Point0,Vector3d(0, 1, 0));
Curve148 = Curve16.copyMirror(Point0,Vector3d(0, 1, 0));
Curve149 = Curve17.copyMirror(Point0,Vector3d(0, 1, 0));
Curve150 = Curve18.copyMirror(Point0,Vector3d(0, 1, 0));
Curve151 = Curve19.copyMirror(Point0,Vector3d(0, 1, 0));
Curve152 = Curve20.copyMirror(Point0,Vector3d(0, 1, 0));
Curve153 = Curve21.copyMirror(Point0,Vector3d(0, 1, 0));
Curve154 = Curve22.copyMirror(Point0,Vector3d(0, 1, 0));
Curve155 = Curve23.copyMirror(Point0,Vector3d(0, 1, 0));
Curve156 = Curve24.copyMirror(Point0,Vector3d(0, 1, 0));
Curve157 = Curve25.copyMirror(Point0,Vector3d(0, 1, 0));
Curve158 = Curve26.copyMirror(Point0,Vector3d(0, 1, 0));
Curve159 = Curve27.copyMirror(Point0,Vector3d(0, 1, 0));
Curve160 = Curve28.copyMirror(Point0,Vector3d(0, 1, 0));
Curve161 = Curve29.copyMirror(Point0,Vector3d(0, 1, 0));
Curve162 = Curve30.copyMirror(Point0,Vector3d(0, 1, 0));
Curve163 = Curve31.copyMirror(Point0,Vector3d(0, 1, 0));
Curve164 = Curve32.copyMirror(Point0,Vector3d(0, 1, 0));
Curve165 = Curve33.copyMirror(Point0,Vector3d(0, 1, 0));
Curve166 = Curve34.copyMirror(Point0,Vector3d(0, 1, 0));
Curve167 = Curve35.copyMirror(Point0,Vector3d(0, 1, 0));
Curve168 = Curve37.copyMirror(Point0,Vector3d(0, 1, 0));
Curve169 = Curve36.copyMirror(Point0,Vector3d(0, 1, 0));
Curve170 = Curve39.copyMirror(Point0,Vector3d(0, 1, 0));
Curve171 = Curve41.copyMirror(Point0,Vector3d(0, 1, 0));
Curve172 = Curve40.copyMirror(Point0,Vector3d(0, 1, 0));
Curve173 = Curve43.copyMirror(Point0,Vector3d(0, 1, 0));
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
Curve192 = Curve66.copyMirror(Point0,Vector3d(0, 1, 0));
Curve193 = Curve63.copyMirror(Point0,Vector3d(0, 1, 0));
Curve194 = Curve64.copyMirror(Point0,Vector3d(0, 1, 0));
Curve195 = Curve65.copyMirror(Point0,Vector3d(0, 1, 0));
Curve196 = Curve67.copyMirror(Point0,Vector3d(0, 1, 0));
Curve197 = Curve68.copyMirror(Point0,Vector3d(0, 1, 0));
Curve198 = Curve69.copyMirror(Point0,Vector3d(0, 1, 0));
Curve199 = Curve70.copyMirror(Point0,Vector3d(0, 1, 0));
Curve200 = Curve71.copyMirror(Point0,Vector3d(0, 1, 0));
Curve201 = Curve72.copyMirror(Point0,Vector3d(0, 1, 0));
Curve202 = Curve91.copyMirror(Point0,Vector3d(0, 1, 0));
Curve203 = Curve73.copyMirror(Point0,Vector3d(0, 1, 0));
Curve204 = Curve74.copyMirror(Point0,Vector3d(0, 1, 0));
Curve205 = Curve75.copyMirror(Point0,Vector3d(0, 1, 0));
Curve206 = Curve76.copyMirror(Point0,Vector3d(0, 1, 0));
Curve207 = Curve77.copyMirror(Point0,Vector3d(0, 1, 0));
Curve208 = Curve78.copyMirror(Point0,Vector3d(0, 1, 0));
Curve209 = Curve79.copyMirror(Point0,Vector3d(0, 1, 0));
Curve210 = Curve80.copyMirror(Point0,Vector3d(0, 1, 0));
Curve211 = Curve62.copyMirror(Point0,Vector3d(0, 1, 0));
Curve212 = Curve92.copyMirror(Point0,Vector3d(0, 1, 0));
Curve213 = Curve93.copyMirror(Point0,Vector3d(0, 1, 0));
Curve214 = Curve94.copyMirror(Point0,Vector3d(0, 1, 0));
Curve215 = Curve95.copyMirror(Point0,Vector3d(0, 1, 0));
Curve216 = Curve96.copyMirror(Point0,Vector3d(0, 1, 0));
Curve217 = Curve97.copyMirror(Point0,Vector3d(0, 1, 0));
Curve218 = Curve98.copyMirror(Point0,Vector3d(0, 1, 0));
Curve219 = Curve99.copyMirror(Point0,Vector3d(0, 1, 0));
Curve220 = Curve100.copyMirror(Point0,Vector3d(0, 1, 0));
Curve221 = Curve101.copyMirror(Point0,Vector3d(0, 1, 0));
Curve222 = Curve102.copyMirror(Point0,Vector3d(0, 1, 0));
Curve223 = Curve103.copyMirror(Point0,Vector3d(0, 1, 0));
Curve224 = Curve104.copyMirror(Point0,Vector3d(0, 1, 0));
Curve225 = Curve105.copyMirror(Point0,Vector3d(0, 1, 0));
Curve226 = Curve106.copyMirror(Point0,Vector3d(0, 1, 0));
Curve227 = Curve107.copyMirror(Point0,Vector3d(0, 1, 0));
Curve228 = Curve108.copyMirror(Point0,Vector3d(0, 1, 0));
Curve229 = Curve109.copyMirror(Point0,Vector3d(0, 1, 0));
Curve230 = Curve110.copyMirror(Point0,Vector3d(0, 1, 0));
Curve231 = Curve111.copyMirror(Point0,Vector3d(0, 1, 0));
Curve232 = Curve112.copyMirror(Point0,Vector3d(0, 1, 0));
Curve233 = Curve113.copyMirror(Point0,Vector3d(0, 1, 0));
Curve234 = Curve114.copyMirror(Point0,Vector3d(0, 1, 0));
Curve235 = Curve115.copyMirror(Point0,Vector3d(0, 1, 0));
Curve236 = Curve116.copyMirror(Point0,Vector3d(0, 1, 0));
Curve237 = Curve117.copyMirror(Point0,Vector3d(0, 1, 0));
Curve238 = Curve118.copyMirror(Point0,Vector3d(0, 1, 0));
Curve239 = Curve119.copyMirror(Point0,Vector3d(0, 1, 0));
Curve240 = Curve120.copyMirror(Point0,Vector3d(0, 1, 0));
Curve241 = Curve121.copyMirror(Point0,Vector3d(0, 1, 0));
Curve242 = Curve122.copyMirror(Point0,Vector3d(0, 1, 0));
Curve245 = Curve123.copyMirror(Point0,Vector3d(0, 1, 0));
Curve246 = Curve124.copyMirror(Point0,Vector3d(0, 1, 0));
Curve247 = Curve125.copyMirror(Point0,Vector3d(0, 1, 0));
Curve248 = Curve126.copyMirror(Point0,Vector3d(0, 1, 0));
Curve249 = Curve127.copyMirror(Point0,Vector3d(0, 1, 0));
Curve250 = Curve243.copyMirror(Point0,Vector3d(0, 1, 0));
Curve251 = Curve244.copyMirror(Point0,Vector3d(0, 1, 0));
Delete(Pl75);
Delete(Pl70);
Delete(Pl100);
Delete(Pl79);
Delete(Pl73);
Delete(Pl105);
Pl116 = SweepCurve(Curve133, Curve203);
Pl118 = SweepCurve(Curve137, Curve203);
Pl122 = SweepCurve(Curve163, Curve203);
Pl118.flipNormal();
Pl139 = SweepCurve(Curve141, Curve160);
Pl140 = SweepCurve(Curve169, Curve160);
Pl75 = SweepCurve(Curve168, Curve209);
Pl75.flipNormal();
Pl2.join(Pl71);
Pl2.simplifyTopology();
Pl31.join(Pl97);
Pl31.simplifyTopology();


MyModelTransformerMap = ObjectNameMap();
MyModelTransformerMap.Add(Pl67, "Pl70");
MyModelTransformerMap.Add(Pl68, "Pl71");
MyModelTransformerMap.Add(Pl69, "Pl73");
ModelTransformer(MyModelTransformerMap).copyMirror(Point0, Vector3d(0, 1, 0));
Curve257 = Curve254.copyMirror(Point0,Vector3d(0, 1, 0));
Curve258 = Curve252.copyMirror(Point0,Vector3d(0, 1, 0));
Curve259 = Curve253.copyMirror(Point0,Vector3d(0, 1, 0));
Curve260 = Curve255.copyMirror(Point0,Vector3d(0, 1, 0));
Curve261 = Curve256.copyMirror(Point0,Vector3d(0, 1, 0));


wetted_surface.add(Pl118);
wetted_surface.add(Pl122);
wetted_surface.add(Pl139);
wetted_surface.add(Pl140);
wetted_surface.add(Pl72);
wetted_surface.add(Pl74);
wetted_surface.add(Pl76);
wetted_surface.add(Pl77);
wetted_surface.add(Pl78);
wetted_surface.add(Pl80);
wetted_surface.add(Pl81);
wetted_surface.add(Pl82);
wetted_surface.add(Pl83);
wetted_surface.add(Pl84);
wetted_surface.add(Pl85);
wetted_surface.add(Pl86);
wetted_surface.add(Pl87);
wetted_surface.add(Pl88);
wetted_surface.add(Pl89);
wetted_surface.add(Pl90);
wetted_surface.add(Pl91);
wetted_surface.add(Pl92);
wetted_surface.add(Pl93);
wetted_surface.add(Pl94);
wetted_surface.add(Pl95);
wetted_surface.add(Pl96);

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

Delete(Pl_fs3);
Delete(Pl_fs5);
Delete(Pl_fs4);
Delete(Pl_fs2);
Delete(Pl_fs6);
Delete(Pl_fs12);
Delete(Pl_fs11);
Delete(Pl_fs8);
Delete(Pl_fs7);
Delete(Pl_fs9);
Delete(Pl_fs18);
Delete(Pl_fs17);
Delete(Pl_fs13);
Delete(Pl_fs15);
Delete(Pl_fs14);
Delete(Pl_fs16);
Delete(Pl_fs19);
Delete(Pl_fs21);
Delete(Pl_fs20);

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

Pl_fs10.front.wetSurface = outer_wet;
Free_surface = Set();
Free_surface.add(Pl_fs10);
Meshing_panels.step(1).subset = Free_surface;

LC1.generateAppliedLoads();
LC1.setFemLoadcase(1);

Meshing_panels.step(1).subset = Free_surface;

Pl_fs10.meshDensity = FS_density;
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

Pl_int_lid1 = CoverCurves(Curve45,Curve46,Curve47,Curve50);
Pl_int_lid2 = CoverCurves(Curve48,Curve49,Curve54);
Pl_int_lid3 = CoverCurves(Curve51,Curve52,Curve53);
Pl_int_lid2.flipNormal();
Pl_int_lid1.flipNormal();

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
LC3.compartment(cm.compartment(Point(-0.9351947381 m,-4.177264884 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC4 = LoadCase();
LC4.compartment(cm.compartment(Point(-0.5328037944 m,0.4861657039 m,10 m))).globalIntensity = DummyHydroPressure();
LC5 = LoadCase();
LC5.compartment(cm.compartment(Point(-1.405334509 m,3.711816231 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC6 = LoadCase();
LC6.compartment(cm.compartment(Point(-3.23127603 m,-2.71058078 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC7 = LoadCase();
LC7.compartment(cm.compartment(Point(-3.771766569 m,1.347247116 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC8 = LoadCase();
LC8.compartment(cm.compartment(Point(-23.90334428 m,-0.8241703011 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC9 = LoadCase();
LC9.compartment(cm.compartment(Point(-23.90334428 m,5.425829699 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC10 = LoadCase();
LC10.compartment(cm.compartment(Point(-33.22228559 m,-5.203457522 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC11 = LoadCase();
LC11.compartment(cm.compartment(Point(-33.22228559 m,1.046542478 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC12 = LoadCase();
LC12.compartment(cm.compartment(Point(-46.34967024 m,0.7570706222 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC13 = LoadCase();
LC13.compartment(cm.compartment(Point(-46.97981953 m,-1.837343063 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC14 = LoadCase();
LC14.compartment(cm.compartment(Point(-49.47732565 m,-5.150199796 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC15 = LoadCase();
LC15.compartment(cm.compartment(Point(-50.60146114 m,-1.417510011 m,-8 m))).globalIntensity = DummyHydroPressure();
LC16 = LoadCase();
LC16.compartment(cm.compartment(Point(-50.99292938 m,5.400329762 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC17 = LoadCase();
LC17.compartment(cm.compartment(Point(-52.2561636 m,0.4618574187 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC18 = LoadCase();
LC18.compartment(cm.compartment(Point(-52.41600474 m,0.6077071298 m,10 m))).globalIntensity = DummyHydroPressure();
LC19 = LoadCase();
LC19.compartment(cm.compartment(Point(-52.73321549 m,-5.512343624 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC20 = LoadCase();
LC20.compartment(cm.compartment(Point(-53.58734306 m,4.770180465 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC21 = LoadCase();
LC21.compartment(cm.compartment(Point(-55.92110884 m,-3.521549667 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC22 = LoadCase();
LC22.compartment(cm.compartment(Point(-56.63464926 m,0.144687735 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC23 = LoadCase();
LC23.compartment(cm.compartment(Point(0.3515554647 m,-0.3895405897 m,-8 m))).globalIntensity = DummyHydroPressure();
LC24 = LoadCase();
LC24.compartment(cm.compartment(Point(0.1018733219 m,4.787375642 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC25 = LoadCase();
LC25.compartment(cm.compartment(Point(1.347247116 m,-3.771766569 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC26 = LoadCase();
LC26.compartment(cm.compartment(Point(3.052761523 m,14.23639477 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC27 = LoadCase();
LC27.compartment(cm.compartment(Point(3.711816231 m,-1.405334509 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC28 = LoadCase();
LC28.compartment(cm.compartment(Point(4.783166336 m,0.09967199432 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC29 = LoadCase();
LC29.compartment(cm.compartment(Point(7.093783329 m,2.435111987 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC30 = LoadCase();
LC30.compartment(cm.compartment(Point(7.728704328 m,-19.63650857 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC31 = LoadCase();
LC31.compartment(cm.compartment(Point(13.0344682 m,-16.32636117 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC32 = LoadCase();
LC32.compartment(cm.compartment(Point(13.2367705 m,23.82126024 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC33 = LoadCase();
LC33.compartment(cm.compartment(Point(17.81195784 m,-37.10121595 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC34 = LoadCase();
LC34.compartment(cm.compartment(Point(18.46943129 m,20.81637196 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC35 = LoadCase();
LC35.compartment(cm.compartment(Point(20.89912542 m,-47.05322727 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC36 = LoadCase();
LC36.compartment(cm.compartment(Point(21.2055304 m,-43.16575482 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC37 = LoadCase();
LC37.compartment(cm.compartment(Point(21.07382814 m,43.68149346 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC38 = LoadCase();
LC38.compartment(cm.compartment(Point(21.80210105 m,41.16727536 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC39 = LoadCase();
LC39.compartment(cm.compartment(Point(22.30474762 m,-48.95055258 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC40 = LoadCase();
LC40.compartment(cm.compartment(Point(22.66802706 m,49.08622516 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC41 = LoadCase();
LC41.compartment(cm.compartment(Point(22.95162617 m,-39.39045045 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC42 = LoadCase();
LC42.compartment(cm.compartment(Point(23.25266509 m,-34.02479733 m,-16.5 m))).globalIntensity = DummyHydroPressure();
LC43 = LoadCase();
LC43.compartment(cm.compartment(Point(23.63901867 m,39.84076142 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC44 = LoadCase();
LC44.compartment(cm.compartment(Point(24.48204899 m,50.73774922 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC45 = LoadCase();
LC45.compartment(cm.compartment(Point(24.74158407 m,-49.62152588 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC46 = LoadCase();
LC46.compartment(cm.compartment(Point(25.20899523 m,45.4245218 m,10 m))).globalIntensity = DummyHydroPressure();
LC47 = LoadCase();
LC47.compartment(cm.compartment(Point(25.20899528 m,-44.20910754 m,10 m))).globalIntensity = DummyHydroPressure();
LC48 = LoadCase();
LC48.compartment(cm.compartment(Point(25.36883642 m,-44.35495725 m,1 m))).globalIntensity = DummyHydroPressure();
LC49 = LoadCase();
LC49.compartment(cm.compartment(Point(25.41361435 m,-44.37807264 m,-8 m))).globalIntensity = DummyHydroPressure();
LC50 = LoadCase();
LC50.compartment(cm.compartment(Point(26.3368574 m,44.31065107 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC51 = LoadCase();
LC51.compartment(cm.compartment(Point(26.07274343 m,44.59074408 m,-8 m))).globalIntensity = DummyHydroPressure();
LC52 = LoadCase();
LC52.compartment(cm.compartment(Point(28.18112177 m,-39.86981454 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC53 = LoadCase();
LC53.compartment(cm.compartment(Point(28.44489112 m,-48.88473468 m,2.5 m))).globalIntensity = DummyHydroPressure();
LC54 = LoadCase();
LC54.compartment(cm.compartment(Point(29.20013198 m,40.42402426 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC55 = LoadCase();
LC55.compartment(cm.compartment(Point(31.3013642 m,41.89344084 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC56 = LoadCase();
LC56.compartment(cm.compartment(Point(31.04358204 m,46.58799032 m,-1.5 m))).globalIntensity = DummyHydroPressure();
LC57 = LoadCase();
LC57.compartment(cm.compartment(Point(31.79593457 m,-46.20976566 m,2.5 m))).globalIntensity = DummyHydroPressure();
compartments = Set();
compartments.add(LC57.compartment(Point(-0.3729626562 m,0.3403159927 m,2.5 m)));
compartments.add(LC57.compartment(Point(-0.5328037944 m,0.4861657039 m,10 m)));
compartments.add(LC57.compartment(Point(-0.9351947381 m,-4.177264884 m,2.5 m)));
compartments.add(LC57.compartment(Point(-1.405334509 m,3.711816231 m,2.5 m)));
compartments.add(LC57.compartment(Point(-23.90334428 m,-0.8241703011 m,-16.5 m)));
compartments.add(LC57.compartment(Point(-23.90334428 m,5.425829699 m,-16.5 m)));
compartments.add(LC57.compartment(Point(-3.23127603 m,-2.71058078 m,2.5 m)));
compartments.add(LC57.compartment(Point(-3.771766569 m,1.347247116 m,-1.5 m)));
compartments.add(LC57.compartment(Point(-33.22228559 m,-5.203457522 m,-16.5 m)));
compartments.add(LC57.compartment(Point(-33.22228559 m,1.046542478 m,-16.5 m)));
compartments.add(LC57.compartment(Point(-46.34967024 m,0.7570706222 m,2.5 m)));
compartments.add(LC57.compartment(Point(-46.97981953 m,-1.837343063 m,2.5 m)));
compartments.add(LC57.compartment(Point(-49.47732565 m,-5.150199796 m,2.5 m)));
compartments.add(LC57.compartment(Point(-50.60146114 m,-1.417510011 m,-8 m)));
compartments.add(LC57.compartment(Point(-50.99292938 m,5.400329762 m,2.5 m)));
compartments.add(LC57.compartment(Point(-52.2561636 m,0.4618574187 m,2.5 m)));
compartments.add(LC57.compartment(Point(-52.41600474 m,0.6077071298 m,10 m)));
compartments.add(LC57.compartment(Point(-52.73321549 m,-5.512343624 m,2.5 m)));
compartments.add(LC57.compartment(Point(-53.58734306 m,4.770180465 m,2.5 m)));
compartments.add(LC57.compartment(Point(-55.92110884 m,-3.521549667 m,2.5 m)));
compartments.add(LC57.compartment(Point(-56.63464926 m,0.144687735 m,2.5 m)));
compartments.add(LC57.compartment(Point(0.3515554647 m,-0.3895405897 m,-8 m)));
compartments.add(LC57.compartment(Point(0.3698946273 m,4.374586472 m,2.5 m)));
compartments.add(LC57.compartment(Point(1.347247116 m,-3.771766569 m,2.5 m)));
compartments.add(LC57.compartment(Point(13.0344682 m,-16.32636117 m,-16.5 m)));
compartments.add(LC57.compartment(Point(13.2367705 m,23.82126024 m,-16.5 m)));
compartments.add(LC57.compartment(Point(17.81195784 m,-37.10121595 m,-16.5 m)));
compartments.add(LC57.compartment(Point(18.46943129 m,20.81637196 m,-16.5 m)));
compartments.add(LC57.compartment(Point(20.89912542 m,-47.05322727 m,2.5 m)));
compartments.add(LC57.compartment(Point(21.07382813 m,43.68149345 m,2.5 m)));
compartments.add(LC57.compartment(Point(21.2055304 m,-43.16575482 m,2.5 m)));
compartments.add(LC57.compartment(Point(21.80210105 m,41.16727536 m,2.5 m)));
compartments.add(LC57.compartment(Point(22.30474762 m,-48.95055258 m,2.5 m)));
compartments.add(LC57.compartment(Point(22.66802705 m,49.08622514 m,2.5 m)));
compartments.add(LC57.compartment(Point(22.95162617 m,-39.39045045 m,2.5 m)));
compartments.add(LC57.compartment(Point(23.25266509 m,-34.02479733 m,-16.5 m)));
compartments.add(LC57.compartment(Point(23.63901867 m,39.84076142 m,-1.5 m)));
compartments.add(LC57.compartment(Point(24.48204899 m,50.73774922 m,-1.5 m)));
compartments.add(LC57.compartment(Point(24.74158407 m,-49.62152588 m,2.5 m)));
compartments.add(LC57.compartment(Point(25.20899523 m,45.4245218 m,10 m)));
compartments.add(LC57.compartment(Point(25.20899528 m,-44.20910754 m,10 m)));
compartments.add(LC57.compartment(Point(25.36883642 m,-44.35495725 m,1 m)));
compartments.add(LC57.compartment(Point(25.41361435 m,-44.37807264 m,-8 m)));
compartments.add(LC57.compartment(Point(26.07274343 m,44.59074408 m,-8 m)));
compartments.add(LC57.compartment(Point(26.3368574 m,44.31065107 m,2.5 m)));
compartments.add(LC57.compartment(Point(28.18112177 m,-39.86981454 m,2.5 m)));
compartments.add(LC57.compartment(Point(28.44489112 m,-48.88473468 m,2.5 m)));
compartments.add(LC57.compartment(Point(29.20013199 m,40.42402426 m,2.5 m)));
compartments.add(LC57.compartment(Point(3.052761523 m,14.23639477 m,-16.5 m)));
compartments.add(LC57.compartment(Point(3.711816231 m,-1.405334509 m,2.5 m)));
compartments.add(LC57.compartment(Point(31.04358204 m,46.58799032 m,2.5 m)));
compartments.add(LC57.compartment(Point(31.3013642 m,41.89344084 m,-1.5 m)));
compartments.add(LC57.compartment(Point(31.79593457 m,-46.20976566 m,2.5 m)));
compartments.add(LC57.compartment(Point(4.783166336 m,0.09967199432 m,-1.5 m)));
compartments.add(LC57.compartment(Point(7.093783329 m,2.435111987 m,-16.5 m)));
compartments.add(LC57.compartment(Point(7.728704328 m,-19.63650857 m,-16.5 m)));
compartments.add(cm.compartment(Point(-0.3729626562 m,0.3403159927 m,2.5 m)));
compartments.add(cm.compartment(Point(-0.5328037944 m,0.4861657039 m,10 m)));
compartments.add(cm.compartment(Point(-0.9351947381 m,-4.177264884 m,2.5 m)));
compartments.add(cm.compartment(Point(-1.405334509 m,3.711816231 m,2.5 m)));
compartments.add(cm.compartment(Point(-23.90334428 m,-0.8241703011 m,-16.5 m)));
compartments.add(cm.compartment(Point(-23.90334428 m,5.425829699 m,-16.5 m)));
compartments.add(cm.compartment(Point(-3.23127603 m,-2.71058078 m,2.5 m)));
compartments.add(cm.compartment(Point(-3.771766569 m,1.347247116 m,-1.5 m)));
compartments.add(cm.compartment(Point(-33.22228559 m,-5.203457522 m,-16.5 m)));
compartments.add(cm.compartment(Point(-33.22228559 m,1.046542478 m,-16.5 m)));
compartments.add(cm.compartment(Point(-46.34967024 m,0.7570706222 m,2.5 m)));
compartments.add(cm.compartment(Point(-46.97981953 m,-1.837343063 m,2.5 m)));
compartments.add(cm.compartment(Point(-49.47732565 m,-5.150199796 m,2.5 m)));
compartments.add(cm.compartment(Point(-50.60146114 m,-1.417510011 m,-8 m)));
compartments.add(cm.compartment(Point(-50.99292938 m,5.400329762 m,2.5 m)));
compartments.add(cm.compartment(Point(-52.2561636 m,0.4618574187 m,2.5 m)));
compartments.add(cm.compartment(Point(-52.41600474 m,0.6077071298 m,10 m)));
compartments.add(cm.compartment(Point(-52.73321549 m,-5.512343624 m,2.5 m)));
compartments.add(cm.compartment(Point(-53.58734306 m,4.770180465 m,2.5 m)));
compartments.add(cm.compartment(Point(-55.92110884 m,-3.521549667 m,2.5 m)));
compartments.add(cm.compartment(Point(-56.63464926 m,0.144687735 m,2.5 m)));
compartments.add(cm.compartment(Point(0.3515554647 m,-0.3895405897 m,-8 m)));
compartments.add(cm.compartment(Point(0.3698946273 m,4.374586472 m,2.5 m)));
compartments.add(cm.compartment(Point(1.347247116 m,-3.771766569 m,2.5 m)));
compartments.add(cm.compartment(Point(13.0344682 m,-16.32636117 m,-16.5 m)));
compartments.add(cm.compartment(Point(13.2367705 m,23.82126024 m,-16.5 m)));
compartments.add(cm.compartment(Point(17.81195784 m,-37.10121595 m,-16.5 m)));
compartments.add(cm.compartment(Point(18.46943129 m,20.81637196 m,-16.5 m)));
compartments.add(cm.compartment(Point(20.89912542 m,-47.05322727 m,2.5 m)));
compartments.add(cm.compartment(Point(21.07382813 m,43.68149345 m,2.5 m)));
compartments.add(cm.compartment(Point(21.2055304 m,-43.16575482 m,2.5 m)));
compartments.add(cm.compartment(Point(21.80210105 m,41.16727536 m,2.5 m)));
compartments.add(cm.compartment(Point(22.30474762 m,-48.95055258 m,2.5 m)));
compartments.add(cm.compartment(Point(22.66802705 m,49.08622514 m,2.5 m)));
compartments.add(cm.compartment(Point(22.95162617 m,-39.39045045 m,2.5 m)));
compartments.add(cm.compartment(Point(23.25266509 m,-34.02479733 m,-16.5 m)));
compartments.add(cm.compartment(Point(23.63901867 m,39.84076142 m,-1.5 m)));
compartments.add(cm.compartment(Point(24.48204899 m,50.73774922 m,-1.5 m)));
compartments.add(cm.compartment(Point(24.74158407 m,-49.62152588 m,2.5 m)));
compartments.add(cm.compartment(Point(25.20899523 m,45.4245218 m,10 m)));
compartments.add(cm.compartment(Point(25.20899528 m,-44.20910754 m,10 m)));
compartments.add(cm.compartment(Point(25.36883642 m,-44.35495725 m,1 m)));
compartments.add(cm.compartment(Point(25.41361435 m,-44.37807264 m,-8 m)));
compartments.add(cm.compartment(Point(26.07274343 m,44.59074408 m,-8 m)));
compartments.add(cm.compartment(Point(26.3368574 m,44.31065107 m,2.5 m)));
compartments.add(cm.compartment(Point(28.18112177 m,-39.86981454 m,2.5 m)));
compartments.add(cm.compartment(Point(28.44489112 m,-48.88473468 m,2.5 m)));
compartments.add(cm.compartment(Point(29.20013199 m,40.42402426 m,2.5 m)));
compartments.add(cm.compartment(Point(3.052761523 m,14.23639477 m,-16.5 m)));
compartments.add(cm.compartment(Point(3.711816231 m,-1.405334509 m,2.5 m)));
compartments.add(cm.compartment(Point(31.04358204 m,46.58799032 m,2.5 m)));
compartments.add(cm.compartment(Point(31.3013642 m,41.89344084 m,-1.5 m)));
compartments.add(cm.compartment(Point(31.79593457 m,-46.20976566 m,2.5 m)));
compartments.add(cm.compartment(Point(4.783166336 m,0.09967199432 m,-1.5 m)));
compartments.add(cm.compartment(Point(7.093783329 m,2.435111987 m,-16.5 m)));
compartments.add(cm.compartment(Point(7.728704328 m,-19.63650857 m,-16.5 m)));

Meshing_panels.step(1).subset = Free_surface;
Meshing_panels.step(1).excludeIncludeMeshSubsetOption = anExcludeMeshSubset;
SimplifyTopology();
SplitPeriodicGeometry();

GenieRules.Meshing.superElementType = 5;

Meshing_panels.execute();
ExportMeshFem().DoExport("T5.FEM");