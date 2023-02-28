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

d1 = 10.0 m; // Diameter of center column
d2 = 12.5 m; // Diamter of radial columns
d_brace = 0.91 m; // Diameter of braces

w_pontoon = 12.5 m; // Horizontal width of pontoon
h_pontoon = 7.00 m;  // Vertical height of pontoon

draft = 20.00 m; // Draft of floater in operating condition
freeboard = 15.00 m; // Freeboard of floater in operating condition

cc_distance = 51.75 m;  // Center-center distance radial to center column

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