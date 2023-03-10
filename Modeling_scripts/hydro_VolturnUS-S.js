// Include following part if scripts is used

Workspace = Reference("/Workspace");
Environment = new EnvironmentFolder(Workspace, "Environment");
Environment.Description = null;
//
Workspace = Reference("/Workspace");
Properties = new CompartmentPropertiesFolder(Workspace, "Properties");
Properties.Description = null;
//
Workspace = Reference("/Workspace");
StabilityProperties = new StabilityPropertiesFolder(Workspace, "StabilityProperties");
StabilityProperties.Description = null;
//
Workspace = Reference("/Workspace");
Models = new ModelsFolder(Workspace, "Models");
Models.Description = null;
//
Workspace = Reference("/Workspace");
HydroModels = new HydroModelsFolder(Workspace, "HydroModels");
HydroModels.Description = null;
//
Workspace = Reference("/Workspace");
Rules = new RulesFolder(Workspace, "Rules");
Rules.Description = null;
//
Workspace = Reference("/Workspace");
Analyses = new AnalysesFolder(Workspace, "Analyses");
Analyses.Description = null;
//
Workspace = Reference("/Workspace");
Plots = new PlotsFolder(Workspace, "Plots");
Plots.Description = null;
//
Workspace = Reference("/Workspace");
Tables = new TablesFolder(Workspace, "Tables");
Tables.Description = null;
//
Workspace = Reference("/Workspace");
Pictures = new PicturesFolder(Workspace, "Pictures");
Pictures.Description = null;
//
Workspace = Reference("/Workspace");
Reports = new ReportsFolder(Workspace, "Reports");
Reports.Description = null;
//
Workspace = Reference("/Workspace");
ViewSettings = new ViewSettingsFolder(Workspace, "ViewSettings");
ViewSettings.Description = null;
//
ViewSettings = Reference("/Workspace/ViewSettings");
ViewSetting1 = new ViewSetting(ViewSettings, "ViewSetting1");
ViewSetting1.BackgroundBottomColor = Color.FromArgb(255, 0, 0, 0);
ViewSetting1.BackgroundColor = Color.FromArgb(255, 0, 0, 0);
ViewSetting1.BackgroundGradient = false;
ViewSetting1.BackgroundTopColor = Color.FromArgb(255, 169, 169, 169);
ViewSetting1.BillboardLabelDepthOffset = 0 m;
ViewSetting1.BillboardLabelHeight = 10;
ViewSetting1.CameraPosition = new Position3D(240 m, -120 m, 80 m);
ViewSetting1.CameraType = CameraType.Perspective;
ViewSetting1.CoordinateSystemLabelBackground = Color.FromArgb(0, 0, 0, 1);
ViewSetting1.CoordinateSystemLabelColor = Color.FromArgb(0, 0, 0, 1);
ViewSetting1.Description = null;
ViewSetting1.FarPlaneDistance = 1000000 m;
ViewSetting1.FieldOfView = 45 deg;
ViewSetting1.GridColor = Color.FromArgb(0, 0, 0, 1);
ViewSetting1.GridLength = 100 m;
ViewSetting1.GridMajorInterval = 10 m;
ViewSetting1.GridMinorInterval = 1 m;
ViewSetting1.GridThickness = 0.04 m;
ViewSetting1.GridWidth = 100 m;
ViewSetting1.InputCoordinateSystemDiameter = 0.5 m;
ViewSetting1.InputCoordinateSystemLabel = "Input system";
ViewSetting1.InputCoordinateSystemSize = 3 m;
ViewSetting1.LabelBackground = null;
ViewSetting1.LabelBorderColor = null;
ViewSetting1.LabelColor = Color.FromArgb(255, 0, 0, 0);
ViewSetting1.LabelPinColor = null;
ViewSetting1.LineThicknessScale = 0.1 m;
ViewSetting1.LookDirection = new Direction3D(-240 m, 120 m, -80 m);
ViewSetting1.MaxNumberOfLabels = 10000;
ViewSetting1.NearPlaneDistance = 0.1 m;
ViewSetting1.OpenWhenLoaded = true;
ViewSetting1.OrthographicWidth = 100 m;
ViewSetting1.ShowCoordinateSystem = true;
ViewSetting1.ShowGrid = false;
ViewSetting1.ShowInputCoordinateSystem = false;
ViewSetting1.SpatialLabelHeight = 1 m;
ViewSetting1.UpDirection = new Vector3D(0, 0, 1);
ViewSetting1.ViewSubtitle = null;
ViewSetting1.ViewTitle = null;
ViewSetting1.ViewTitleColor = Color.FromArgb(0, 0, 0, 1);
//
ViewSetting1 = Reference("/Workspace/ViewSettings/ViewSetting1");
HydroViewSetting1 = new HydroViewSetting(ViewSetting1, "HydroViewSetting1");
HydroViewSetting1.BaseLineColor = Color.FromArgb(0, 0, 0, 1);
HydroViewSetting1.BeamColor = Color.FromArgb(255, 128, 128, 128);
HydroViewSetting1.BeamDivisions = 18;
HydroViewSetting1.BeamThickness = 1;
HydroViewSetting1.BuoyancyModelColor = Color.FromArgb(255, 247, 150, 70);
HydroViewSetting1.CompartmentDamageColor = Color.FromArgb(255, 255, 99, 71);
HydroViewSetting1.CompartmentDeckTankColor = Color.FromArgb(255, 255, 99, 71);
HydroViewSetting1.CompartmentEdgeColor = Color.FromArgb(0, 0, 0, 1);
HydroViewSetting1.CompartmentEdgeThickness = 2;
HydroViewSetting1.CompartmentFloodedByOpeningColor = Color.FromArgb(255, 255, 99, 71);
HydroViewSetting1.CompartmentPlateBackColor = Color.FromArgb(255, 0, 255, 127);
HydroViewSetting1.CompartmentPlateColor = Color.FromArgb(255, 216, 191, 216);
HydroViewSetting1.CompartmentSharpEdgeLimit = 45 deg;
HydroViewSetting1.DeckEdgePointSize = 1 m;
HydroViewSetting1.Description = null;
HydroViewSetting1.GlobalCoordinateSystemDiameter = 0.5 m;
HydroViewSetting1.GlobalCoordinateSystemSize = 3 m;
HydroViewSetting1.HydroPressureLoadArrowColor = Color.FromArgb(255, 173, 216, 230);
HydroViewSetting1.HydroPressureLoadArrowDiameter = 0.5 m;
HydroViewSetting1.HydroPressureLoadArrowDivisions = 4;
HydroViewSetting1.HydroPressureLoadArrowLength = 3 m;
HydroViewSetting1.KeelLineColor = Color.FromArgb(0, 0, 0, 1);
HydroViewSetting1.KeelLineDiameter = 0.1 m;
HydroViewSetting1.LoadCrossSectionArrowDiameter = 2 m;
HydroViewSetting1.LoadCrossSectionArrowLength = 10 m;
HydroViewSetting1.LoadCrossSectionFrameColor = Color.FromArgb(255, 0, 0, 0);
HydroViewSetting1.LoadCrossSectionFrameThickness = 1;
HydroViewSetting1.LoadCrossSectionPointSize = 3 m;
HydroViewSetting1.MarkerColor = Color.FromArgb(255, 255, 99, 71);
HydroViewSetting1.MarkerDiameter = 3 m;
HydroViewSetting1.OpeningDiameter = 3 m;
HydroViewSetting1.PanelEdgeColor = Color.FromArgb(255, 0, 0, 0);
HydroViewSetting1.PanelEdgeThickness = 1;
HydroViewSetting1.PanelShrinkFactor = new Fraction(0.1);
HydroViewSetting1.PerpendicularHeight = 30 m;
HydroViewSetting1.PlateBackColor = Color.FromArgb(255, 255, 165, 0);
HydroViewSetting1.PlateColor = Color.FromArgb(255, 178, 34, 34);
HydroViewSetting1.PointMassColor = Color.FromArgb(255, 147, 137, 83);
HydroViewSetting1.PointMassDiameter = 5 m;
HydroViewSetting1.PointMassDivisions = 16;
HydroViewSetting1.SolidBackColor = Color.FromArgb(255, 0, 0, 255);
HydroViewSetting1.SolidColor = Color.FromArgb(255, 0, 128, 0);
HydroViewSetting1.SpecifyCompartmentPlateBackColor = false;
HydroViewSetting1.SpecifyPlateBackColor = true;
HydroViewSetting1.SpecifySolidBackColor = true;
HydroViewSetting1.ThrusterForceColor = Color.FromArgb(255, 255, 165, 0);
HydroViewSetting1.ThrusterForceRadius = 3 m;
HydroViewSetting1.UnprotectedOpeningColor = Color.FromArgb(255, 255, 165, 0);
HydroViewSetting1.WatertightOpeningColor = Color.FromArgb(255, 0, 128, 0);
HydroViewSetting1.WeathertightOpeningColor = Color.FromArgb(255, 255, 0, 0);
//
Workspace = Reference("/Workspace");
Wizards = new WizardsFolder(Workspace, "Wizards");
Wizards.Description = null;

//


FEM_path = ""

Mass_type = "Including_full_ballast";  // "Including_full_ballast" / "Including_fixed_ballast" / "Excluding_ballast"

Drag_type = "Morison";  // "Critical" / "Morison"

mean_drift = true;  // true / false

diff_freq = true;  // true / false

execute = false;  // true / false

tank_method = "Quasi-static";  // "Dynamic" / "Quasi-static"

anchor_formulation = "Point-mass"; // "Point-mass" / "Anchor-element"

Hs = 3 m;
Tp = 10 s;

Directions1 = new Directions(Environment, "Directions1");
DirectionSet1 = new DirectionSet(Directions1, "DirectionSet1");
DirectionSet1.Items.Add(new AngleInterval());
DirectionSet1.Items[0].From = 0 deg;
DirectionSet1.Items[0].To = 180 deg;
DirectionSet1.Items[0].Step = 45 deg;

Water1 = new Water(Environment, "Water1");

FrequencySet1 = new FrequencySet(Water1, "FrequencySet1");
FrequencySet1.WavePeriodItems.Add(new TimeInterval());
FrequencySet1.WavePeriodItems[0].From = 4 s;
FrequencySet1.WavePeriodItems[0].To = 10 s;
FrequencySet1.WavePeriodItems[0].Step = 2 s;

FrequencySet1.WavePeriodItems.Add(new TimeInterval());
FrequencySet1.WavePeriodItems[1].From = 10 s;
FrequencySet1.WavePeriodItems[1].To = 30 s;
FrequencySet1.WavePeriodItems[1].Step = 1 s;

FrequencySet1.WavePeriodItems.Add(new TimeInterval());
FrequencySet1.WavePeriodItems[2].From = 30 s;
FrequencySet1.WavePeriodItems[2].To = 50 s;
FrequencySet1.WavePeriodItems[2].Step = 2 s;

FrequencySet1.WavePeriodItems.Add(new TimeInterval());
FrequencySet1.WavePeriodItems[3].From = 50 s;
FrequencySet1.WavePeriodItems[3].To = 120 s;
FrequencySet1.WavePeriodItems[3].Step = 5 s;

Location1 = new Location(Environment, "Location1");
Location1.Gravity = 9.80665 m/s^2;
Location1.WaterDepth = 300 m;

FrequencyDomainCondition1 = new FrequencyDomainCondition(Location1, "FrequencyDomainCondition1");
FrequencyDomainCondition1.DirectionSet = DirectionSet1;
FrequencyDomainCondition1.FrequencySet = FrequencySet1;

BretschneiderSpectrum1 = new BretschneiderSpectrum(Water1, "BretschneiderSpectrum1");
BretschneiderSpectrum1.Hs = Hs;
BretschneiderSpectrum1.Tp = Tp;

SeaState1 = new SeaState(Location1, "SeaState1");
SeaState1.Duration = 3600 s;
SeaState1.SpectrumItems.Add(new SpectrumItem());
SeaState1.SpectrumItems[0].Spectrum = BretschneiderSpectrum1;

HydroModel1 = new HydroModel(HydroModels, "HydroModel1");
HydroModel1.FPx = 1 m;

Panel = new ElementModel(Models, "Panel");
Panel.FileName = FEM_path+"T1.FEM";

Morison = new ElementModel(Models, "Morison");
Morison.FileName = FEM_path+"T6.FEM";

Compartment = new ElementModel(Models, "Compartment");
Compartment.FileName = FEM_path+"T5.FEM";

Free_surface = new ElementModel(Models, "Free_surface");
Free_surface.FileName = FEM_path+"T73.FEM";

Internal_lid = new ElementModel(Models, "Internal_lid");
Internal_lid.FileName = FEM_path+"T23.FEM";


PanelModel1 = new PanelModel(HydroModel1, "PanelModel1");
PanelModel1.ElementModel = Panel;
PanelModel1.SymmetryXZ = true;
PanelModel1.ShowHydroPressureArrows = true;


Morison2DProperties = new Morison2DProperties(Properties, "Morison2DProperties");

Col1 = new Morison2DProperty(Morison2DProperties, "Col1");
Col1.Cdy = 1;
Col1.Cdz = 1;
Col1.Cay = 0;
Col1.Caz = 0;

Col2 = new Morison2DProperty(Morison2DProperties, "Col2");
Col2.Cdy = 1;
Col2.Cdz = 1;
Col2.Cay = 0;
Col2.Caz = 0;

Pon = new Morison2DProperty(Morison2DProperties, "Pon");
Pon.Cdy = 1.68;
Pon.Cdz = 2.55;
Pon.Cay = 0;
Pon.Caz = 0;

if (anchor_formulation == "Anchor-element"){
    AnchorProperties1 = new AnchorProperties(Properties, "AnchorProperties1");
    AnchorProperty1 = new AnchorProperty(AnchorProperties1, "AnchorProperty1");
    AnchorProperty1.HorizontalStiffness = 0 N / m;
    AnchorProperty1.VerticalStiffness = 0 N / m;
    AnchorProperty1.PreTension = 1272440 N;
    AnchorProperty1.AngleSeaSurface = 36.24 deg;
}

LoadingConditions1 = new LoadingConditions(HydroModel1, "LoadingConditions1");
LoadingCondition1 = new LoadingCondition(LoadingConditions1, "LoadingCondition1");
LoadingCondition1.Location = Location1;
LoadingCondition1.WaterLineZ = 0 m;
LoadingCondition1.TrimAngle = 0 deg;
LoadingCondition1.HeelAngle = 0 deg;
LoadingCondition1.BalancingTolerance = new Fraction(0.0001);

MorisonModel1 = new MorisonModel(HydroModel1, "MorisonModel1");
MorisonModel1.ElementModel = Morison;
MorisonModel1.UpdateSections();
MorisonModel1.Sections[0].Morison2DProperty = Col1;
MorisonModel1.Sections[1].Morison2DProperty = Pon;
MorisonModel1.Sections[2].Morison2DProperty = Col2;
MorisonModel1.Sections[0].DragOnly = true;
MorisonModel1.Sections[1].DragOnly = true;
MorisonModel1.Sections[2].DragOnly = true;

if (anchor_formulation == "Anchor-element"){
    AnchorElement1 = new AnchorElement(MorisonModel1, "AnchorElement1");
    AnchorElement1.FairleadNode = 2;
    AnchorElement1.WindlassNode = 2;
    AnchorElement1.AnchorProperty = AnchorProperty1;
    AnchorElement1.XaxisAngle = 240 deg;

    AnchorElement2 = new AnchorElement(MorisonModel1, "AnchorElement2");
    AnchorElement2.FairleadNode = 14;
    AnchorElement2.WindlassNode = 14;
    AnchorElement2.AnchorProperty = AnchorProperty1;
    AnchorElement2.XaxisAngle = 0 deg;

    AnchorElement3 = new AnchorElement(MorisonModel1, "AnchorElement3");
    AnchorElement3.FairleadNode = 21;
    AnchorElement3.WindlassNode = 21;
    AnchorElement3.XaxisAngle = 120 deg;
    AnchorElement3.AnchorProperty = AnchorProperty1;
}
StructureReductions1 = new StructureReductions(Properties, "StructureReductions1");
StructureReduction1 = new StructureReduction(StructureReductions1, "StructureReduction1");
StructureReduction1.Value = new Fraction(0);
Permeabilities1 = new Permeabilities(Properties, "Permeabilities1");
Permeability1 = new Permeability(Permeabilities1, "Permeability1");
Permeability1.Value = new Fraction(1);

Contents1 = new Contents(Properties, "Contents1");
Seawater = new Fluid(Contents1, "Seawater");
Seawater.Density = 1025 kg/m^3;
Seawater.Color = Color.FromArgb(255, 0, 191, 255);
Iron_ore_concrete = new Fluid(Contents1, "Iron_ore_concrete");
Iron_ore_concrete.Density = 2560 kg/m^3;
Iron_ore_concrete.Color = Color.FromArgb(255, 255, 165, 0);

FillingFractions1 = new FillingFractions(Properties, "FillingFractions1");
Fluid_ballast = new FillingFraction(FillingFractions1, "Fluid_ballast");
Fluid_ballast.Value = new Fraction(0);
Fixed_ballast = new FillingFraction(FillingFractions1, "Fixed_ballast");
Fixed_ballast.Value = new Fraction(0);

CompartmentModel1 = new CompartmentModel(HydroModel1, "CompartmentModel1");
CompartmentModel1.FillingFractions.Add(new FillingFractionInterval());
CompartmentModel1.FillingFractions[0].ToFraction = new Fraction(1);
CompartmentModel1.FillingFractions[0].Step = new Fraction(0.1);
CompartmentModel1.ElementModel = Compartment;
CompartmentModel1.UpdateCompartments();
CompartmentModel1.Compartments[41].Group = "P1";
CompartmentModel1.Compartments[42].Group = "P1";
CompartmentModel1.Compartments[43].Group = "P1";
CompartmentModel1.Compartments[44].Group = "P1";

CompartmentModel1.Compartments[2].Group = "P2";
CompartmentModel1.Compartments[7].Group = "P2";
CompartmentModel1.Compartments[8].Group = "P2";
CompartmentModel1.Compartments[11].Group = "P2";

CompartmentModel1.Compartments[6].Group = "P3";
CompartmentModel1.Compartments[9].Group = "P3";
CompartmentModel1.Compartments[10].Group = "P3";
CompartmentModel1.Compartments[15].Group = "P3";

CompartmentModel1.Compartments[21].Group = "C";
CompartmentModel1.Compartments[23].Group = "C";
CompartmentModel1.Compartments[53].Group = "C";

CompartmentModel1.UpdateCompartments();

MassModel1 = new MassModel(HydroModel1, "MassModel1");
if (Mass_type=="Including_full_ballast") {
    MassModel1.UserSpecifiedMass = 20448.18 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = -2.0820 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 43.696 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 43.602 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 25.777 m;
    Fluid_ballast.Value = new Fraction(0);
    Fixed_ballast.Value = new Fraction(0);
} else if (Mass_type=="Including_fixed_ballast") {
    MassModel1.UserSpecifiedMass = 9224.29 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = 15.78 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 61.827 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 61.68 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 38.379 m;
    Fluid_ballast.Value = new Fraction(0.9261296);
    Fixed_ballast.Value = new Fraction(0);
} else if (Mass_type=="Excluding_ballast") {
    MassModel1.UserSpecifiedMass = 6684.29 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = 28.86 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 62.521 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 68.676 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 31.859 m;
    Fluid_ballast.Value = new Fraction(0.9261296);
    Fixed_ballast.Value = new Fraction(0.1591847);
}

if (anchor_formulation == "Point-mass") {
    MassModel2 = new MassModel(HydroModel1, "MassModel2");
    MassModel2.MassModelType = MassModelType.PointMassTable;
    MassModel2.PointMasses.Add(new PointMass());
    MassModel2.PointMasses[0].Mass = 76697.4451 kg;
    MassModel2.PointMasses[0].X = 58.0 m;
    MassModel2.PointMasses[0].Y = 0.0 m;
    MassModel2.PointMasses[0].Z = -14 m;
    MassModel2.PointMasses.Add(new PointMass());
    MassModel2.PointMasses[1].Mass = 76697.4451 kg;
    MassModel2.PointMasses[1].X = -29.0 m;
    MassModel2.PointMasses[1].Y = -50.229 m;
    MassModel2.PointMasses[1].Z = -14 m;
    MassModel2.PointMasses.Add(new PointMass());
    MassModel2.PointMasses[2].Mass = 76697.4451 kg;
    MassModel2.PointMasses[2].X = -29.0 m;
    MassModel2.PointMasses[2].Y = 50.229 m;
    MassModel2.PointMasses[2].Z = -14 m;
}

SecondOrderFreeSurfaceModel1 = new SecondOrderFreeSurfaceModel(LoadingCondition1, "SecondOrderFreeSurfaceModel1");
SecondOrderFreeSurfaceModel1.ElementModel = Free_surface;

CompartmentContents1 = new CompartmentContents(LoadingCondition1, "CompartmentContents1");
CompartmentModel1.UpdateCompartments();

CompartmentContents1.Contents[44].IntactFluid = Seawater;
CompartmentContents1.Contents[42].IntactFluid = Seawater;
CompartmentContents1.Contents[41].IntactFluid = Seawater;
CompartmentContents1.Contents[43].IntactFluid = Seawater;
CompartmentContents1.Contents[7].IntactFluid = Seawater;
CompartmentContents1.Contents[8].IntactFluid = Seawater;
CompartmentContents1.Contents[2].IntactFluid = Seawater;
CompartmentContents1.Contents[11].IntactFluid = Seawater;
CompartmentContents1.Contents[15].IntactFluid = Seawater;
CompartmentContents1.Contents[6].IntactFluid = Seawater;
CompartmentContents1.Contents[9].IntactFluid = Seawater;
CompartmentContents1.Contents[10].IntactFluid = Seawater;

CompartmentContents1.Contents[23].IntactFluid = Iron_ore_concrete;
CompartmentContents1.Contents[21].IntactFluid = Iron_ore_concrete;
CompartmentContents1.Contents[53].IntactFluid = Iron_ore_concrete;



CompartmentContents1.Contents[23].FillingFraction = Fixed_ballast;
CompartmentContents1.Contents[21].FillingFraction = Fixed_ballast;
CompartmentContents1.Contents[53].FillingFraction = Fixed_ballast;

CompartmentContents1.Contents[44].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[42].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[41].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[43].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[7].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[8].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[2].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[11].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[15].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[6].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[9].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[10].FillingFraction = Fluid_ballast;


AdditionalMatrices1 = new AdditionalMatrices(LoadingCondition1, "AdditionalMatrices1");
if (Drag_type == "Critical"){
    CriticalDampingMatrix1 = new CriticalDampingMatrix(AdditionalMatrices1, "CriticalDampingMatrix1");
    CriticalDampingMatrix1.Heave = 0.05;
    CriticalDampingMatrix1.Roll = 0.05;
    CriticalDampingMatrix1.Pitch = 0.05;
}

RestoringMatrix1 = new RestoringMatrix(AdditionalMatrices1, "RestoringMatrix1");
RestoringMatrix1.RestoringMatrixTable[0,0] = 79300.7933 N/m;
RestoringMatrix1.RestoringMatrixTable[1,1] = 79315.33389 N/m;
RestoringMatrix1.RestoringMatrixTable[5,5] = 196590315.3 N*m;

FixedLids1 = new FixedLids(LoadingCondition1, "FixedLids1");
FixedLid1 = new FixedLid(FixedLids1, "FixedLid1");
FixedLid1.ElementModel = Internal_lid;

WadamAnalysis1 = new WadamAnalysis(Analyses, "WadamAnalysis1");
WadamAnalysis1.DataPreparation = true;
WadamAnalysis1.ExecuteAnalysis = true;
WadamAnalysis1.ToleranceWaterline = new Fraction(0.05);
WadamAnalysis1.ToleranceCOG = new Fraction(0.05);
WadamAnalysis1.CharacteristicLength = 100 m;
WadamAnalysis1.AnnulusSize = 500 m;
WadamAnalysis1.UseFreeSurfaceIntegral = true;
WadamAnalysis1.TranslationalConvergenceCriteria = new Fraction(0.001);
WadamAnalysis1.RotationalConvergenceCriteria = new Fraction(0.001);
WadamAnalysis1.MaxLinearizationIterations = 20;
WadamAnalysis1.UseStochasticLinearization = true;
WadamAnalysis1.MaxMatrixSize = 50000;
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems.Add(new MaximumRollAngleItem());
WadamAnalysis1.MaximumRollAngleItems[0].MaxRollAngle = 0.1 deg;
WadamAnalysis1.MaximumRollAngleItems[1].Heading = 45 deg;
WadamAnalysis1.MaximumRollAngleItems[1].MaxRollAngle = 2.5 deg;
WadamAnalysis1.MaximumRollAngleItems[2].Heading = 90 deg;
WadamAnalysis1.MaximumRollAngleItems[2].MaxRollAngle = 4 deg;
WadamAnalysis1.MaximumRollAngleItems[3].Heading = 135 deg;
WadamAnalysis1.MaximumRollAngleItems[3].MaxRollAngle = 2.5 deg;
WadamAnalysis1.MaximumRollAngleItems[4].Heading = 180 deg;
WadamAnalysis1.MaximumRollAngleItems[4].MaxRollAngle = 0.1 deg;
WadamAnalysis1.MaximumRollAngleItems[5].Heading = 225 deg;
WadamAnalysis1.MaximumRollAngleItems[5].MaxRollAngle = 2.5 deg;
WadamAnalysis1.MaximumRollAngleItems[6].Heading = 270 deg;
WadamAnalysis1.MaximumRollAngleItems[6].MaxRollAngle = 4 deg;
WadamAnalysis1.MaximumRollAngleItems[7].Heading = 315 deg;
WadamAnalysis1.MaximumRollAngleItems[7].MaxRollAngle = 2.5 deg;
WadamAnalysis1.MaximumRollAngleItems[8].Heading = 360 deg;
WadamAnalysis1.MaximumRollAngleItems[8].MaxRollAngle = 0.1 deg;
WadamAnalysis1.CalculateEigenvalues = true;
WadamAnalysis1.PrintType = PrintType.NormalPrint;
WadamAnalysis1.PanelPlateOutOfPlaneTolerance = new Fraction(0.8);
WadamAnalysis1.PanelToPlateAngularDifference = 20 deg;
WadamAnalysis1.LoadTransferPressureReductionZoneExtension = 1 m;
WadamAnalysis1.LoadTransferIncludeStaticLoad = true;
WadamAnalysis1.LoadTransferIncludeStaticGravity = true;
WadamAnalysis1.LoadTransferIncludeInertia = true;
WadamAnalysis1.LoadTransferIncludeFluctuatingGravity = true;
WadamAnalysis1.LoadTransferIncludeFluctuatingPressure = true;
WadamAnalysis1.LoadTransferIncludeMooringLoad = true;
WadamAnalysis1.CalculateDriftForces = true;
WadamAnalysis1.FarFieldIntegration = true;
WadamAnalysis1.UseRollGZCurve = true;
WadamAnalysis1.SaveRestartFileName = "WADAM";
WadamAnalysis1.HydroModel = HydroModel1;
WadamAnalysis1.LoadingCondition = LoadingCondition1;
WadamAnalysis1.UseStochasticLinearization = false;
WadamAnalysis1.UseRollGZCurve = false;
WadamAnalysis1.UseMorisonModel = true;
WadamAnalysis1.RemoveIrregularFrequencies = true;
WadamAnalysis1.EnvironmentCondition = FrequencyDomainCondition1;
WadamAnalysis1.SaveWamitFiles = true;

if (Drag_type == "Critical"){
    WadamAnalysis1.UseMorisonModel = false;
} else if (Drag_type == "Morison"){
    WadamAnalysis1.UseMorisonModel = true;
    WadamAnalysis1.SeaState = SeaState1;
    WadamAnalysis1.DragMethod = DragMethod.Stochastic;
    WadamAnalysis1.WaveType = WaveType.IncidentWave;
}

if (mean_drift){
    WadamAnalysis1.CalculateDriftForces = true;
    WadamAnalysis1.IncludeBiDirectionalWaves = false;
    WadamAnalysis1.PressureIntegration = true;
    WadamAnalysis1.PressureIntegrationControlSurfaceType = PressureIntegrationControlSurfaceType.Circular;
    WadamAnalysis1.FarFieldIntegration = true;
}else{
    WadamAnalysis1.CalculateDriftForces = false;
}

if (diff_freq){
    WadamAnalysis1.DifferenceFrequencies = true;
    WadamAnalysis1.CombiDiffType = CombinationDiffType.Selected;
    WadamAnalysis1.SecondOrderQuadraticForce = true;
    WadamAnalysis1.MaxFreqDiff = 0.25 rad/s;
    WadamAnalysis1.MaxDirDiff = 0 deg;
    WadamAnalysis1.UseFreeSurfaceIntegral = true;
    WadamAnalysis1.UseThreeLayerFreeSurface = true;
    if (Drag_type == "Morison") {
        WadamAnalysis1.IncludeDampingFromMorisonModel = true;
    }
}else{
    WadamAnalysis1.DifferenceFrequencies = false;
}

if (tank_method == "Dynamic"){
    WadamAnalysis1.IncludeDynamicsOfInternalFluid = true;
    WadamAnalysis1.RemoveIrregularFrequencies = false;
} else if (tank_method == "Quasi-static"){
    WadamAnalysis1.IncludeDynamicsOfInternalFluid = false;
}


WadamAnalysis1.IncludeLimitingFrequencies = true;
WadamAnalysis1.IncludeForwardSpeed = false;

if (execute){
    WadamAnalysis1.Execute();
}