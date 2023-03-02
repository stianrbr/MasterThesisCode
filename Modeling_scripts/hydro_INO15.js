// HydroD V6.1-02 started 2023-02-27 09:44:57


FEM_path = "S:\\Master\\Linear convergence study\\INO\\Mesh_1_m\\Geo_1m\\"
Mass_type = "Including_full_ballast"



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

Col1 = new Morison2DProperty(Morison2DProperties, "Col");
Col1.Cdy = 1;
Col1.Cdz = 1;
Col1.Cay = 1;
Col1.Caz = 1;


Pon = new Morison2DProperty(Morison2DProperties, "Pon");
Pon.Cdy = 1;
Pon.Cdz = 1;
Pon.Cay = 1;
Pon.Caz = 1;

AnchorProperties1 = new AnchorProperties(Properties, "AnchorProperties1");
AnchorProperty1 = new AnchorProperty(AnchorProperties1, "AnchorProperty1");
AnchorProperty1.HorizontalStiffness = 0 N/m;
AnchorProperty1.VerticalStiffness = 0 N/m;
AnchorProperty1.PreTension = 1272440 N;
AnchorProperty1.AngleSeaSurface = 36.24 deg;


LoadingConditions1 = new LoadingConditions(HydroModel1, "LoadingConditions1");
LoadingCondition1 = new LoadingCondition(LoadingConditions1, "LoadingCondition1");
LoadingCondition1.Location = Location1;
LoadingCondition1.WaterLineZ = 0 m;
LoadingCondition1.TrimAngle = 0 deg;
LoadingCondition1.HeelAngle = 0 deg;

MorisonModel1 = new MorisonModel(HydroModel1, "MorisonModel1");
MorisonModel1.ElementModel = Morison;
MorisonModel1.UpdateSections();
MorisonModel1.Sections[0].Morison2DProperty = Col;
MorisonModel1.Sections[1].Morison2DProperty = Pon;
MorisonModel1.Sections[0].DragOnly = true;
MorisonModel1.Sections[1].DragOnly = true;

AnchorElement1 = new AnchorElement(MorisonModel1, "AnchorElement1");
AnchorElement1.FairleadNode = 2;
AnchorElement1.WindlassNode = 2;
AnchorElement1.AnchorProperty = AnchorProperty1;
AnchorElement1.XaxisAngle = 0 deg;

AnchorElement2 = new AnchorElement(MorisonModel1, "AnchorElement2");
AnchorElement2.FairleadNode = 9;
AnchorElement2.WindlassNode = 9;
AnchorElement2.AnchorProperty = AnchorProperty1;
AnchorElement2.XaxisAngle = 120 deg;

AnchorElement3 = new AnchorElement(MorisonModel1, "AnchorElement3");
AnchorElement3.FairleadNode = 16;
AnchorElement3.WindlassNode = 16;
AnchorElement3.XaxisAngle = 240 deg;
AnchorElement3.AnchorProperty = AnchorProperty1;

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
CompartmentModel1.Compartments[3].Group = "P1";
CompartmentModel1.Compartments[4].Group = "P1";
CompartmentModel1.Compartments[33].Group = "P1";
CompartmentModel1.Compartments[35].Group = "P1";

CompartmentModel1.Compartments[1].Group = "P2";
CompartmentModel1.Compartments[2].Group = "P2";
CompartmentModel1.Compartments[34].Group = "P2";
CompartmentModel1.Compartments[36].Group = "P2";

CompartmentModel1.Compartments[11].Group = "P3";
CompartmentModel1.Compartments[8].Group = "P3";
CompartmentModel1.Compartments[25].Group = "P3";
CompartmentModel1.Compartments[23].Group = "P3";

CompartmentModel1.Compartments[19].Group = "C";
CompartmentModel1.Compartments[20].Group = "C";
CompartmentModel1.Compartments[41].Group = "C";

CompartmentModel1.UpdateCompartments();

MassModel1 = new MassModel(HydroModel1, "MassModel1");
if (Mass_type=="Including_full_ballast") {
    MassModel1.UserSpecifiedMass = 13244.78 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = 7.5579 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 55.507 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 55.510 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 41.336 m;
    Fluid_ballast.Value = new Fraction(0);
    Fixed_ballast.Value = new Fraction(0);
} else if (Mass_type=="Including_fixed_ballast") {
    MassModel1.UserSpecifiedMass = 12458.98 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = 9.28 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 59.971 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 62.295 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 45.829 m;
    Fluid_ballast.Value = new Fraction(0.44385);
    Fixed_ballast.Value = new Fraction(0);
} else if (Mass_type=="Excluding_ballast") {
    MassModel1.UserSpecifiedMass = 7097.64 tonne;
    MassModel1.UserSpecifiedCenterOfGravityX = 16.86 m;
    MassModel1.UserSpecifiedCenterOfGravityY = 0 m;
    MassModel1.UserSpecifiedCenterOfGravityZ = 28.56 m;
    MassModel1.UserSpecifiedRadiusOfGyrationX = 60.596 m;
    MassModel1.UserSpecifiedRadiusOfGyrationY = 72.6380 m;
    MassModel1.UserSpecifiedRadiusOfGyrationZ = 41.070 m;
    Fluid_ballast.Value = new Fraction(0.44385);
    Fixed_ballast.Value = new Fraction(0.159153);
}


SecondOrderFreeSurfaceModel1 = new SecondOrderFreeSurfaceModel(LoadingCondition1, "SecondOrderFreeSurfaceModel1");
SecondOrderFreeSurfaceModel1.ElementModel = Free_surface;

CompartmentContents1 = new CompartmentContents(LoadingCondition1, "CompartmentContents1");
CompartmentModel1.UpdateCompartments();

CompartmentContents1.Contents[33].IntactFluid = Seawater;
CompartmentContents1.Contents[4].IntactFluid = Seawater;
CompartmentContents1.Contents[3].IntactFluid = Seawater;
CompartmentContents1.Contents[35].IntactFluid = Seawater;
CompartmentContents1.Contents[36].IntactFluid = Seawater;
CompartmentContents1.Contents[34].IntactFluid = Seawater;
CompartmentContents1.Contents[2].IntactFluid = Seawater;
CompartmentContents1.Contents[1].IntactFluid = Seawater;
CompartmentContents1.Contents[25].IntactFluid = Seawater;
CompartmentContents1.Contents[11].IntactFluid = Seawater;
CompartmentContents1.Contents[8].IntactFluid = Seawater;
CompartmentContents1.Contents[23].IntactFluid = Seawater;


CompartmentContents1.Contents[20].IntactFluid = Iron_ore_concrete;
CompartmentContents1.Contents[19].IntactFluid = Iron_ore_concrete;
CompartmentContents1.Contents[41].IntactFluid = Iron_ore_concrete;



CompartmentContents1.Contents[20].FillingFraction = Fixed_ballast;
CompartmentContents1.Contents[19].FillingFraction = Fixed_ballast;
CompartmentContents1.Contents[41].FillingFraction = Fixed_ballast;

CompartmentContents1.Contents[33].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[4].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[3].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[35].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[36].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[34].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[2].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[1].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[25].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[11].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[8].FillingFraction = Fluid_ballast;
CompartmentContents1.Contents[23].FillingFraction = Fluid_ballast;


AdditionalMatrices1 = new AdditionalMatrices(LoadingCondition1, "AdditionalMatrices1");
CriticalDampingMatrix1 = new CriticalDampingMatrix(AdditionalMatrices1, "CriticalDampingMatrix1");
CriticalDampingMatrix1.Heave = 0.05;
CriticalDampingMatrix1.Roll = 0.05;
CriticalDampingMatrix1.Pitch = 0.05;

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

WadamAnalysis1.UseMorisonModel = false;
WadamAnalysis1.CalculateDriftForces = false;
WadamAnalysis1.DifferenceFrequencies = false;
WadamAnalysis1.IncludeLimitingFrequencies = true;
WadamAnalysis1.IncludeForwardSpeed = false;