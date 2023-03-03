
path = "C:\\Users\\stian\\OneDrive - NTNU\\5. klasse\\Masteroppgave\\MasterThesisCode\\Modeling_scripts\\"

vessel = "VolturnUS-S"  // "VolturnUS-S" / "INO15"

geom_file = path+"geometry_"+vessel+".js"

hydro_file = path+"hydro_"+vessel+".js"


job = new Job();

GeniEActivity1 = new GeniEActivity();

this.Activities.Add(GeniEActivity1);

HydroDActivity1 = new HydroDActivity();

this.Activities.Add(HydroDActivity1);

GeniEActivity1.CmdInputFile = geom_file;
Copy(geom_file, GeniEActivity1.Workspace + "geometry_"+vessel+".js");
GeniEActivity1.CmdInputFile = GeniEActivity1.Workspace + "geometry_"+vessel+".js"

GeniEActivity1.UseLocalCmdInputFile = true;
GeniEActivity1.InputMode = InputMode.Background;

GeniEActivity1.SelectedLicenses.Remove("RefineMesh");

HydroDActivity1.CmdInputFile = hydro_file;
HydroDActivity1.UseLocalCmdInputFile = true;
Copy(hydro_file, HydroDActivity1.Workspace + "hydro_"+vessel+".js");
HydroDActivity1.CmdInputFile = HydroDActivity1.Workspace + "hydro_"+vessel+".js"


this.Execute();

