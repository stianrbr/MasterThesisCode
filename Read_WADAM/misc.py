import re
import numpy as np
import warnings
from io import StringIO
regexp = "[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?"

class Radiation_res:
    def __init__(self,matrix=None, period=None):
        self.Matrix = matrix
        self.period = period

class Dimensions:
    def __init__(self, rho=None, g=None, vol=None, l=None, wa=None):
        self.rho = rho
        self.g = g
        self.vol = vol
        self.l  = l
        self.wa = wa

class Point:
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z

class Mass_matrix:
    def __init__(self, M, XG, YG, ZG, XRAD, YRAD, ZRAD, XYRAD, XZRAD, YZRAD):
        self.M = M
        self.XG, self.YG, self.ZG = XG, YG, ZG
        self.COG = Point(XRAD, YRAD, ZRAD)
        self.XYRAD, self.XZRAD, self.YZRAD = XYRAD, XZRAD, YZRAD


class Hydro_static:
    def __init__(self, vol_disp, t_moor_vertical, wp_area, xb, yb, zb):
        self.vol_disp = vol_disp
        self.t_moor_vertical = t_moor_vertical
        self.wp_area = wp_area
        self.COB = Point(xb, yb, zb)

class Environmental_data:
    def __init__(self, water_depth, wave_periods, wave_headings):
        self.water_depth = water_depth
        self.wave_periods = wave_periods
        self.wave_headings = wave_headings

class Tank:
    def __init__(self, tank_numb=None, mass=None, cog=None, added_mass=None, damping=None, restoring=None):
        self.tank_numb = tank_numb
        self.mass = mass
        self.cog = cog
        self.added_mass = added_mass
        self.damping = damping
        self.restoring = restoring

class RAO:

    def __init__(self, amp, ph=None, period=None, heading=None):
        self.amplitude = amp
        self.phase = ph
        self.period = period
        self.heading = heading

class RAO_2nd_order:

    def __init__(self, amp, ph, pi, pj, hi, hj):
        self.amplitude = amp
        self.phase = ph
        self.period_i = pi
        self.period_j = pj
        self.heading_i = hi
        self.heading_j = hj

class Dimensionalizing_factors:
    def __init__(self, ro, g, vol, l, wa):
        self.RO = ro
        self.G = g
        self.VOL = vol
        self.L = l
        self.WA = wa


class WADAM_res:
    def __init__(self, filename):
        self.filename = filename
        self.compartment_flag = False
        self.mooring_flag = False
        # Added mass
        self.Total_added_mass = []
        self.Hull_added_mass = None
        self.Zero_freq_added_mass = None
        self.Inf_freq_added_mass = None
        # Damping
        self.Total_damping = []
        self.Total_pot_damping = None
        self.Hull_damping = None
        # Restoring
        self.Hydrostat_restoring = None
        self.Mooring_restoring = None
        self.Total_restoring = None
        # Load RAO
        self.LoadRAO = []
        # Response RAO
        self.MotionRAO = []
        # Tanks
        self.tanks = None
        # Mean drift (far-field)
        self.mean_drift_horizontal = []
        # Mean drift (PI/CS)
        self.mean_drift_full = []
        # QTF
        self.diff_freq_load_QTF = []
        self.diff_freq_motion_QTF = []
        # Factors for dimensionalizing
        self.dimensions = None

        self.read_file()
        self.dimensionalize()

    def read_file(self):
        """

        :return:
        """
        """
        Extracting data from LIS-file
        """
        with open(self.filename) as f:
            self.datafile = f.readlines()
            for i, line in enumerate(self.datafile):
                if not line.isspace():
                    if line == '1\n':
                        if self.datafile[i+1] == ' 4.3 GLOBAL HYDRODYNAMIC RESULTS\n':
                            page_start = i+3
                        else:
                            page_start = i
                    if "COMPARTMENT" in line:
                        self.compartment_flag = True
                        self.tanks = []
                        self.Hull_added_mass = []
                        self.Hull_damping = []
                        self.Total_pot_damping = []

                    if "STATIC RESTORING COEFFICIENT MATRIX FOR TLP AND MOORING ELEMENTS" in line:
                        self.mooring_flag = True

                    if "NUMBER OF BASIC PANELS         " in line:
                        self.Num_panels = int(re.findall(regexp, line)[0])

                    if "MASS PROPERTIES AND STRUCTURAL DATA:" in line:
                        try:
                            section = self.datafile[i:i+20]
                            m = float(re.findall(regexp,[s for s in section if "MASS OF THE STRUCTURE                      M" in s][0])[0])
                            xg = float(re.findall(regexp,[s for s in section if "CENTRE OF GRAVITY                          XG" in s][0])[0])
                            yg = float(re.findall(regexp, [s for s in section if "                                           YG" in s][0])[0])
                            zg = float(re.findall(regexp, [s for s in section if "                                           ZG" in s][0])[0])
                            xrad = float(re.findall(regexp, [s for s in section if "ROLL  RADIUS OF GYRATION                   XRAD" in s][0])[0])
                            yrad = float(re.findall(regexp, [s for s in section if "PITCH RADIUS OF GYRATION                   YRAD" in s][0])[0])
                            zrad = float(re.findall(regexp, [s for s in section if "YAW   RADIUS OF GYRATION                   ZRAD" in s][0])[0])
                            xyrad = float(re.findall(regexp, [s for s in section if "ROLL-PITCH CENTRIFUGAL MOMENT              XYRAD" in s][0])[0])
                            xzrad = float(re.findall(regexp, [s for s in section if "ROLL-YAW   CENTRIFUGAL MOMENT              XZRAD" in s][0])[0])
                            yzrad = float(re.findall(regexp, [s for s in section if "PITCH-YAW  CENTRIFUGAL MOMENT              YZRAD" in s][0])[0])
                            self.Massmatrix = Mass_matrix(m, xg, yg, zg, xrad, yrad, zrad, xyrad, xzrad, yzrad)
                        except:
                            raise Exception("Failed reading mass properties")
                    elif "HYDROSTATIC DATA:" in line:
                        try:
                            section = self.datafile[i:i + 20]
                            vol_disp = float(re.findall(regexp,[s for s in section if "DISPLACED VOLUME                           VOL" in s][0])[0])
                            try:
                                t_moor_vertical = float(re.findall(regexp,[s for s in section if "VERTICAL PRETENSION OF MOORING ELEMENTS" in s][0])[0])
                            except:
                                t_moor_vertical = None
                                warnings.warn("No vertical mooring tension found.")
                            wp_area = float(re.findall(regexp,[s for s in section if "WATER PLANE AREA                           WPLA" in s][0])[0])
                            xb = float(re.findall(regexp,[s for s in section if "CENTRE OF BUOYANCY                         XCB" in s][0])[0])
                            yb = float(re.findall(regexp,[s for s in section if "                                           YCB" in s][0])[0])
                            zb = float(re.findall(regexp,[s for s in section if "                                           ZCB" in s][0])[0])
                            self.Hydro_static = Hydro_static(vol_disp, t_moor_vertical, wp_area,  xb, yb, zb)
                        except:
                            raise Exception("Failed reading hydrostatic properties")
                    elif "2.8 ENVIRONMENTAL DATA:" in line:
                        try:
                            section = self.datafile[i:i+10]
                            water_depth = float(re.findall(regexp,[s for s in section if "WATER DEPTH                                          =" in s][0])[0])
                            n_wave_per = int(re.findall(regexp,[s for s in section if "NUMBER OF WAVE LENGTHS                               =" in s][0])[0])
                            n_wave_head = int(re.findall(regexp,[s for s in section if "NUMBER OF HEADING ANGLES                             =" in s][0])[0])

                        except:
                            raise Exception("Failed reading environmental data")
                    elif "WAVE DESCRIPTION:" in line:
                        try:
                            section = self.datafile[i+5:i+n_wave_per+5]
                            wave_per = [float(re.findall(regexp, w)[3]) for w in section]
                        except:
                            raise Exception("Failed reading wave periods")
                    elif "HEADING ANGLES (ANGLE BETWEEN POS. X-AXIS AND DIRECTION" in line:
                        try:
                            section = self.datafile[i+5:i+n_wave_head+5]
                            wave_head = [float(re.findall(regexp, w)[1]) for w in section]
                            self.Environment_data = Environmental_data(water_depth, np.array(wave_per), np.array(wave_head))
                        except:
                            raise Exception("Failed reading wave headings")

                    """
                    Extracting restoring data
                    """

                    if "HYDROSTATIC RESTORING COEFFICIENT MATRIX" in line:
                        c_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(np.longdouble)
                        self.Hydrostat_restoring = c_matrix
                    if "TOTAL RESTORING COEFFICIENT MATRIX" in line:
                        c_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(np.longdouble)
                        self.Total_restoring = c_matrix
                    if self.compartment_flag:
                        if "STATIC RESTORING COEFFICIENT MATRIX FOR FLUID IN TANK" in line:
                            self.tank_numb = int(re.findall(regexp, line)[0])
                            self.c_t_matrix = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(np.longdouble)
                            self.tanks.append(Tank(tank_numb=self.tank_numb, added_mass=[], damping=[], restoring=self.c_t_matrix))
                    if self.mooring_flag:
                        if "STATIC RESTORING COEFFICIENT MATRIX FOR TLP AND MOORING ELEMENTS" in line:
                            self.Mooring_restoring = np.array([re.findall(regexp, self.datafile[i+4+c_i])[1:] for c_i in range(0,6)]).astype(np.longdouble)

                    """
                    Extracting added mass data
                    """

                    if "ADDED MASS MATRIX             " in line:
                        a_matrix = np.array([re.findall(regexp, self.datafile[i+4+a_i])[1:] for a_i in range(0,6)]).astype(np.longdouble)
                        wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                        self.Total_added_mass.append(Radiation_res(matrix=a_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))

                    if "POTENTIAL ADDED MASS MATRIX FOR INFINITE WAVE PERIOD" in line:
                        self.Inf_freq_added_mass = np.array([re.findall(regexp, self.datafile[i+4+a_i])[1:] for a_i in range(0,6)]).astype(np.longdouble)

                    if "POTENTIAL ADDED MASS MATRIX FOR ZERO WAVE PERIOD" in line:
                        self.Zero_freq_added_mass = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(np.longdouble)

                    if self.compartment_flag:
                        if "ADDED MASS MATRIX OF OUTER HULL" in line:
                            a_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            self.Hull_added_mass.append(Radiation_res(matrix=a_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))
                        if "ADDED MASS MATRIX OF TANK " in line:
                            tank_numb = int(re.findall(regexp, line)[0])
                            a_t_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + a_i])[1:] for a_i in range(0, 6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            for t in list(self.tanks):
                                if t.tank_numb == tank_numb:
                                    temp = Radiation_res(matrix=a_t_matrix, period=self.Environment_data.wave_periods[wave_len_number-1])
                                    t.added_mass.append(temp)
                                    break



                    """
                    Extracting damping data
                    """
                    if self.compartment_flag:
                        if "TOTAL POTENTIAL DAMPING MATRIX" in line:
                            b_matrix = np.array([re.findall(regexp, self.datafile[i+4+b_i])[1:] for b_i in range(0,6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            self.Total_pot_damping.append(Radiation_res(matrix=b_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))
                        if "POTENTIAL DAMPING MATRIX OF OUTER HULL" in line:
                            b_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + b_i])[1:] for b_i in range(0, 6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            self.Hull_damping.append(Radiation_res(matrix=b_matrix,period=self.Environment_data.wave_periods[wave_len_number - 1]))
                        if "POTENTIAL DAMPING MATRIX OF TANK " in line:
                            tank_numb = int(re.findall(regexp, line)[0])
                            b_t_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + b_i])[1:] for b_i in range(0, 6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            for t in self.tanks:
                                if t.tank_numb == tank_numb:
                                    t.damping.append(Radiation_res(matrix=b_t_matrix, period=self.Environment_data.wave_periods[wave_len_number-1]))
                    else:
                        if "TOTAL DAMPING MATRIX" in line:
                            b_matrix = np.array([re.findall(regexp, self.datafile[i + 4 + b_i])[1:] for b_i in range(0, 6)]).astype(np.longdouble)
                            wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                            self.Total_damping.append(Radiation_res(matrix=b_matrix, period=self.Environment_data.wave_periods[wave_len_number - 1]))

                    """
                    Extracting 1st order load RAO and motion RAO
                    """
                    if "EXCITING FORCES AND MOMENTS" in line:
                        wave_len_number = int(re.findall(regexp, self.datafile[page_start+2])[1])
                        heading_number = int(re.findall(regexp, self.datafile[page_start+4])[0])
                        amp_load = np.array([re.findall(regexp, self.datafile[i + 5 + 2*f_i])[3] for f_i in range(0,6)]).astype(np.longdouble)
                        ph_load = np.array([re.findall(regexp, self.datafile[i + 5 + 2 * f_i])[4] for f_i in range(0, 6)]).astype(np.longdouble)
                        if len(self.LoadRAO) < heading_number:
                            self.LoadRAO.append([])
                        self.LoadRAO[heading_number - 1].append(RAO(amp=amp_load, ph=ph_load, period=self.Environment_data.wave_periods[wave_len_number - 1], heading=self.Environment_data.wave_headings[heading_number - 1]))

                        assert("MOTION" in self.datafile[i+18])
                        amp_motion = np.array([re.findall(regexp, self.datafile[i + 23 + 2*x_i])[3] for x_i in range(0,6)]).astype(np.longdouble)
                        ph_motion = np.array([re.findall(regexp, self.datafile[i + 23 + 2 * x_i])[4] for x_i in range(0, 6)]).astype(np.longdouble)
                        if len(self.MotionRAO) < heading_number:
                            self.MotionRAO.append([])
                        self.MotionRAO[heading_number - 1].append(RAO(amp=amp_motion, ph=ph_motion, period=self.Environment_data.wave_periods[wave_len_number - 1], heading=self.Environment_data.wave_headings[heading_number - 1]))

                    if "HORIZONTAL MEAN DRIFT FORCES AND MOMENT" in line:
                        d_matrix = np.array([re.findall(regexp, self.datafile[i+4 + d_i])[0] for d_i in range(0,3)]).astype(np.longdouble)
                        wave_len_number = int(re.findall(regexp, self.datafile[page_start + 2])[1])
                        heading_number = int(re.findall(regexp, self.datafile[page_start + 4])[1])
                        if len(self.mean_drift_horizontal) < heading_number:
                            self.mean_drift_horizontal.append([])
                        self.mean_drift_horizontal[heading_number-1].append(RAO(amp=d_matrix, period=self.Environment_data.wave_periods[wave_len_number - 1], heading=self.Environment_data.wave_headings[heading_number-1]))

                    if "MEAN DRIFT FORCES AND MOMENTS" in line and not("HORIZONTAL" in line):
                        d_matrix = np.array([re.findall(regexp, self.datafile[i+4 + d_i])[0] for d_i in range(0,6)]).astype(np.longdouble)
                        wave_len_number = int(re.findall(regexp, self.datafile[page_start + 2])[1])
                        heading_number = int(re.findall(regexp, self.datafile[page_start + 4])[1])
                        if len(self.mean_drift_full) < heading_number:
                            self.mean_drift_full.append([])
                        self.mean_drift_full[heading_number-1].append(RAO(amp=d_matrix, period=self.Environment_data.wave_periods[wave_len_number - 1], heading=self.Environment_data.wave_headings[heading_number-1]))

                    if "QUADRATIC SECOND-ORDER DIFFERENCE-FREQUENCY FORCES AND MOMENTS" in line:
                        wave_len_number_i =    int(re.findall(regexp, self.datafile[page_start + 2])[1])
                        wave_len_number_j = int(re.findall(regexp, self.datafile[page_start + 2])[2])
                        heading_number_i = int(re.findall(regexp, self.datafile[page_start +4])[0])
                        heading_number_j = int(re.findall(regexp, self.datafile[page_start +4])[1])
                        amp_load = np.array([re.findall(regexp, self.datafile[i + 5 + 2 * f_i])[3] for f_i in range(0, 6)]).astype(np.longdouble)
                        ph_load = np.array([re.findall(regexp, self.datafile[i + 5 + 2 * f_i])[4] for f_i in range(0, 6)]).astype(np.longdouble)
                        if len(self.diff_freq_load_QTF) < heading_number_i:
                            self.diff_freq_load_QTF.append([])
                        if len(self.diff_freq_load_QTF[heading_number_i - 1]) < heading_number_j:
                            self.diff_freq_load_QTF[heading_number_i - 1].append([])
                        self.diff_freq_load_QTF[heading_number_i - 1][heading_number_j - 1].append(RAO_2nd_order(amp=amp_load, ph=ph_load, pi=self.Environment_data.wave_periods[wave_len_number_i - 1],
                                                                                                                 pj= self.Environment_data.wave_periods[wave_len_number_j-1],
                                                                                                                 hi= self.Environment_data.wave_headings[heading_number_i-1],
                                                                                                                 hj= self.Environment_data.wave_headings[heading_number_j-1]))
                        assert("QUADRATIC SECOND-ORDER DIFFERENCE-FREQUENCY MOTIONS" in self.datafile[i+19])
                        amp_motion = np.array([re.findall(regexp, self.datafile[i + 24 + 2 * x_i])[3] for x_i in range(0, 6)]).astype(np.longdouble)
                        ph_motion = np.array([re.findall(regexp, self.datafile[i + 24 + 2 * x_i])[4] for x_i in range(0, 6)]).astype(np.longdouble)
                        if len(self.diff_freq_motion_QTF) < heading_number_i:
                            self.diff_freq_motion_QTF.append([])
                        if len(self.diff_freq_motion_QTF[heading_number_i - 1]) < heading_number_j:
                            self.diff_freq_motion_QTF[heading_number_i - 1].append([])
                        self.diff_freq_motion_QTF[heading_number_i - 1][heading_number_j - 1].append(RAO_2nd_order(amp=amp_motion, ph=ph_motion, pi=self.Environment_data.wave_periods[wave_len_number_i - 1],
                                                                                                                 pj= self.Environment_data.wave_periods[wave_len_number_j-1],
                                                                                                                 hi= self.Environment_data.wave_headings[heading_number_i-1],
                                                                                                                 hj= self.Environment_data.wave_headings[heading_number_j-1]))

                    if "THE OUTPUT IS NON-DIMENSIONALIZED USING -" in line:
                        ro = float(re.findall(regexp, self.datafile[i+8])[0])
                        g = float(re.findall(regexp, self.datafile[i+9])[0])
                        vol = float(re.findall(regexp, self.datafile[i+10])[0])
                        l = float(re.findall(regexp, self.datafile[i+11])[0])
                        wa = float(re.findall(regexp, self.datafile[i+12])[0])
                        self.dimensions = Dimensionalizing_factors(ro, g, vol, l, wa)


    def dimensionalize(self):
        """

        :return:
        """
        """
        Establishing factors for re-dimensionalizing
        """
        added_mass_dimension = np.zeros((6, 6), dtype=np.longdouble)
        added_mass_dimension[0:3, 0:3] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL
        added_mass_dimension[3:6, 0:3] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.L
        added_mass_dimension[0:3, 3:6] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.L
        added_mass_dimension[3:6, 3:6] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.L ** 2

        damping_dimensions = np.zeros((6, 6), dtype=np.longdouble)
        damping_dimensions[0:3, 0:3] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * np.sqrt(self.dimensions.G / self.dimensions.L)
        damping_dimensions[3:6, 0:3] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * np.sqrt(self.dimensions.G * self.dimensions.L)
        damping_dimensions[0:3, 3:6] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * np.sqrt(self.dimensions.G * self.dimensions.L)
        damping_dimensions[3:6, 3:6] = np.ones((3, 3), dtype=np.longdouble) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.L * np.sqrt(self.dimensions.G * self.dimensions.L)

        restoring_dimensions = np.zeros((6, 6))
        restoring_dimensions[0:3, 0:3] = np.ones((3, 3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G / self.dimensions.L
        restoring_dimensions[3:6, 0:3] = np.ones((3, 3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G
        restoring_dimensions[0:3, 3:6] = np.ones((3, 3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G
        restoring_dimensions[3:6, 3:6] = np.ones((3, 3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G * self.dimensions.L

        excitation_dimensions = np.zeros((6))
        excitation_dimensions[0:3] = np.ones((3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G * self.dimensions.WA / self.dimensions.L
        excitation_dimensions[3:6] = np.ones((3)) * self.dimensions.RO * self.dimensions.VOL * self.dimensions.G * self.dimensions.WA

        motions_dimensions = np.zeros((6))
        motions_dimensions[0:3] = np.ones((3)) * self.dimensions.WA
        motions_dimensions[3:6] = np.ones((3)) * self.dimensions.WA / self.dimensions.L * 180 / np.pi

        mean_drift_dimensions = np.zeros((6))
        mean_drift_dimensions[0:3] = np.ones((3)) * self.dimensions.RO * self.dimensions.G * self.dimensions.L * self.dimensions.WA**2
        mean_drift_dimensions[3:6] = np.ones((3)) * self.dimensions.RO * self.dimensions.G * self.dimensions.L**2 * self.dimensions.WA**2

        QTF_force_dimensions = np.zeros((6))
        QTF_force_dimensions[0:3] = np.ones((3)) * self.dimensions.RO * self.dimensions.G * self.dimensions.L
        QTF_force_dimensions[3:6] = np.ones((3)) * self.dimensions.RO * self.dimensions.G * self.dimensions.L**2

        QTF_motion_dimensions = np.zeros((6))
        QTF_motion_dimensions[0:3] = np.ones((3)) * 1 / self.dimensions.L
        QTF_motion_dimensions[3:6] = np.ones((3)) * 1 / self.dimensions.L**2 *  180 / np.pi

        """
        Performing multiplication
        """
        # Added mass
        for elem in self.Total_added_mass: elem.Matrix *= added_mass_dimension
        if self.compartment_flag:
            for elem in self.Hull_added_mass: elem.Matrix *= added_mass_dimension
        if self.Zero_freq_added_mass is not None:
            self.Zero_freq_added_mass *= added_mass_dimension
        if self.Inf_freq_added_mass is not None:
            self.Inf_freq_added_mass *= added_mass_dimension

        # Damping
        for elem in self.Total_damping: elem.Matrix *= damping_dimensions
        if self.compartment_flag:
            for elem in self.Total_pot_damping: elem.Matrix *= damping_dimensions
            for elem in self.Hull_damping: elem.Matrix *= damping_dimensions

        # Restoring
        self.Hydrostat_restoring *= restoring_dimensions
        if self.Total_restoring is not None:
            self.Total_restoring *= restoring_dimensions
        if self.mooring_flag:
            self.Mooring_restoring *= restoring_dimensions

        # Load RAO
        for l_RAO_h in self.LoadRAO:
            for elem in l_RAO_h: elem.amplitude *= excitation_dimensions

        # Displacement RAO
        for d_rao_h in self.MotionRAO:
            for elem in d_rao_h: elem.amplitude *= motions_dimensions

        if self.compartment_flag:
            for tank in self.tanks:
                for a in tank.added_mass: a.Matrix *= added_mass_dimension
                for b in tank.damping: b.Matrix *= damping_dimensions
                tank.restoring *= restoring_dimensions

        if len(self.mean_drift_horizontal) >0:
            for m_drift_h in self.mean_drift_horizontal:
                for elem in m_drift_h: elem.amplitude *= np.array([mean_drift_dimensions[0], mean_drift_dimensions[1],
                                                                   mean_drift_dimensions[5]])
        if len(self.mean_drift_full) > 0:
            for m_drift_h in self.mean_drift_full:
                for elem in m_drift_h: elem.amplitude *= mean_drift_dimensions

        if len(self.diff_freq_load_QTF) > 0:
            for m_drift_hi in self.diff_freq_load_QTF:
                for m_drift_hij in m_drift_hi:
                    for elem in m_drift_hij: elem.amplitude *= QTF_force_dimensions

        if len(self.diff_freq_motion_QTF) >0:
            for m_drift_hi in self.diff_freq_motion_QTF:
                for m_drift_hij in m_drift_hi:
                    for elem in m_drift_hij: elem.amplitude *= QTF_motion_dimensions


    def get_Total_added_mass(self, i, j):
        t_am = np.array([t.Matrix[i,j] for t in self.Total_added_mass])
        p = np.array([t.period for t in self.Total_added_mass])
        return t_am, p

    def get_Hull_added_mass(self, i, j):
        if self.Hull_added_mass is not None:
            h_am = np.array([t.Matrix[i,j] for t in self.Hull_added_mass])
            p = np.array([t.period for t in self.Hull_added_mass])
            return h_am, p
        else:
            warnings.warn("Hull added mass is None")
            return self.Hull_added_mass

    def get_Total_damping(self, i, j):
        t_damp = np.array([t.Matrix[i,j] for t in self.Total_damping])
        p = np.array([t.period for t in self.Total_damping])
        return t_damp, p

    def get_Total_pot_damping(self, i, j):
        if self.Total_pot_damping is not None:
            t_damp = np.array([t.Matrix[i,j] for t in self.Total_pot_damping])
            p = np.array([t.period for t in self.Total_pot_damping])
            return t_damp, p
        else:
            warnings.warn("Total potential damping is None")
            return self.Total_pot_damping

    def get_Hull_damping(self, i, j):
        if self.Hull_damping is not None:
            h_damp = np.array([t.Matrix[i,j] for t in self.Hull_damping])
            p = np.array([t.period for t in self.Hull_damping])
            return h_damp, p
        else:
            warnings.warn("Hull damping is None")
            return self.Hull_damping


    def get_Motion_RAO(self, i, heading):
        h_RAO = self.MotionRAO[heading]
        x_amp = np.array([x.amplitude[i] for x in h_RAO])
        x_ph = np.array([x.phase[i] for x in h_RAO])
        period = np.array([x.period for x in h_RAO])
        return x_amp, x_ph, period

    def get_Load_RAO(self, i, heading):
        h_RAO = self.LoadRAO[heading]
        f_amp = np.array([x.amplitude[i] for x in h_RAO])
        f_ph = np.array([x.phase[i] for x in h_RAO])
        period = np.array([x.period for x in h_RAO])
        return f_amp, f_ph, period

    def get_tanks_added_mass(self, i, j):
        if self.tanks is not None:
            periods = self.Environment_data.wave_periods
            temp_array = np.zeros_like(periods)
            for tank in self.tanks:
                temp_tank = np.array([t.Matrix[i,j] for t in tank.added_mass])
                temp_array += temp_tank
            return temp_array, periods
        else:
            warnings.warn("Tanks is None")
            return self.tanks

    def get_tanks_damping(self, i, j):
        if self.tanks is not None:
            periods = self.Environment_data.wave_periods
            temp_array = np.zeros_like(periods)
            for tank in self.tanks:
                temp_tank = np.array([t.Matrix[i,j] for t in tank.damping])
                temp_array += temp_tank
            return temp_array, periods
        else:
            warnings.warn("Tanks is None")
            return self.tanks

    def get_tanks_restoring(self):
        if self.tanks is not None:
            periods = self.Environment_data.wave_periods
            temp_array = np.zeros((6,6))
            for tank in self.tanks:
                temp_array += tank.restoring
            return temp_array, periods
        else:
            warnings.warn("Tanks is None")
            return self.tanks

    def get_QTF_load(self, i, heading1, heading2):
        if len(self.diff_freq_load_QTF):
            temp = self.diff_freq_load_QTF[heading1][heading2]
            period_i = np.unique(np.array([qtf.period_i for qtf in temp]))
            period_j = np.unique(np.array([qtf.period_j for qtf in temp]))
            qtf = np.full((len(period_i), len(period_j)), np.nan)
            for elem in temp:
                ii = np.argwhere(period_i==elem.period_i)[0][0]
                jj = np.argwhere(period_j ==elem.period_j)[0][0]
                qtf[ii][jj] = elem.amplitude[i]
            return qtf, period_i, period_j

    def get_Mean_drift_full(self, i, heading):
        if len(self.mean_drift_full):
            temp = self.mean_drift_full[heading]
            period = np.array([elem.period for elem in temp])
            amp = np.array([elem.amplitude[i] for elem in temp])
            return amp, period

    def get_Mean_drift_horizontal(self, i, heading):
        if len(self.mean_drift_horizontal):
            if i in [0, 1, 5]:
                if i==5: i=2
                temp = self.mean_drift_horizontal[heading]
                period = np.array([elem.period for elem in temp])
                amp = np.array([elem.amplitude[i] for elem in temp])
                return amp, period
            else:
                warnings.warn("Only horizontal dofs.")



Addedmassunits = np.array([[r'kg',r'kg',r'kg',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m'],
                  [r'kg',r'kg',r'kg',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m'],
                  [r'kg',r'kg',r'kg',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m'],
                  [r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot m^2$',r'kg$\cdot m^2$',r'kg$\cdot m^2$'],
                  [r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot m^2$',r'kg$\cdot m^2$',r'kg$\cdot m^2$'],
                  [r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot$m',r'kg$\cdot m^2$',r'kg$\cdot m^2$',r'kg$\cdot m^2$']])

Dampingunits = np.array([[r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$'],
                  [r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$'],
                  [r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$'],
                  [r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$'],
                  [r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$'],
                  [r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$',r'$\frac{kg\cdot m^2}{s}$']])

Excitation_units = np.array([r'N/m', r'N/m', r'N/m', r'N $\cdot$ m/m', r'N $\cdot$ m/m', r'N $\cdot$ m/m'])

Mean_drift_units = np.array([r'N/m$^2$', r'N/m$^2$', r'N/m$^2$', r'N $\cdot$ m/m$^2$', r'N $\cdot$ m/m$^2$', r'N $\cdot$ m/m$^2$'])

motionRAO_units = np.array([r'm/m', r'm/m', r'm/m', r'deg/m', r'deg/m', r'deg/m'])
















