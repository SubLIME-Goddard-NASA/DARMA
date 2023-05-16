import os
import math
import csv
import datetime
import matplotlib.pyplot as plt
from numpy import array, trapz, std
#from sklearn.linear_model import LinearRegression
from scipy.integrate import simps
import numpy as np
import pandas as pd

class SpectraAnalysis:
    # Constructor to initialize variables for frequency and intensity of spectrum and power scan
    def __init__(self):
        self.frequency = []
        self.intensity = []
        self.freq_power = []
        self.int_power = []
        

    # Import data file
    # 'file_name' is data file name (.txt)
    # 'sens' is sensitivity setting in volts
    # 'delim' is the dilimeter which is defaulted to a single space, but can be changed to tab (\t), comma (,), etc.
    def import_data(self, file_name, sens, delim=' '):
        # Open file and read data file lines
        
        data = np.genfromtxt(file_name, skip_header = True, names = ['Frequency', 'Flux'])
        
        self.frequency = data['Frequency']
        
         # Convert intensity units to volts 
            
        self.intensity = (data['Flux']*sens)/10000
        
        # Prints True if frequency and intensity lists were made correctly
        if len(self.frequency) > 0 and len(self.intensity) > 0:
            return True
        else:
            return False
            
    # Import power scan file and convert intensity to volts
    # 'file_name' is data file name (.txt)
    # 'sens' is lock-in sensitivity setting in volts
    def import_power(self, file_name, sens, delim=' '):
        # Open file and read data file lines
        with open(file_name) as data_file:
            data_file_lines = data_file.readlines()
        
        self.freq_power = []
        self.int_power = []
        
        # Loop through all data points and save to lists
        for line in data_file_lines:
            values = line.split(delim)
            self.freq_power.append(float(values[0]))
            
            # Convert intensity units to volts
            volts_power = (float(values[len(values)-1])*sens)/10000
            
            # Append voltage values to intensity list
            self.int_power.append(volts_power)
        
        # Prints True if frequency and intensity lists were made correctly        
        if len(self.freq_power) > 0 and len(self.int_power) > 0:
            return True
        else:
            return False
        
    # Power correct data (i.e. divide spectrum by power scan)
    # 'freq_data' is frequency list of data
    # 'int_data' is intensity list of data
    # 'freq_power' is frequency list of power scan
    # 'int_power' is intensity list of power scan
    def power_correct(self, freq_data, int_data):
        int_correct = []
    
        for i in range(len(int_data)):
            for j in range(len(self.freq_power)):
                if freq_data[i] >= self.freq_power[j] and freq_data[i] < self.freq_power[j+1]:
                    if abs(freq_data[i] - self.freq_power[j]) < abs(freq_data[i] - self.freq_power[j+1]):
                        int_correct.append(int_data[i]/self.int_power[j])
                        
                    else:
                        int_correct.append(int_data[i]/self.int_power[j+1])
        
        return int_correct

    # Convert intensity units (arb units from lock-in amp) to volts (V)
    # 'int_value' is the intensity value
    # 'sens' is sensitivity setting in volts
    def int_volts(self, int_value, sens):
        return (int_value*sens)/10000

    # Convert intensity units to Power (Watts)
    # 'int_value' is the intensity value
    # 'sens' is sensitivity setting in volts
    # 'optical_responsivity' is the Optical responsivity of the detector. Look up in manual.
    def int_watts(self, int_value, sens, optical_responsivity):
        opt_resp = optical_responsivity * 1000   # Optical responsivity of detector (V)
        return (int_value*sens)*(1/10000)*(1/opt_resp)
        
       
    # Convert intensity units to Kelvin
    # 'int_value' is the intensity value
    # 'sens' is sensitivity setting in volts
    # 'opt_resp' is optical responsivity of detector. Look up in manual.
    # 'bandwidth' is the detector bandwidth. Look up in manual.
    def int_Kelvin(self, int_value, sens, opt_resp, bandwidth):
        kb = 1.3806*(10**-23)                    # Boltzmann constant (J/K)
        return  (int_value*sens)*(1/10000)*(1/opt_resp)*(1/(kb*bandwidth))
        
    # Calculate FWHM in frequency space (2f-line shape)
    # 'center' is center frequency of rotational transition from spectrum
    def get_peak_width(self, center, distance_from_center=30):
        
        # Find index of closest point in self.frequency to given center frequency
        index = [abs(y-center) for y in self.frequency].index(min(abs(x-center) for x in self.frequency))
        
        #distance_from_center = 30
        
        print(distance_from_center)
        
        left_list = [self.intensity[(index - distance_from_center) + u] for u in range(distance_from_center)]
        right_list = [self.intensity[index + v] for v in range(distance_from_center)]
        total_peak_list = left_list + right_list
        
        freq_list = [self.frequency[(index - len(left_list)) + w] for w in range(len(total_peak_list))]
        
        # Find indexs of local minimum intensity to the right and left of the given center
        left_min_index = total_peak_list.index(min(left_list))
        right_min_index = total_peak_list.index(min(right_list))
        
        # Get difference in frequency from local minima points. (i.e. frequency width from minimum of left lobe to minimum of right lobe)
        width = freq_list[right_min_index] - freq_list[left_min_index]
        
        # Calculate FWHM in frequency units
        w = width/(2*math.sqrt(3))
        FWHM_freq = 2*math.sqrt(2*math.log(2))*w
        
        print('Peak width (MHz): ' +str(width))
        print('w: ' +str(w))
        print('FWHM (MHz): ' +str(FWHM_freq))
        
        return FWHM_freq
    
    # Get peak width from file. This is used when you need peak width(s) for multiple peaks that are not in the same data file.
    # 'filename' is the data file name (.txt)
    # 'sens' is sensitivity setting in volts
    # 'center' is center frequency of the rotational transition in spectrum
    # 'delim=' is the deliminator; default is space ' '; can change to tab '\t' or comma ','
    def get_peak_width_from_file(self, filename, sens, center, delim=' '):
        if self.import_data(filename, sens, delim=delim):
            return self.get_peak_width(center)
        else:
            return 'Error: peak width calculation failed'
    
    # Calculate standard deviation from x1 to x2
    # When calculating standard deviation of noise: 'x1' is starting frequency value and 'x2' is ending frequency value.
    def get_std_dev(self, x1, x2):
        y_set = []
        
        for i, x_val in enumerate(self.frequency):
            if x_val >= x1 and x_val <= x2:
                y_set.append(self.intensity[i])
            
        std_dev = std(y_set)
        
        return std_dev
            
    # Calculate peak area from x1 to x2
    # 'center_freq' is center frequency value of peak being integrated
    # 'x' is the number of MHz from the center frequency to integrate across
    # 'boltz' is set as 1 if need peak area in intensity units of Kelvin (e.g., for Boltzmann analysis). Leave blank if you do not need to convert intensity to Kelvin.
    def get_peak_area(self, center_freq, x, boltz=None):
        # Initialize variables
        x_set = []
        y_set = []
        neg_y_set = []
        
        # Calculate integration range (x1 to x2) from given x
        x1 = center_freq - x
        x2 = center_freq + x
        
        # Set integration range from x1 to x2
        # x1 should start when the first negative lobe of 2f line starts to go below zero
        # x2 should be where the second negative lobe ends
        for i, x_val in enumerate(self.frequency):

            if x_val >= x1 and x_val <= x2:
                x_set.append(x_val)
                
                # y_set includes all positive values for 2f line
                if self.intensity[i] >= 0:
                    y_set.append(self.intensity[i])
                    neg_y_set.append(0)
                # neg_y_set includes absolute value of negative lobes of 2f line
                else:
                    y_set.append(0)
                    neg_y_set.append(abs(self.intensity[i]))
        
        # Power correct data
        y_set = self.power_correct(x_set, y_set)
        neg_y_set = self.power_correct(x_set, neg_y_set)        
        
        if boltz is not None:
            # Constants
            kb = 1.3806*(10**-23)          #Boltzmann (J/K)
            bandwidth = 2.4*(10**12)       #Detector bandwidth (Hz)
            opt_resp = 170000              #Optical responsivity of detector (V/W)
            
            y_set_Kelvin = []
            neg_y_set_Kelvin = []
            
            # Convert intensity to Kelvin and append to 'y_set_Kelvin'
            for int_val in y_set:
                int_val_Kelvin = (int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                y_set_Kelvin.append(int_val_Kelvin)
            
            # Convert intensity to Kelvin and append to 'neg_y_set_Kelvin'
            for neg_int_val in neg_y_set:
                neg_int_val_Kelvin = (neg_int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                neg_y_set_Kelvin.append(neg_int_val_Kelvin)
            
            
            # Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set_Kelvin)
            area_neg = self.integrate(x_set, neg_y_set_Kelvin)
            
            # Calculate total peak area = sum of positive lobe and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot
            
        else:
            #Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set)
            area_neg = self.integrate(x_set, neg_y_set)
            
            #Calculate total peak area = sum of positive lobes and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot

    # Calculate peak area from x1 to x2
    # 'center_freq' is center frequency value of peak being integrated
    # 'x' is the number of MHz from the center frequency to integrate across
    # 'boltz' is set as 1 if need peak area in intensity units of Kelvin (e.g., for Boltzmann analysis). Leave blank if you do not need to convert intensity to Kelvin.
    def get_peak_area_2(self, x1, x2, boltz=None):
        # Initialize variables
        x_set = []
        y_set = []
        neg_y_set = []
        
        # Calculate integration range (x1 to x2) from given x
        #x1 = center_freq - x
        #x2 = center_freq + x
        
        # Set integration range from x1 to x2
        # x1 should start when the first negative lobe of 2f line starts to go below zero
        # x2 should be where the second negative lobe ends
        for i, x_val in enumerate(self.frequency):
            if x_val >= x1 and x_val <= x2:
                x_set.append(x_val)
                
                # y_set includes all positive values for 2f line
                if self.intensity[i] >= 0:
                    y_set.append(self.intensity[i])
                    neg_y_set.append(0)
                # neg_y_set includes absolute value of negative lobes of 2f line
                else:
                    y_set.append(0)
                    neg_y_set.append(abs(self.intensity[i]))
        
        # Power correct data
        #y_set = self.power_correct(x_set, y_set)
        #neg_y_set = self.power_correct(x_set, neg_y_set)        
        
        if boltz is not None:
            # Constants
            kb = 1.3806*(10**-23)          #Boltzmann (J/K)
            bandwidth = 2.4*(10**12)       #Detector bandwidth (Hz)
            opt_resp = 170000              #Optical responsivity of detector (V/W)
            
            y_set_Kelvin = []
            neg_y_set_Kelvin = []
            
            # Convert intensity to Kelvin and append to 'y_set_Kelvin'
            for int_val in y_set:
                int_val_Kelvin = (int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                y_set_Kelvin.append(int_val_Kelvin)
            
            # Convert intensity to Kelvin and append to 'neg_y_set_Kelvin'
            for neg_int_val in neg_y_set:
                neg_int_val_Kelvin = (neg_int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                neg_y_set_Kelvin.append(neg_int_val_Kelvin)
            
            
            # Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set_Kelvin)
            area_neg = self.integrate(x_set, neg_y_set_Kelvin)
            
            # Calculate total peak area = sum of positive lobe and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot
            
        else:
            #Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set)
            area_neg = self.integrate(x_set, neg_y_set)
            
            #Calculate total peak area = sum of positive lobes and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot


    # Calculate peak area from x1 to x2
    # 'center_freq' is center frequency value of peak being integrated
    # 'x' is the number of MHz from the center frequency to integrate across
    # 'boltz' is set as 1 if need peak area in intensity units of Kelvin (e.g., for Boltzmann analysis). Leave blank if you do not need to convert intensity to Kelvin.
    def get_peak_area_3(self, center_freq, x, boltz=None):
        # Initialize variables
        x_set = []
        y_set = []
        neg_y_set = []
        
        # Calculate integration range (x1 to x2) from given x
        x1 = center_freq - x
        x2 = center_freq + x
        
        # Set integration range from x1 to x2
        # x1 should start when the first negative lobe of 2f line starts to go below zero
        # x2 should be where the second negative lobe ends
        for i, x_val in enumerate(self.frequency):
            if x_val >= x1 and x_val <= x2:
                x_set.append(x_val)
                
                # y_set includes all positive values for 2f line
                if self.intensity[i] >= 0:
                    y_set.append(self.intensity[i])
                    neg_y_set.append(0)
                # neg_y_set includes absolute value of negative lobes of 2f line
                else:
                    y_set.append(0)
                    neg_y_set.append(abs(self.intensity[i]))
        
        # Power correct data
        #y_set = self.power_correct(x_set, y_set)
        #neg_y_set = self.power_correct(x_set, neg_y_set)        
        
        if boltz is not None:
            # Constants
            kb = 1.3806*(10**-23)          #Boltzmann (J/K)
            bandwidth = 2.4*(10**12)       #Detector bandwidth (Hz)
            opt_resp = 170000              #Optical responsivity of detector (V/W)
            
            y_set_Kelvin = []
            neg_y_set_Kelvin = []
            
            # Convert intensity to Kelvin and append to 'y_set_Kelvin'
            for int_val in y_set:
                int_val_Kelvin = (int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                y_set_Kelvin.append(int_val_Kelvin)
            
            # Convert intensity to Kelvin and append to 'neg_y_set_Kelvin'
            for neg_int_val in neg_y_set:
                neg_int_val_Kelvin = (neg_int_val)*(1/opt_resp)*(1/(kb*bandwidth))
                neg_y_set_Kelvin.append(neg_int_val_Kelvin)
            
            
            # Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set_Kelvin)
            area_neg = self.integrate(x_set, neg_y_set_Kelvin)
            
            # Calculate total peak area = sum of positive lobe and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot
            
        else:
            #Calculate area of positive lobe and negative lobes
            area_pos = self.integrate(x_set, y_set)
            area_neg = self.integrate(x_set, neg_y_set)
            
            #Calculate total peak area = sum of positive lobes and negative lobes
            area_tot = area_pos + area_neg
            
            print(area_tot)
            return area_tot

        
    # Calculate peak area of peaks in different data files. This is used when you need peak area(s) for multiple peaks that are not in the same data file.
    #'filename' is data file name (.txt)
    # 'sens' is sensitivity setting in volts
    # 'center_freq' is center frequency of rotational transition
    # 'x' is the number of MHz from the center frequency to integrate across
    # 'delim=' is the deliminator; default is space ' '; can change to tab '\t' or comma ','
    # 'boltz' equals 1 if need peak area in intensity units of Kelvin (e.g., for Boltzmann analysis). Leave blank if you do not need to convert intensity to Kelvin.
    def get_peak_area_from_file(self, filename, sens, center_freq, x, delim=' ', boltz=None):
        if self.import_data(filename, sens, delim=delim):
            return self.get_peak_area(center_freq, x, boltz=boltz)
        else:
            return 'Error: peak area calculation failed'
            
            
    # Calculate peak area of peaks in different data files. This is used when you need peak area(s) for multiple peaks that are not in the same data file.
    #'filename' is data file name (.txt)
    # 'sens' is sensitivity setting in volts
    # 'center_freq' is center frequency of rotational transition
    # 'x' is the number of MHz from the center frequency to integrate across
    # 'boltz' equals 1 if need peak area in intensity units of Kelvin (e.g., for Boltzmann analysis). Leave blank if you do not need to convert intensity to Kelvin.
    def get_peak_area_from_file_2(self, filename, sens, center_freq, x, delim=' ', boltz=None):
        if self.import_data(filename, sens, delim=delim):
            return self.get_peak_area_3(center_freq, x, boltz=boltz)
        else:
            return 'Error: peak area calculation failed'
    
    # Integration function that uses the Composite Simpson's rule
    # 'freq' is frequency list
    # 'int' is intensity list
    def integrate(self, freq, int):
        # The y and x values
        y = array(int)
        x = array(freq)
    
        # Compute the area using the composite Simpson's rule.
        peak_area_s = simps(y,x=x)
    
        return peak_area_s

 
class BoltzmannAnalysis:
    # Initialize variables
    def __init__(self):
        self.__cat_file = ''
        self.__Q_300K = 0.0
        self.log_int = []
        self.E_l = []
        self.g_u = []

    # Set partition function value (Q) for catalog file(s). Temperature dependent value. Look up in spectral catalog.    
    def set_Q(self, q_val):
        self.__Q_300K = float(q_val)
    
    # Set catalog (.cat) file
    # 'file_name' is catalog file name (.cat)
    def set_cat_file(self, file_name):
        if os.path.isfile(file_name):
            self.__cat_file = file_name
            return True
        else:
            return False 

    # Looks up Log(Intensity) ('log_int'), lower state energy ('E_l'), and upper state degeneracy ('g_u') for each listed center frequency in catalog (.cat) file
    # 'centers' is list of center frequencies used in the analysis
    def __cat_lookup(self, centers):
        self.log_int = []
        self.E_l = []
        self.g_u = []
        
        c = centers[:]
        found_centers = []
        
        with open(self.__cat_file) as cat_f:
            lines = cat_f.readlines()
            for cent in c:
                for line in lines:
                    
                    freq_c = float(line[0:13].replace(' ', ''))
                    if cent == freq_c and not freq_c in found_centers:
                        print('found')
                        found_centers.append(freq_c)
                        self.log_int.append(float(line[21:29].replace(' ', '')))
                        self.E_l.append(float(line[31:41].replace(' ', '')))
                        self.g_u.append(float(line[41:44].replace(' ', '')))
        
                print(str(found_centers))
                print(cent)
                print(str(c))
        
        # Prints True if values were found and appended to lists
        if len(self.log_int) > 0 and len(self.E_l) > 0 and len(self.g_u) > 0:
            return True
        else:
            return False
                    
    # Looks up Log(Intensity) ('log_int'), lower state energy ('E_l'), and upper state degeneracy ('g_u') for each center frequency in .csv file
    # 'centers' is list of center frequencies used in the analysis
    def __cat_lookup_csv(self, centers):
        # Converts .cat file to .csv
        self.log_int = []
        self.E_l = []
        self.g_u = []
        
        with open(self.cat_file) as file:
            reader = csv.reader(file, delimiter=',')
            for row in reader:
                if float(row[0]) in centers:
                    self.log_int.append(float(vals[2]))
                    self.E_l.append(float(vals[4]))
                    self.g_u.append(float(vals[5]))
                    centers.remove(float(row[0]))

        # Prints True if values were found and appended to lists
        if len(self.log_int) > 0 and len(self.E_l) > 0 and len(self.g_u) > 0:
            return True
        else:
            return False

    # Exports data analysis results to .csv file
    def export_csv(self, data, save_path, lwa_folder, molecule):
        timestamp = datetime.datetime.now().strftime('%j%H%M%S')
        filename = save_path+'boltzmann'+lwa_folder+'_'+ molecule +'_'+timestamp+'.csv'
        
        with open (filename, "w") as cfile:
            writer = csv.writer(cfile, delimiter = ',', lineterminator = '\n')
            writer.writerow(['Eu (K)', 'y_set', 'Center (MHz)', 'Center (Hz)', 'Center (K)', 'Intensity', 'Peak Area (V)', 'El (K)', 'A Coefficient', 'B Coefficient', 'Q'])
            for i, x_item in enumerate(data['x_set']):
                writer.writerow([x_item, data['y_set'][i], data['v_mhz'][i], data['v_hz'][i], data['v_k'][i], data['int_set'][i], data['peak_area_set'][i], data['El_K'][i], data['A_coeff'][i], data['B_coeff'][i], data['Q'][i]])
        print('Created file: '+filename)

    # Boltzmann analysis procedure
    # 'centers' is list of center frequencies used in the analysis
    # 'peak_areas' is calculated peak areas used in analysis
    # 'peak widths' is calculated peak widths used in analysis
    # 'N_init' is initial guess for density. Leave as None if not considering optical depth correction.
    # 'T_init' is initial guess for temperature. Leave as None if not considering optical depth correction.
    def double_peak_check(self, center, freq_offset = 0.15):
        molecular_catalogue = pd.read_fwf(self.__cat_file, widths = [13,8,8,2,10,3,7,4,2,2,2,8,2,2], header = None)
        cat_freqs = molecular_catalogue[0]
        
        w = abs(cat_freqs - center) < freq_offset
        
        peak_freqs = cat_freqs[w]
        log_ints = molecular_catalogue[2][w]
        E_lower = molecular_catalogue[4][w]
        g_u = molecular_catalogue[5][w]
        
        return list(peak_freqs), list(log_ints), list(E_lower), list(g_u)
    

    def boltzmann_analysis(self, centers, peak_areas, peak_widths, N_init=None, T_init=None, freq_offset = 0.15):
        x_set = []
        y_set = []
        center_mhz = []
        center_hz = []
        center_k = []
        int_set = []
        peak_area_set = []
        E_l_kelvin_set = []
        A_coefficient_set = []
        B_coefficient_set = []
        Q_set = []
        
        # Prints values from catalog file
        if not self.__Q_300K == 0.0:
            if not self.__cat_file == '':
            
                if self.__cat_lookup(centers):
                    print('Upper state degeneracy: ' +str(self.g_u))
                    print('Log Intensity: ' +str(self.log_int))
                    print('Centers: ' + str(centers))
                    print('Peak widths: ' +str(peak_widths))
                    
                    # Initialize constants
                    h  = 6.6267*(10**-34)    #Planck (J s)
                    kb = 1.3806*(10**-23)    #Boltzmann (J/K)
                    c  = 299792458           #Speed of light (m/s)
                
                    for i, center_freq in enumerate(centers):
                        
                        if center_freq is None: 
                            continue
                        
                        peak_freqs, log_ints, E_lower, g_u = self.double_peak_check(float(center_freq), freq_offset)

                        if len(peak_freqs) == 2:
                            print(f'Double peak found at {peak_freqs[0]} and {peak_freqs[1]} MHz')
                            
                            # Conversions
                            center_freq = peak_freqs[0]
                            center_freq_Hz = center_freq*(10**6)  # Convert MHz to Hz
                            intensity_val = 10**(log_ints[0]) # Convert log intensity to intensity (nm^2 MHz)
                            E_l_kelvin1 = E_lower[0]/0.695028     # Convert lower state energy from cm-1 to K
                            E_trans1 = (center_freq_Hz*h)/kb       # Convert transition energy from Hz to K
                            E_u_kelvin1 = E_l_kelvin1 + E_trans1     # Upper state energy (K)

                            # Calculate Einstein A coefficient (spontaneous emission)
                            A1 = (2.7964*(10**-16))*(intensity_val)*(center_freq**2)*((10**self.__Q_300K)/g_u[0])*(1/((math.exp(-E_l_kelvin1/300))-(math.exp(-E_u_kelvin1/300))))

                            # Conversions
                            center_freq = peak_freqs[1]
                            center_freq_Hz = center_freq*(10**6)  # Convert MHz to Hz
                            intensity_val = 10**(log_ints[1]) # Convert log intensity to intensity (nm^2 MHz)
                            E_l_kelvin2 = E_lower[1]/0.695028     # Convert lower state energy from cm-1 to K
                            E_trans2 = (center_freq_Hz*h)/kb       # Convert transition energy from Hz to K
                            E_u_kelvin2 = E_l_kelvin2 + E_trans2     # Upper state energy (K)
                            
                            # Calculate Einstein A coefficient (spontaneous emission)
                            A2 = (2.7964*(10**-16))*(intensity_val)*(center_freq**2)*((10**self.__Q_300K)/g_u[1])*(1/((math.exp(-E_l_kelvin2/300))-(math.exp(-E_u_kelvin2/300))))
                            
                            E_u_kelvin = (E_u_kelvin1 + E_u_kelvin2)/2
                            E_l_kelvin = (E_l_kelvin1 + E_l_kelvin2)/2
                            E_trans = (E_trans1 + E_trans2)/2
                            
                            g_u = np.sum(g_u)
                            A_coefficient = A1 + A2

                        else:
                            # Conversions
                            g_u = g_u[0]
                            center_freq_Hz = center_freq*(10**6)  # Convert MHz to Hz
                            intensity_val = 10**(log_ints[0]) # Convert log intensity to intensity (nm^2 MHz)
                            E_l_kelvin = E_lower[0]/0.695028     # Convert lower state energy from cm-1 to K
                            E_trans = (center_freq_Hz*h)/kb       # Convert transition energy from Hz to K
                            E_u_kelvin = E_l_kelvin + E_trans     # Upper state energy (K)


                            # Calculate Einstein A coefficient (spontaneous emission)
                            A_coefficient = (2.7964*(10**-16))*(intensity_val)*(center_freq**2)*(10**self.__Q_300K/g_u)*(1/((math.exp(-E_l_kelvin/300))-(math.exp(-E_u_kelvin/300))))
       
                        # Convert Einstein A (spontaneous emission) to Einstein B (stimulated absorption)
                        B_coefficient = (A_coefficient*(c**3))/(math.pi*8*h*(center_freq_Hz**3))
                        
                        # Values appended to .csv file rows
                        center_hz.append(center_freq_Hz)
                        center_mhz.append(center_freq)
                        center_k.append(E_trans)
                        int_set.append(intensity_val)
                        peak_area_set.append(peak_areas[i])
                        E_l_kelvin_set.append(E_l_kelvin)
                        A_coefficient_set.append(A_coefficient)
                        B_coefficient_set.append(B_coefficient)
                        Q_set.append(self.__Q_300K)

                        # Calculate x-axis list ('x_set', i.e. upper state energy in Kelvin (E_u_Kelvin)) 
                        # and y-axis list ('y_set', i.e. gamma*peak area*optical depth correction) for Boltzmann/rotation diagram
                        x_set.append(E_u_kelvin)
                        
                        if N_init is not None and T_init is not None:                        
                            # Calculate FWHM in velocity units (cm/s)
                            FWHM_vel = (peak_widths[i]/center_freq)*c
                            
                            # Calculate upper state density
                            N_u = (N_init/10**self.__Q_300K)*g_u*math.exp(-E_u_kelvin/T_init)
                            
                            # Calculate optical depth
                            op_depth = (h/FWHM_vel)*N_u*B_coefficient*(math.exp(E_trans/T_init)-1)
                            
                            # Calculate optical depth correction factor
                            op_depth_corr = op_depth/(1-math.exp(-op_depth))

                            print('ODC: ' +str(op_depth_corr))
                            
                            y_set.append(math.log(((8*math.pi*kb*(center_freq_Hz**2)*peak_areas[i])/(h*(c**3)*A_coefficient*g_u))*(op_depth_corr)))

                        else:
                            y_set.append(math.log((8*math.pi*kb*(center_freq_Hz**2)*peak_areas[i])/(h*(c**3)*A_coefficient*g_u)))

                    # Optical depth correction: This function takes the initial guess values, plugs into the equations, 
                    # calculates temp and density from line of best fit, then uses those values for next guess values, and repeats until convergence.
                    if N_init is not None and T_init is not None:
                        x = array(x_set).reshape((-1, 1))
                        y = array(y_set)

                        model = LinearRegression().fit(x, y)
                        
                        y_int = model.intercept_
                        slope = model.coef_
                        print(y_int)
                        print(slope)

                        N_calc = math.exp(y_int)*self.__Q_300K
                        T_calc = -1/slope
                        
                        print('N_calc ='+str(N_calc))
                        print('T_calc ='+str(T_calc))
                        
                        if N_calc == N_init and T_calc == T_init:
                            return {'x_set':x_set, 'y_set':y_set, 'v_hz':center_hz, 'v_mhz':center_mhz, 'v_k':center_k, 'int_set':int_set, 'peak_area_set':peak_area_set, 'El_K':E_l_kelvin_set, 'A_coeff':A_coefficient_set, 'B_coeff':B_coefficient_set}

                        else:
                            return self.boltzmann_analysis(centers, peak_areas, peak_widths, N_init=N_calc, T_init=T_calc)
                    else:
                        return {'x_set':x_set, 'y_set':y_set, 'v_hz':center_hz, 'v_mhz':center_mhz, 'v_k':center_k, 'int_set':int_set, 'peak_area_set':peak_area_set, 'El_K':E_l_kelvin_set, 'A_coeff':A_coefficient_set, 'B_coeff':B_coefficient_set, 'Q' : Q_set}
                    
                else:
                    print('Error: __cat_lookup failed')
            else:
                print('Error: __cat_file not set')
        else:
            print('Error: Q_300K not set')
        
        return None
        


        
        #Plot Boltzmann diagram to determine rotational temperature and density
        #plt.plot(x, y)
        #plt.xlabel('Upper State Energy (K)')
        #plt.ylabel('Ln(I(T) (kb/$h^2 g_u v B_01$')