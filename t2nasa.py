###################################################################################### 
# Program to fit ab-initio thermodynamic data from JANEF to NASA-7 (CHEMKIN) format
# using the scipy least square fitting method 
# Author: Bastian Schnieder (bastian.schnieder@rub.de)
# Version: v1.0.0 (07.04.2025)
###################################################################################### 

try: 
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt
except ImportError:
    print("try: pip install [missing package]")

###################################################################################### 
# REQUIRED INPUT INFORMATION  
######################################################################################  
                          
DeltaH_f0   = -286.75E3                        # exp. value for heat of formation at T_zero (in cal)
S_R_zero    = 46.31                            # exp. value for entropy at T_zero (in cal)
T_zero      = 298.15                           # reference temperature for exp. values (in Kelvin)
H_RT_zero   = -DeltaH_f0 / 8.314472 / 298.15   # dimensionless enthalpy 
Rgas        = 1.9872                           # ideal gas constant in cal/K/mol

Low_temp_threshold  = 200                      # lower temperature range (in Kelvin)
High_temp_threshold = 4000                     # upper temperature range (in Kelvin)
Searching_T_step    = 10                       # temperature step (in Kelvin)

# NASA-7 output format details 

column_1_24  = 'CO2 L 7/88'                    # species name printed in output format 
column_25_44 = 'C 1O 2 0 0'                    # sumformula printed in output format
column_45    = 'G'                             # G "gas", S "solid", L "liquid"
f            = '%+10.8e\n'                     # decimal points for NASA output format  

######################################################################################    
# JANAF (Joint Army-Navy-Air Force) input format 
# More information: http://combustion.berkeley.edu/gri_mech/data/nasa_plnm.html
# Temperature T, Isobaric Heat Capacity Cp, Entropy S, Enthalpy H
# T(Klevin) | Cp(cal/(mol K)) | S (cal/(mol K)) | H (kcal/mol)
###################################################################################### 

# this part can be changed that the script reads the input data from an external file 

# insert your JANAF table here (delete exmaple values):
data = [
    [200.00, 10, 80, 65],
    ...
    ]
DATA = np.array(data)   

###################################################################################### 
# DO NOT CHANGE ANYTHING BELOW HERE 
###################################################################################### 

print("########################################################")
print("                 START OF THE PROGRAM                   ")
print("########################################################")

print("-------------- input information --------------")
print("lowest temperature:", Low_temp_threshold, "Kelvin")
print("highest temperature:", High_temp_threshold, "Kelvin")
print("temperature step:", Searching_T_step, "Kelvin")
print("-------------- output information --------------")


###################################################################################### 
# CONVERSION to dimensionless form and relative to its standard state 
###################################################################################### 

T_exp     = DATA[:, 0]
Cp_R_exp  = DATA[:, 1] / Rgas
CP_R_zero = interp1d(T_exp, Cp_R_exp, kind='linear')(T_zero)
S_R_exp   = DATA[:, 2] / Rgas
s_temp    = interp1d(T_exp, S_R_exp, kind='linear')(T_zero)
S_R_exp   = S_R_exp + (S_R_zero - s_temp)
h_temp    = interp1d(T_exp, DATA[:, 3], kind='linear')(T_zero)
H_R_exp   = (DATA[:, 3] - h_temp + H_RT_zero * Rgas * T_zero / 1000.) * 1000 / Rgas / T_exp

###################################################################################### 
# DEFINITION OF POLYNOMIALS 
###################################################################################### 

global a1, a2, a3, a4, a5, Cp_R_Cutting, Cp_R_slope_Cutting, Cutting_temperature

def dist_Cp(x, t, y_exp):
    a1, a2, a3, a4, a5 = x
    p = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4
    r = p - y_exp
    return r

def dist_enthalpie(x, t, y_exp):
    global a1, a2, a3, a4, a5
    a6 = x[0]
    p = a1 + a2 * t / 2 + a3 * t**2 / 3 + a4 * t**3 / 4 + a5 * t**4 / 5 + a6 / t
    r = p - y_exp
    return r

def dist_entropie(x, t, y_exp):
    global a1, a2, a3, a4, a5
    a7 = x[0]
    p = a1 * np.log(t) + a2 * t + a3 * t**2 / 2 + a4 * t**3 / 3 + a5 * t**4 / 4 + a7
    r = p - y_exp
    return r

def dist_Cp_HT(x, t, y_exp):
    global Cp_R_Cutting, Cp_R_slope_Cutting, Cutting_temperature
    a1, a2, a3, a4, a5 = x
    p = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4
    r = p - y_exp
    Cp_cutting = a1 + a2 * Cutting_temperature + a3 * Cutting_temperature**2 + a4 * Cutting_temperature**3 + a5 * Cutting_temperature**4
    Cp_cutting_plus = a1 + a2 * (Cutting_temperature + 1) + a3 * (Cutting_temperature + 1)**2 + a4 * (Cutting_temperature + 1)**3 + a5 * (Cutting_temperature + 1)**4
    Cp_R_slope_Cutting_local = Cp_cutting_plus - Cp_cutting
    return np.concatenate([r, [10 * np.abs(Cp_R_Cutting - p[0]), 10 * np.abs(Cp_R_slope_Cutting_local - Cp_R_slope_Cutting)]])

###################################################################################### 
# FINDING CUTTING TEMPERATURE
###################################################################################### 

res = []

for Cutting_temperature in range(Low_temp_threshold, High_temp_threshold + 1, Searching_T_step):
    cutting_pos = np.where(T_exp == Cutting_temperature)[0][0]

    # LOW TEMPERATURE LIMIT
    T_BT_exp = T_exp[:cutting_pos]
    H_R_BT_exp = H_R_exp[:cutting_pos]
    S_R_BT_exp = S_R_exp[:cutting_pos]
    Cp_R_BT_exp = Cp_R_exp[:cutting_pos]
    x_0 = [10, 1e-4, 1e-7, 1e-11, 1e-13]
    res0 = least_squares(dist_Cp, x_0, args=(T_BT_exp, Cp_R_BT_exp))
    a1, a2, a3, a4, a5 = res0.x
    res1 = least_squares(dist_enthalpie, [1], args=(T_BT_exp, H_R_BT_exp))
    a6_BT = res1.x[0]
    res2 = least_squares(dist_entropie, [1], args=(T_BT_exp, S_R_BT_exp))
    a7_BT = res2.x[0]

    # HIGH TEMPERATURE LIMIT
    T_HT_exp = T_exp[cutting_pos:]
    H_R_HT_exp = H_R_exp[cutting_pos:]
    S_R_HT_exp = S_R_exp[cutting_pos:]
    Cp_R_HT_exp = Cp_R_exp[cutting_pos:]
    T = Cutting_temperature
    Cp_R_Cutting = a1 + a2 * T + a3 * T**2 + a4 * T**3 + a5 * T**4
    T = T + 1
    Cp_R_Cutting_plus = a1 + a2 * T + a3 * T**2 + a4 * T**3 + a5 * T**4
    Cp_R_slope_Cutting = Cp_R_Cutting_plus - Cp_R_Cutting
    res3 = least_squares(dist_Cp_HT, x_0, args=(T_HT_exp, Cp_R_HT_exp))
    a1, a2, a3, a4, a5 = res3.x
    res4 = least_squares(dist_enthalpie, [1], args=(T_HT_exp, H_R_HT_exp))
    a6_HT = res4.x[0]
    res5 = least_squares(dist_entropie, [1], args=(T_HT_exp, S_R_HT_exp))
    a7_HT = res5.x[0]

    res.append(res0.cost + res1.cost + res2.cost + res3.cost + res4.cost + res5.cost)

I = np.argmin(res)
A = np.arange(Low_temp_threshold, High_temp_threshold + 1, Searching_T_step)
best_cutting_temperature = A[I]

# DEBUG (START)
# print("cutting temperature (K):",best_cutting_temperature)
# DEBUG (END)

#########################################################################################  
# MAIN ROUTINE START
print("this program uses a nonlinear least-square fit (scipy module)")
#########################################################################################  

######################################################################################### 
# Solve a nonlinear least-squares problem with bounds on the variables
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
######################################################################################### 

Low_temp_threshold = T_exp[0]
Cutting_temperature = best_cutting_temperature
High_temp_threshold = T_exp[-1]
cutting_pos = np.where(T_exp == Cutting_temperature)[0][0]

# LOW TEMPERATURE LIMIT ( T < cutting temperature)
T_BT_exp = T_exp[:cutting_pos]
H_R_BT_exp = H_R_exp[:cutting_pos]
S_R_BT_exp = S_R_exp[:cutting_pos]
Cp_R_BT_exp = Cp_R_exp[:cutting_pos]
x_0 = [10, 1e-4, 1e-7, 1e-11, 1e-13]
res0 = least_squares(dist_Cp, x_0, args=(T_BT_exp, Cp_R_BT_exp))
a1_BT, a2_BT, a3_BT, a4_BT, a5_BT = res0.x
a1, a2, a3, a4, a5 = res0.x
res1 = least_squares(dist_enthalpie, [1], args=(T_BT_exp, H_R_BT_exp))
a6_BT = res1.x[0]
res2 = least_squares(dist_entropie, [1], args=(T_BT_exp, S_R_BT_exp))
a7_BT = res2.x[0]

# HIGH TEMPERATURE LIMIT ( T > cutting temperature)
T_HT_exp = T_exp[cutting_pos:]
H_R_HT_exp = H_R_exp[cutting_pos:]
S_R_HT_exp = S_R_exp[cutting_pos:]
Cp_R_HT_exp = Cp_R_exp[cutting_pos:]
T = Cutting_temperature
Cp_R_Cutting = a1_BT + a2_BT * T + a3_BT * T**2 + a4_BT * T**3 + a5_BT * T**4
T = T + 1
Cp_R_Cutting_plus = a1_BT + a2_BT * T + a3_BT * T**2 + a4_BT * T**3 + a5_BT * T**4
Cp_R_slope_Cutting = Cp_R_Cutting_plus - Cp_R_Cutting
res3 = least_squares(dist_Cp_HT, x_0, args=(T_HT_exp, Cp_R_HT_exp))
a1_HT, a2_HT, a3_HT, a4_HT, a5_HT = res3.x
a1, a2, a3, a4, a5 = res3.x
res4 = least_squares(dist_enthalpie, [1], args=(T_HT_exp, H_R_HT_exp))
a6_HT = res4.x[0]
res5 = least_squares(dist_entropie, [1], args=(T_HT_exp, S_R_HT_exp))
a7_HT = res5.x[0]

###################################################################################### 
# MAIN ROUTINE END
###################################################################################### 

###################################################################################### 
# VISUALIZATION OF OPTIMIZED PARAMETERS 
###################################################################################### 

T_BT = np.arange(Low_temp_threshold, Cutting_temperature + 1)
H_R_BT = a1_BT + a2_BT * T_BT / 2 + a3_BT * T_BT**2 / 3 + a4_BT * T_BT**3 / 4 + a5_BT * T_BT**4 / 5 + a6_BT / T_BT
Cp_R_BT = a1_BT + a2_BT * T_BT + a3_BT * T_BT**2 + a4_BT * T_BT**3 + a5_BT * T_BT**4
S_R_BT = a1_BT * np.log(T_BT) + a2_BT * T_BT + a3_BT * T_BT**2 / 2 + a4_BT * T_BT**3 / 3 + a5_BT * T_BT**4 / 4 + a7_BT

T_HT = np.arange(Cutting_temperature, High_temp_threshold + 1)
H_R_HT = a1_HT + a2_HT * T_HT / 2 + a3_HT * T_HT**2 / 3 + a4_HT * T_HT**3 / 4 + a5_HT * T_HT**4 / 5 + a6_HT / T_HT
Cp_R_HT = a1_HT + a2_HT * T_HT + a3_HT * T_HT**2 + a4_HT * T_HT**3 + a5_HT * T_HT**4
S_R_HT = a1_HT * np.log(T_HT) + a2_HT * T_HT + a3_HT * T_HT**2 / 2 + a4_HT * T_HT**3 / 3 + a5_HT * T_HT**4 / 4 + a7_HT

plt.figure(figsize=(12, 10))
plt.subplot(2, 2, 1)
plt.plot(A, np.log10(res), 'o-', markersize=8, label='Residuals')
plt.plot(best_cutting_temperature, np.log10(res[I]), 'ro', markersize=10, label='Best Cutting Temperature')
plt.title('Residuals vs T cutting')
plt.ylabel('Log10 residuals')
plt.xlabel('Temperature (K)')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(T_BT, H_R_BT * T_BT * Rgas, 'r-', label='Fitted data (BT)')
plt.plot(T_HT, H_R_HT * T_HT * Rgas, 'b-', label='Fitted data (HT)')
plt.plot(T_exp, H_R_exp * T_exp * Rgas, '--g', label='Input data')
plt.title('Enthalpy vs T')
plt.ylabel('H (J/mol)')
plt.xlabel('Temperature (K)')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(T_BT, Cp_R_BT * Rgas, 'r-', label='Fitted data (BT)')
plt.plot(T_HT, Cp_R_HT * Rgas, 'b-', label='Fitted data (HT)')
plt.plot(T_exp, Cp_R_exp * Rgas, '--g', label='Input data')
plt.title('Heat capacity vs T')
plt.ylabel('Cp (J/(mol.K))')
plt.xlabel('Temperature (K)')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(T_BT, S_R_BT * Rgas, 'r-', label='Fitted data (BT)')
plt.plot(T_HT, S_R_HT * Rgas, 'b-', label='Fitted data (HT)')
plt.plot(T_exp, S_R_exp * Rgas, '--g', label='Input data')
plt.title('Entropy vs T')
plt.ylabel('S (J/(mol.K))')
plt.xlabel('Temperature (K)')
plt.legend()

plt.tight_layout()

###################################################################################### 
# SAVING NASA-7 OUTPUT FORMAT 
###################################################################################### 

with open('nasa.dat', 'w') as file:
    file.write(f'{column_1_24}{column_25_44}{column_45} {Low_temp_threshold:.2f} {High_temp_threshold:.2f} {Cutting_temperature:.2f} 1\n')
    file.write(f'{a1_HT:+10.8e} {a2_HT:+10.8e} {a3_HT:+10.8e} {a4_HT:+10.8e} {a5_HT:+10.8e} 2\n')
    file.write(f'{a6_HT:+10.8e} {a7_HT:+10.8e} {a1_BT:+10.8e} {a2_BT:+10.8e} {a3_BT:+10.8e} 3\n')
    file.write(f'{a4_BT:+10.8e} {a5_BT:+10.8e} {a6_BT:+10.8e} {a7_BT:+10.8e} 4\n')

print("optimized data is saved in NASA-7 format in nasa.dat")

visualize = input("Do you want to visualize the optimized parameters? (yes/no)")
if visualize.strip().lower() == "yes":
    plt.show()
    print("########################################################")
    print("                  END OF THE PROGRAM                    ")
    print("########################################################")
else: 
    print("########################################################")
    print("                  END OF THE PROGRAM                    ")
    print("########################################################")

# PRINT TO CONSOLE (START) 
#print(f'{column_1_24}{column_25_44}{column_45} {Low_temp_threshold:.2f} {High_temp_threshold:.2f} {Cutting_temperature:.2f} 1')
#print(f'{a1_HT:+10.8e} {a2_HT:+10.8e} {a3_HT:+10.8e} {a4_HT:+10.8e} {a5_HT:+10.8e} 2')
#print(f'{a6_HT:+10.8e} {a7_HT:+10.8e} {a1_BT:+10.8e} {a2_BT:+10.8e} {a3_BT:+10.8e} 3')
#print(f'{a4_BT:+10.8e} {a5_BT:+10.8e} {a6_BT:+10.8e} {a7_BT:+10.8e} 4')
# PRINT TO CONSOLE (END)

######################################################################################################  
# This code was inspired by the Matlab Converter from Raphael Boichot.
# The Matlab code is available at: 
# https://github.com/Raphael-Boichot/JANAF-thermochemical-tables-to-NASA-Glenn-coefficients-converter
######################################################################################################  
