# bb_cal
To run this code you will need fastmie.m from https://github.com/OceanOptics/MieTheory
The program 'Theoretical_slope.m' is used to compute the theoretical slope for bead calibrations of backscattering meters.
VSF and beam c are computed for beads distributed normally with center value 'D0', uncertainty in center, 'err_D0' and
standard deviation 'delta_D0' (taken straight from the calibration bead bottle - note, the manufacturer write that it is a k=2 uncertainty, hence we assume that err_D0 is TWICE the standard deviation of the mean value).
Values for VSF at wavelength with centroid wavelength 'wl' standard deviation 'delta_wl' as measured, for example, with a spectro-radiometer.
Program assumes that beam attenuation is measured with an AC instrument measuring the beam attenuation (acceptance angle 0.93degrees) at wavelength
'c_wl' each with uncertainty 'delta_c_wl' (wavelength does NOT have to match that of the backscattering sensor).
Backscatterin instrument angular reponse is given by a Gaussian with mean 'theta' and standard deviation 'd_theta'. Uncertainty in the mean values is given by 'err_theta'.
NB: Wavelength and bead-size should be in either in units of microns or nanometers and angular parameters are assumed in degrees.

Output: 'VSF' (mean VSF for 1bead/m^3), 'd_VSF' (uncertainty of VSF for 1bead/m^3), 'beam_c' (for 1bead/m^3), 'd_beam_c' (uncertainty in beam
attenaution for 1bead/m^3), 'ratio'-ratio of the two  and 'unc_ratio' (uncertainty in the ratio).

An example of how to call it:
[VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio]= Theoretical_slope([0.1 0.2 0.7],[0.004 0.006 0.007],[0.020 0.003 0.040],[0.470 0.555],[0.015 0.010],0.532,0.01,120,5,16);

for three beads, 0.1, 0.2 and 0.7um, uncertainty in the mean value 4, 6 and 7nm, dispersion around the meand value of 20, 3 and 40nm, for two backscatteriing wavelengths, 470 and 555nm with standard deviation of 15 and 10nm. For beam attenuation measured at 532nm with an ac-9 (assume stdev=10nm).
Backscattering angle assumed 120 +/-5 with a dispersion of 16degrees.

For any comments/questions please contact: emmanuel.boss@maine.edu
Thanks for Giorgio Dall'Olmo for suggesting many improvements!
