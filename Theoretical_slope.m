function [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio, ratio_nom, ang, VSF_nom, beam_c_nom] = Theoretical_slope(D0,err_D0,delta_D0,wl,delta_wl,c_wl,delta_c_wl,theta,err_theta,d_theta);
% To run this code you will need fastmie.m from https://github.com/OceanOptics/MieTheory The program 'Theoretical_slope.m' 
% is used to compute the theoretical slope for bead calibrations of backscattering meters. 
% 
% The slope is the ratio of volume scattering function (VSF, [m^-1 Sr^-1]) to beam atenuation [m^-1] by a given single-angle
% backscattering sensor such as the eco series from Seabird.
%
%Together with lab measurements of beam attenuation at a specific
%wavelength and VSF at a specific (possibly other) wavelength it provides
%means to compute the calibration coefficient of a sensor
%
% VSF and beam c are computed for beads distributed with a normal distribution with center value 'D0', uncertainty in center, 'err_D0' and standard deviation 'delta_D0' 
% (taken straight from the polystyrene calibration bead bottle). 
% 
% Note that the manufacturer specify that err_D0 is a k=2 estimate, hence
% we treat it as 2*standard deviation of the mean.
% 
% Values for VSF at wavelength with centroid wavelength 'wl' in vacuo (we
% correct for value in water in the code).
% deviation 'delta_wl' as measured, for example, with a spectro-radiometer pointed at the sensor. These also are corrected. 
% 
% Program assumes that beam attenuation is measured with a Seabird Scientific AC instrument measuring the beam attenuation 
% (acceptance angle 0.93degrees) at wavelength 'c_wl' each with uncertainty 'delta_c_wl' (wavelength does NOT have to match that of the backscattering sensor). 
% Backscattering instrument angular reponse is given by a Gaussian with mean 'theta' and standard deviation 'd_theta'. Uncertainty in 
% the mean values is given by 'err_theta'. 
%
% NB: UNITS: Wavelength and bead-size should be in either in units of microns or nanometers and angular parameters are assumed in degrees.
% Output: 'VSF' (mean VSF for 1 bead m^-3 of water), 'd_VSF' (uncertainty of VSF for 1 bead m^-3), 'beam_c' (for 1bead m^-3), 'd_beam_c' 
% (uncertainty in beam attenaution for 1bead m^-3), 'ratio'-ratio of VSF/beam-c and 'unc_ratio' (uncertainty in the ratio).
%
%Bead index of refraction is taken from Jones et al., 2013
%
% An example of how to call it: [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio, ratio_nom, ang, VSF_nom, beam_c_nom]= Theoretical_slope([0.1 0.2 0.7],[0.004 0.006 0.007],[0.020 0.003 0.040],[0.470 0.555],[0.015 0.010],0.532,0.01,120,5,16);
% 
% For three beads of size 0.1, 0.2 and 0.7um, uncertainty in the mean value 4, 6 and 7nm, dispersion around the meand value of 20, 3 and 40nm, 
% for two backscatteriing wavelengths, 470 and 555nm with standard deviation of 15 and 10nm. For beam attenuation measured at 532nm with 
% an ac-9 (assume stdev=10nm). Backscattering angle assumed 120 +/-5 with a dispersion of 16degrees.
% 
% For any comments/questions please contact: emmanuel.boss@maine.edu
% Emmanuel Boss 2023-04-26 + help from Giorgio Dall'Olmo

% check is Octave is being run, and if so, load needed packages
	if exist('OCTAVE_VERSION',  'builtin')
		pkg load statistics
		pkg load nan
	end

	%%% used in all Mie computations
	accept_angle = 0.93; % acceptance angle of transmissometer
	nang = 100*3 + 1; %number of angles for mie computation
	dang = pi/2/(nang - 1); %angular resolution in radians
	ang = [0: dang: pi]*180/pi; %angles in degrees
	T = 20; %temperature of calibration water, for index of refraction computations

	NN = length(D0); %number of different size beads
	KK = length(wl); %number of different wavelengths of beta 
	JJ = length(c_wl); %number of wavelengths of beam_c
	NNN = 200; %number of size parameter discritization wanted for a given size bead and a given wavelength
	N = 5000; % number of random realizations
	now_1 = now; % start timer

	% initialise arrays used to collect output of for loops (this is key to speed up for loops in octave)
	rho = nan(NNN + 1, 1);
	beta__ = nan(NNN + 1, length(ang));
	VSF = nan(NN, KK, length(ang));
	VSF_nom = nan(NN, KK, length(ang));
	d_VSF = nan(NN, KK, length(ang));
	Q_ext = nan(NNN + 1, 1);
	r = nan(N, 1);
	c = nan(N, 1);
	
	beam_c = nan(NN, JJ);
	d_beam_c = nan(NN, JJ);
	beam_c_nom = nan(NN, JJ);

    beta_ = nan(N, 1);
	mean_beta_ = nan(KK, JJ);
    d_beta_ = nan(KK, JJ);
    ratio = nan(NN, KK, JJ);
    unc_ratio = nan(NN, KK, JJ);

	for nn = 1 : NN %number of bead sizes
	    for k = 1 : KK %number of wl for beta
	        %do the necessary Mie computations.

	        max_wl = (wl(k) + 3*delta_wl(k)); %99th percentile, in water.
	        min_wl = (wl(k) - 3*delta_wl(k)); %99th percentile, in water.
	        max_D = D0(nn) + err_D0(nn) + 3*delta_D0(nn); %99th percentile.
	        min_D = D0(nn) - err_D0(nn) - 3*delta_D0(nn); %99th percentile.

	        %parameter space for which we will do Mie computations.
	        min_D_over_lambda = min_D/max_wl;
	        max_D_over_lambda = max_D/min_wl;
	        dd = (max_D_over_lambda - min_D_over_lambda)/NNN;

	        for i = 1 : NNN + 1
	            [n(i), nm] = IoR(wl(k), T);   %compute index of refraction of beads relative to water
	            rho(i) = pi*nm*(min_D_over_lambda + (i - 1)*dd); %rho
	            [S1 S2 Qb Qc Qback] = fastmie(rho(i), n(i), nang);
	            S11 = 0.5*((abs(S1)).^2  +  (abs(S2)).^2);
	            S11_rad = S11.*sin(ang'/180*pi)*2*pi;
	            S11_int = integrate(S11_rad, dang, length(S11_rad));
	            beta__(i, :) = Qb*S11/S11_int; %table of VSF for each wl
            end
            %compute nominal value
            [n_nom, nm] = IoR(wl(k), T);   %compute index of refraction of beads relative to water
	        rho_nom = pi*nm*(D0(nn)/wl(k)); 
            [S1 S2 Qb Qc Qback] = fastmie(rho_nom, n_nom, nang);
            S11 = 0.5*((abs(S1)).^2  +  (abs(S2)).^2);
            S11_rad = S11.*sin(ang'/180*pi)*2*pi;
            S11_int = integrate(S11_rad, dang, length(S11_rad));
	        beta_nom(:) = Qb*S11/S11_int; %nominal VSF for each wl
            
			
	        %sample randomly from wavelength and bead size space based on bead and wavelength distribution.
	        for j = 1 : N %number of random realization
	            DD = normrnd(normrnd(D0(nn), err_D0(nn)/2), delta_D0(nn)); %choose the random bead size
				if any(DD<=0)
					keyboard()
				end	
	            lambda = normrnd(wl(k), delta_wl(k)); %choose the random wavelength
				if any(lambda<=0)
					keyboard()
                end	
                [nnn, nm] = IoR(lambda, T); 
	            rr = pi*nm*DD/lambda;
 
	            beta(j, :) = pi*DD^2/4*interp1(rho, beta__, rr, 'linear','extrap'); % allow for extrapolation 
	        end
			
	        VSF(nn, k, :) = nanmean(beta);
	        d_VSF(nn, k, :) = nanstd(beta);
            VSF_nom(nn,k,:)=pi*D0(nn).^2/4*beta_nom;
	    end
	end
	 
	for nn = 1: NN %number of bead sizes
	    for jj = 1: JJ %number of wavelengths of beam-c
	        %do the necessary Mie computations.
	        max_wl = c_wl(jj) + 3*delta_c_wl(jj); %99th percentile.
	        min_wl = c_wl(jj) - 3*delta_c_wl(jj); %99th percentile.
	        max_D = D0(nn) + err_D0(nn) + 3*delta_D0(nn); %99th percentile.
	        min_D = D0(nn) - err_D0(nn) - 3*delta_D0(nn); %99th percentile.
            min_D_over_lambda = min_D/max_wl;
	        max_D_over_lambda = max_D/min_wl;
	        dd = (max_D_over_lambda - min_D_over_lambda)/NNN;

	        for i = 1: NNN + 1
	            [n(i),  nm] = IoR(c_wl(jj), T);
	            rho(i) = pi*nm*(min_D_over_lambda + (i - 1)*dd); %rho
	            [S1 S2 Qb Qc Qback] = fastmie(rho(i), n(i), nang);
	            S11 = 0.5*((abs(S1)).^2  +  (abs(S2)).^2);
	            S11_rad = S11.*sin(ang'/180*pi)*2*pi;
	            S11_int = integrate(S11_rad, dang, length(S11_rad));
	            %computations for beam attenuation: correcting for acceptance angle of ac - 9
	            delta = (pi - accept_angle*pi/180)/length(ang);
	            ang_ = [accept_angle*pi/180: delta: pi]*180/pi;
	            nan_ = length(ang_);
	            S11_ = interp1(ang, S11, ang_, 'pchip')';
	            S11_rad_ = S11_.*sin(ang_'*pi/180)*2*pi; %taking care of azymuthal weighing
	            S11_int_ = integrate(S11_rad_, delta, nan_);
	            Qb_corr = Qb*S11_int_/S11_int; % correct Qb for acceptance angle.
	            Qa = Qc - Qb;
	            Q_ext(i) = Qa + Qb_corr;
	        end

            %compute nominal value
            [n_nom,  nm] = IoR(c_wl(jj), T);  %compute index of refraction of beads relative to water
	        rho_nom = pi*real(n_nom)*(D0(nn)/c_wl(jj)); 
            [S1 S2 Qb Qc Qback] = fastmie(rho_nom, n_nom, nang);
		    beam_c_nom(nn, jj)=Qc*pi*D0(nn)^2/4;
            

	        %sample randomly from wavelength and bead size space.
	        for j = 1: N
	            DD = normrnd(normrnd(D0(nn), err_D0(nn)), delta_D0(nn)); 
				if any(DD<=0)
					keyboard()
				end	
	            lambda = normrnd(c_wl(jj), delta_c_wl(jj));
				if any(lambda<=0)
					keyboard()
                end	
                [nnn, nm] = IoR(lambda, T); 
	            rr = pi*real(nnn)*DD/lambda;
	            c(j) = pi*DD^2/4*interp1(rho, Q_ext, rr, 'linear','extrap');
            end
	        beam_c(nn, jj) = mean(c);
	        d_beam_c(nn, jj) = std(c);
	    end
	end


	for nn = 1: NN
	    for k = 1: KK
	        for jj = 1: JJ
	          for j = 1: N %this is where the angular distribution of the scattering angle is taken into account
	            theta_ = normrnd(normrnd(theta, err_theta), d_theta);
	            if theta_<0
	                theta_ = theta_ + 180;
	            end
	            if theta_>180
	                theta_ = theta_ - 180;
	            end
				
				% create random VSF
				rnd_VSF = normrnd(squeeze(VSF(nn, k, :)), squeeze(d_VSF(nn, k, :)));
				
	            beta_(j) = interp1(ang, rnd_VSF, theta_, 'linear');
	          end
			  beta_nom=interp1(ang, squeeze(VSF_nom(nn,k,:)), theta, 'linear');
	          mean_beta_(k, jj) = nanmean(beta_);
	          d_beta_(k, jj) = nanstd(beta_);
	          ratio(nn, k, jj) = mean_beta_(k, jj)./beam_c(nn, jj);
	          unc_ratio(nn, k, jj) = sqrt((d_beta_(k, jj)/mean_beta_(k, jj)).^2 + (d_beam_c(nn, jj)/beam_c(nn, jj))^2)*ratio(nn, k, jj);
              ratio_nom(nn, k, jj)=beta_nom./beam_c(nn, jj);
	        end
	    end
	end

	now_2 = now;
	time_lapse = (now_2 - now_1)*24*3600 %in seconds
	
	return
end



function [n, nm] = IoR(wl, T)
%function to compute the index of refraction of beads relative to pure water of
%temperature T.
	np = 1.5718  +  0.008412/(wl^2)  +  0.000235/(wl^4); %Jones et al.
	ni = 0.0003;  %imaginary part of index of refraction
	n0 = 1.31405; n1 = 1.779e-4; n2 = -1.05e-6; n3 = 1.6e-8; n4 = -2.02e-6; n5 = 15.868; n6 = 0.01155; n7 = -0.00423; n8 = -4382; n9 = 1.1455e6; %index of refraction of water Quan and Fry,  1995
	nm = n0 + n4*T^2 + (n5 + n7*T)/(wl*1000) + n8/(wl*1000)^2 + n9/(wl*1000)^3; %checked
	n = (np + ni*sqrt(-1))/nm; %index of refraction relative to water
	return
end



function intgr = integrate(y, dx, N)
%function to perform an integral of y with respect to x,  with dx the
%increment and N the total number of y's to integrate.
	intgr = 3/8*(y(1) + y(N)) + 7/6*(y(2) + y(N - 1)) + 23/24*(y(3) + y(N - 2));
	intgr = (intgr + sum(y(4: N - 3)))*dx;
end
