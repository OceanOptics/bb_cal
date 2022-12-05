function [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio] = Theoretical_slope(D0,err_D0,delta_D0,wl,delta_wl,c_wl,delta_c_wl,theta,err_theta,d_theta);
% Program to compute the theoretical slope for bead calibrations.
% VSF and beam c are computed for beads distributed normally with center value D0, uncertainty in center, err_D0 and
% standard deviation delta_D0 (taken straight from the calibration bead bottle.
% Values for VSF at wavelenght with centroid wavelength wl standard
% deviation delta_wl as measured with a radiometer.
% Program assumes that beam attenuation is measured with an AC instrument measuring the beam attenuation (acceptance angle 0.93degrees) at wavelength
% c_wl each with uncertainty delta_c_wl
% BB instrument angular reponse is given by a Gaussian with mean theta and
% standard deviation d_theta. Uncertainty in the mean values is given by err_theta.
% NB: everything should be in either micron or nanometer (wavelength, bead
% size) and theta parameters re assumed in degrees.
% Output: VSF (mean VSF for 1bead/m^3), D_VSF (uncertainty of VSF for
% 1bead/m^3), beam_c (for 1bead/m^3), d_beam_c (uncertainty in beam
% attenaution for 1bead/m^3), ratio-ratio of the two  and unc_ratio
% (uncertainty in the ratio.
% An example of how to call it:
% [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio]= Theoretical_slope([0.1 0.2 0.7],[0.004 0.006 0.007],[0.020 0.003 0.040],[0.470 0.555],[0.015 0.010],0.532,0.01,120,5,16);
% for three beads, 0.1, 0.2 and 0.7um, uncertainty in the mean value 4, 6
% and 7nm, dispersion around the meand value of 20, 3 and 40nm, for two backscatteriing wavelengths, 470 and 555nm with standard deviation of 15 and 10nm. 
% For beam attenuation measured at 532nm with an ac-9 (assume stdev=10nm).
% Backscattering angle assume 120 +/-5 with a dispersion of 16degrees.
%Emmanuel Boss 11-30-2022

%%% used in all Mie computations
accept_angle=0.93; % acceptance angle of transmissometer
nang=100*3+1; %number of angles for mie computation
dang=pi/2/(nang-1); %angular resolution in radians
nan=2*nang-1; %number of angles for computation
ang=[0:dang:pi]*180/pi; %angles in degrees
T=24; %temperature of calibration water, for index of refraction computations

NN=length(D0); %number of beads
KK=length(wl); %number of wavelengths of beta 
JJ=length(c_wl); %number of wavelengths of beam_c
now_1=now;

for nn=1:NN
    for k=1:KK
        %do the necessary Mie computations.
        max_wl=wl(k)+3*delta_wl(k); %99th percentile.
        min_wl=wl(k)-3*delta_wl(k); %99th percentile.
        max_D=D0(nn)+err_D0(nn)+3*delta_D0(nn); %99th percentile.
        min_D=D0(nn)-err_D0(nn)-3*delta_D0(nn); %99th percentile.

        %parameter space for which we will do Mie computations.
        min_D_over_lambda=min_D/max_wl;
        max_D_over_lambda=max_D/min_wl;
        NNN=500; %number of size parameter discritization wanted
        dd=(max_D_over_lambda-min_D_over_lambda)/NNN;

        for i=1:NNN+1
            [n, nm]=IoR(wl(k),T);
            rho(i)=pi*real(n)*(min_D_over_lambda+(i-1)*dd); %rho
            [S1 S2 Qb Qc Qback]=fastmie(rho(i),n,nang);
            S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
            S11_rad=S11.*sin(ang'/180*pi)*2*pi;
            S11_int=integrate(S11_rad,dang,length(S11_rad));
            beta__(i,:)=Qb*S11/S11_int; %table VSF
        end
        %sample randomly from wavelength and bead size space based on bead and wavelength distribution.
        for j=1:20000
            DD=normrnd(D0(nn),err_D0(nn))+normrnd(0,delta_D0(nn));
            lambda=normrnd(wl(k),delta_wl(k));
            rr(j)=pi*nm*DD/lambda;
            beta(j,:)=pi*DD^2/4*interp1(rho,beta__,rr(j),'linear','extrap');
        end
        VSF(nn,k,:)=mean(beta);
        d_VSF(nn,k,:)=std(beta);
    end
end
clear rr
 %number of wavelengths of beta
for nn=1:NN
    for jj=1:JJ
        %do the necessary Mie computations.
        max_wl=c_wl(jj)+3*delta_c_wl(jj); %99th percentile.
        min_wl=c_wl(jj)-3*delta_c_wl(jj); %99th percentile.
        max_D=D0(nn)+err_D0(nn)+3*delta_D0(nn); %99th percentile.
        min_D=D0(nn)-err_D0(nn)-3*delta_D0(nn); %99th percentile.
        dd=(max_D_over_lambda-min_D_over_lambda)/900;

        for i=1:NNN+1
            [n, nm]=IoR(c_wl(jj),T);
            rho(i)=pi*real(n)*(min_D_over_lambda+(i-1)*dd); %rho
            [S1 S2 Qb Qc Qback]=fastmie(rho(i),n,nang);
            S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
            S11_rad=S11.*sin(ang'/180*pi)*2*pi;
            S11_int=integrate(S11_rad,dang,length(S11_rad));
            %computations for beam attenuation:correcting for acceptance angle of ac-9
            delta=(pi-accept_angle*pi/180)/length(ang);
            ang_=[accept_angle*pi/180:delta:pi]*180/pi;
            nan_=length(ang_);
            S11_=interp1(ang,S11,ang_,'pchip')';
            S11_rad_=S11_.*sin(ang_'*pi/180)*2*pi; %taking care of azymuthal weighing
            S11_int_=integrate(S11_rad_,delta,nan_);
            Qb_corr=Qb*S11_int_/S11_int; % correct Qb for acceptance angle.
            Qa=Qc-Qb;
            Q_ext(i)=Qa+Qb_corr;
        end
        %sample randomly from wavelength and bead size space.
        for j=1:20000
            DD=normrnd(D0(nn),err_D0(nn))+normrnd(0,delta_D0(nn));
            lambda=normrnd(c_wl(jj),delta_c_wl(jj));
            r(j)=pi*nm*DD/lambda;
            c(j,:)=pi*DD^2/4*interp1(rho,Q_ext,r(j),'linear','extrap');
        end
         beam_c(nn,jj)=mean(c);
         d_beam_c(nn,jj)=std(c);
    end
end
for nn=1:NN
    for k=1:KK
        for jj=1:JJ
          for j=1:20000 %this is where the angular distribution of theta is taken into account
            theta_=normrnd(theta,err_theta)+normrnd(0,d_theta);
            if theta_>180
                theta_=theta_-180;
            end
            beta_(j)=interp1(ang,normrnd(squeeze(VSF(nn,k,:)),squeeze(d_VSF(nn,k,:))),theta_,'linear');
          end
          mean_beta_(k,jj)=nanmean(beta_);
          d_beta_(k,jj)=nanstd(beta_);
          ratio(nn,k,jj)=mean_beta_(k,jj)./beam_c(nn,jj)
          unc_ratio(nn,k,jj)=sqrt((d_beta_(k,jj)/mean_beta_(k,jj)).^2+(d_beam_c(nn,jj)/beam_c(nn,jj))^2)*ratio(nn,k,jj);
        end
    end
end

now_2=now;
time_lapse=(now_2-now_1)*24*3600 %in seconds
return
end

function [n, nm]=IoR(wl,T)
%function to compute the index of refraction of beads relative to pure water of
%temperature T.
np=1.5718 + 0.008412/(wl^2) + 0.000235/(wl^4); %Jones et al.
ni=0.0003;  %imaginary part of index of refraction
n0=1.31405; n1=1.779e-4; n2=-1.05e-6; n3=1.6e-8; n4=-2.02e-6; n5=15.868; n6=0.01155; n7=-0.00423; n8=-4382; n9=1.1455e6; %index of refraction of water Quan and Fry, 1995
nm=n0+n4*T^2+(n5+n7*T)/(wl*1000)+n8/(wl*1000)^2+n9/(wl*1000)^3; %checked
n=(np+ni*sqrt(-1))/nm; %index of refraction relative to water
return
end

function intgr=integrate(y,dx,N)
%function to perform an integral of y with respect to x, with dx the
%increment and N the total number of y's to integrate.
intgr=3/8*(y(1)+y(N))+7/6*(y(2)+y(N-1))+23/24*(y(3)+y(N-2));
intgr=(intgr+sum(y(4:N-3)))*dx;
end