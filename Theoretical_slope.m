function [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio] = Theoretical_slope(D0,err_D0,delta_D0,wl,delta_wl,c_wl,delta_c_wl,theta,err_theta,d_theta);
% To run this code you will need fastmie.m from https://github.com/OceanOptics/MieTheory The program 'Theoretical_slope.m' 
% is used to compute the theoretical slope for bead calibrations of backscattering meters. VSF and beam c are computed for 
% beads distributed normally with center value 'D0', uncertainty in center, 'err_D0' and standard deviation 'delta_D0' 
% (taken straight from the calibration bead bottle. 
% 
% Note that the manufacturer specify that err_D0 is a k=2 estimate, hence
% we treat it as 2*standard deviation of the mean.
% 
% Values for VSF at wavelength with centroid wavelength 'wl' standard 
% deviation 'delta_wl' as measured, for example, with a spectro-radiometer. Program assumes that beam attenuation is 
% measured with an AC instrument measuring the beam attenuation (acceptance angle 0.93degrees) at wavelength 'c_wl' 
% each with uncertainty 'delta_c_wl' (wavelength does NOT have to match that of the backscattering sensor). Backscatterin 
% instrument angular reponse is given by a Gaussian with mean 'theta' and standard deviation 'd_theta'. Uncertainty in 
% the mean values is given by 'err_theta'. NB: Wavelength and bead-size should be in either in units of microns or nanometers 
% and angular parameters are assumed in degrees.
% Output: 'VSF' (mean VSF for 1bead/m^3), 'd_VSF' (uncertainty of VSF for 1bead/m^3), 'beam_c' (for 1bead/m^3), 'd_beam_c' 
% (uncertainty in beam attenaution for 1bead/m^3), 'ratio'-ratio of the two and 'unc_ratio' (uncertainty in the ratio).
% 
% An example of how to call it: [VSF,d_VSF, beam_c, d_beam_c,ratio,unc_ratio]= Theoretical_slope([0.1 0.2 0.7],[0.004 0.006 0.007],[0.020 0.003 0.040],[0.470 0.555],[0.015 0.010],0.532,0.01,120,5,16);
% 
% for three beads, 0.1, 0.2 and 0.7um, uncertainty in the mean value 4, 6 and 7nm, dispersion around the meand value of 20, 3 and 40nm, 
% for two backscatteriing wavelengths, 470 and 555nm with standard deviation of 15 and 10nm. For beam attenuation measured at 532nm with 
% an ac-9 (assume stdev=10nm). Backscattering angle assumed 120 +/-5 with a dispersion of 16degrees.
% 
% For any comments/questions please contact: emmanuel.boss@maine.edu
%  Emmanuel Boss 12-07-2022 with help from Giorgio Dall'Olmo

% check is Octave is being run, and if so, load needed packages
if exist('OCTAVE_VERSION',  'builtin')
    pkg load statistics
    pkg load nan
end

%%% used in all Mie computations
accept_angle=0.93; % acceptance angle of transmissometer
nang=100*3+1; %number of angles for mie computation
dang=pi/2/(nang-1); %angular resolution in radians
ang=[0:dang:pi]*180/pi; %angles in degrees
T=24; %temperature of calibration water, for index of refraction computations
NNNN=100000; %number of sampling tried for each MC computation
NNN=500; %number of size parameter discritization wanted
NN=length(D0); %number of different beads
KK=length(wl); %number of different wavelengths of beta 
JJ=length(c_wl); %number of different wavelengths of beam_c
now_1=now;

% initialise arrays used to collect output of for loops (this is key to speed up for loops in octave)
%for lookup tables
rho_bb = nan(NN,KK,NNN + 1);
rho_cp = nan(NN,JJ,NNN + 1);
beta__ = nan(NN,KK,NNN + 1, length(ang));
Q_ext = nan(NN,JJ,NNN + 1);

%computations
DD =  nan(NNNN,1);
theta_ = nan(NNNN,1);
lambda_bb = nan(NNNN,1);
rr = nan(NNNN,1);
VSF__ = nan(NNNN,length(ang));
VSF_ = nan(NNNN,1);
lambda_cp = nan(NNNN,1);
cp = nan(NNNN,1);

%outputs
VSF = nan(NN, KK, JJ);
d_VSF = nan(NN, KK, JJ);
beam_c = nan(NN, KK, JJ);
d_beam_c = nan(NN, KK, JJ);
ratio = nan(NN, KK, JJ);
unc_ratio = nan(NN, KK, JJ);

for nn=1:NN
    max_D=D0(nn)+1.5*err_D0(nn)+3*delta_D0(nn); %99th percentile.
    min_D=D0(nn)-1.5*err_D0(nn)-3*delta_D0(nn); %99th percentile.
    if min_D<=0
          disp('cannot have negative D')
          return
    end
    for k=1:KK %loop over bb wavelengths
        max_wl_bb=wl(k)+3*delta_wl(k); %99th percentile.
        min_wl_bb=wl(k)-3*delta_wl(k); %99th percentile.
        min_D_over_lambda_bb=min_D/max_wl_bb;
        max_D_over_lambda_bb=max_D/min_wl_bb;
        if min_D_over_lambda_bb<=0
            disp('cannot have negative D/lambda_bb')
            return
        end
        dd=(max_D_over_lambda_bb-min_D_over_lambda_bb)/NNN;
        for i=1:NNN+1 %compute VSF/cross-section for all discretized values
            [n, nm]=IoR(wl(k),T);
            rho_bb(nn,k,i)=pi*nm*(min_D_over_lambda_bb+(i-1)*dd); %rho
            [S1 S2 Qb Qc Qback]=fastmie(rho_bb(nn,k,i),n,nang);
            S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
            S11_rad=S11.*sin(ang'/180*pi)*2*pi;
            S11_int=integrate(S11_rad,dang,length(S11_rad));
            beta__(nn,k,i,:)=Qb*S11/S11_int; %table of VSF/cross-section
            if prod(S11)<0
                S11
                disp('cannot have negative VSF')
                return
            end
        end
    end
    for jj=1:JJ %loop over beam_c wavelengths
        max_wl_cp=c_wl(jj)+3*delta_c_wl(jj); %99th percentile.
        min_wl_cp=c_wl(jj)-3*delta_c_wl(jj); %99th percentile.
        min_D_over_lambda_cp=min_D/max_wl_cp;
        max_D_over_lambda_cp=max_D/min_wl_cp;
        if min_D_over_lambda_cp<=0
            disp('cannot have negative D/lambda_cp')
            return
        end
        dd=(max_D_over_lambda_cp-min_D_over_lambda_cp)/NNN;
        for i=1:NNN+1 %compute beam-c/cross-section for all discretized values
            [n, nm]=IoR(c_wl(jj),T);
            rho_cp(nn,jj,i)=pi*nm*(min_D_over_lambda_cp+(i-1)*dd); %rho
            [S1 S2 Qb Qc Qback]=fastmie(rho_cp(nn,jj,i),n,nang);
            S11=0.5*((abs(S1)).^2 + (abs(S2)).^2);
            S11_rad=S11.*sin(ang'/180*pi)*2*pi;
            S11_int=integrate(S11_rad,dang,length(S11_rad));
            %computations for beam attenuation:correcting for acceptance angle of ac-9
            delta=(pi-accept_angle*pi/180)/length(ang);
            ang_=[accept_angle*pi/180:delta:pi]*180/pi;
            nang_=length(ang_);
            S11_=interp1(ang,S11,ang_,'pchip')';
            S11_rad_=S11_.*sin(ang_'*pi/180)*2*pi; %taking care of azymuthal weighing
            S11_int_=integrate(S11_rad_,delta,nang_);
            Qb_corr=Qb*S11_int_/S11_int; % correct Qb for acceptance angle.
            Qa=Qc-Qb;
            Q_ext(nn,jj,i)=Qa+Qb_corr;
            if Q_ext(nn,jj,i)<=0
                disp('cannot have negative cp')
                return
            end
        end
    end
end
disp('finished computing the lookup tables')
% now that we have the lookup tables we can use them:
% since angle distribution is the same for all wavelengths, no need to
% compute more than once.
theta_=normrnd(normrnd(theta,err_theta),d_theta,1,NNNN);
I=find(theta_>180);
theta_(I)=theta_(I)-180;
if prod(theta_)<=0
    disp('cannot have negative theta')
    return
end
for nn=1:NN
    DD=normrnd(normrnd(D0(nn),err_D0(nn)/2,1,NNNN),delta_D0(nn),1,NNNN);
    if min(DD)<=0
        disp('cannot have negative DD')
        return
    end
    for k=1:KK %loop on bb wavelengths
        lambda_bb=normrnd(wl(k),delta_wl(k),1,NNNN);
        if min(lambda_bb)<=0
            disp('cannot have negative lambda_bb')
            return
        end
        [n, nm]=IoR(lambda_bb,T);
        rr=pi*nm.*DD./lambda_bb;
        VSF__=interp1(squeeze(rho_bb(nn,k,:))',squeeze(beta__(nn,k,:,:)),rr,'linear','extrap');
        %angular averaging
        for i=1:NNNN
            VSF_(i)=pi*DD(i)^2/4*mean(interp1(ang',VSF__(i,:),theta_,'pchip')); %VSF value averaged over the random angles chosen
        end
        for jj=1:JJ
            %sample randomly from wavelength, bead size and angle space based on bead, angle and wavelength distribution.
            lambda_cp=normrnd(c_wl(jj),delta_c_wl(jj),1,NNNN);
            if min(lambda_cp)<=0
                disp('cannot have negative lambda_cp')
                return
            end
            [n, nm]=IoR(lambda_cp,T);
            rr=pi*nm.*DD./lambda_cp;
            cp=(pi*DD.^2/4.*interp1(squeeze(rho_cp(nn,jj,:)),squeeze(Q_ext(nn,jj,:)),rr,'linear','extrap'))'; %compute actual beam_c
            VSF(nn,k,jj)=mean(VSF_);
            d_VSF(nn,k,jj)=std(VSF_);
            beam_c(nn,k,jj)=mean(cp);
            d_beam_c(nn,k,jj)=std(cp);
            ratio(nn,k,jj)=mean(VSF_./cp);
            unc_ratio(nn,k,jj)=std(VSF_./cp);
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
np=1.5718 + 0.008412./(wl.^2) + 0.000235./(wl.^4); %Jones et al.
ni=0.0003;  %imaginary part of index of refraction
n0=1.31405; n1=1.779e-4; n2=-1.05e-6; n3=1.6e-8; n4=-2.02e-6; n5=15.868; n6=0.01155; n7=-0.00423; n8=-4382; n9=1.1455e6; %index of refraction of water Quan and Fry, 1995
nm=n0+n4*T^2+(n5+n7*T)./(wl*1000)+n8./(wl*1000).^2+n9./(wl*1000).^3; %checked
n=(np+ni*sqrt(-1))/nm; %index of refraction relative to water
return
end

function intgr=integrate(y,dx,N)
%function to perform an integral of y with respect to x, with dx the
%increment and N the total number of y's to integrate.
intgr=3/8*(y(1)+y(N))+7/6*(y(2)+y(N-1))+23/24*(y(3)+y(N-2));
intgr=(intgr+sum(y(4:N-3)))*dx;
end
