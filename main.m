clear; clc; close all;
% Lidar Final Project 2020 - UV fluorescence lidar for the detection of
% Biological aerosols
% Loo Ting Tan
% lota0178@colorado.edu/ tanlooting@hotmail.com

pi=3.14159265;
h=6.626070040E-34;
c=2.99792458e8;                 % light speed (m/s)
Oss=0.641;                      % Oscillator strength for Na D2 line
kB=1.3806508e-23;               % Boltzmann constant (Unit:  J/K)
Qe=1.60217733e-19 ;

%laser constants
l_e = 60E-6; %mJ
laser_wl = 355E-9; %nm
d_wl = 5; %filter bandwidth [nm]
r_A= pi()*(0.81/2)^2; %telescope area [m]
fov = 0.5; %mrad
sigma =1E-12; % average fluorescence cross section
T_trans=0.8;
R_trans =0.9*0.85*0.9*0.9;
phi = 0.1; %quantum fluorescence yield
csp = 500; % [s-1] dark noise current
qde= 0.4;

r_step=0.5;
range = 1:r_step:1000; %2km
peak_loc= 100:100:900;
for po = 1:length(peak_loc)
    % number density
    dep= 20; %20m

    rho= 10E6; % cloud density 
    start_i=find(range==peak_loc(po));
    n_d= zeros(size(range));
    n_d(start_i:start_i+dep/r_step -1)= rho/(dep/r_step); % assume gaussian

    step = 0.5;
    wl = 100:step:1000;
    wl2= 300:step:800;
    % linear function of UV fingerprint

    d_w=10;
    %wle=450;
    %l_wl= 2/(d_w)*sqrt(log(2)/pi())*exp(-4*log(2)*((wl - wle)/d_w).^2);

    %% normalized fluorescence spectra - trytophan
    vplus = 1e7/23.86e3;
    vminus = 1e7/19.29e3;
    vmax = 1e7/21.80e3;
    p=(vmax-vminus)/(vplus-vmax);
    hf = vplus-vminus;
    a= vmax+hf*p/(p^2-1);
    l_wl = exp(-(log(2)/(log(p)^2))*abs(log((a-wl)/(a-vmax))).^2);

    %%
    for j=1:length(wl2) %emission wavelength

        cum_i=find(wl==wl2(j));
        cum_l_wl(:,j)=sum(l_wl(cum_i-d_w/2/step:cum_i+d_w/2/step))*d_w;
        %plot(wl, l_wl);

        a_m_l = 8*pi()/3*1.54e-3*exp(-range/7)*((532/(laser_wl*1E9))^4);
        a_am_l = 50*((2.47*1E-3*exp(-range/2)+5.13E-6*exp(-(range-20).^2/36)))*532/(laser_wl*1E9);
        a_m_e = 8*pi()/3*1.54e-3*exp(-range/7)*((532/(wl2(j)))^4);
        a_am_e = 50*((2.47*1E-3*exp(-range/2)+5.13E-6*exp(-(range-20).^2/36)))*532/(wl2(j));
        total_u = a_m_l+a_am_l;
        total_d=a_m_e + a_am_e;
        for k = 1: length(range)
                trans_t_u(k)= exp(-sum(total_u(1:k))*r_step);
                trans_t_d(k,j)= exp(-sum(total_d(1:k))*r_step);
        end
        %extinction spectra for NADH
    end

    a_c_l= exp(-16900*1E2*rho*dep/(6.022E23));

    %% Total Photon Count for each wavelength
    dcs = (phi*sigma*cum_l_wl).*(n_d)';
    ph_c= ((l_e*laser_wl/h/c)*(r_A./(4*pi()*range.^2))'*(T_trans*R_trans*qde).*(dcs).*trans_t_u').*(trans_t_d)*a_c_l;

    for i = 1:length(wl2)
        tot_ph_c(po,i)= nansum(ph_c(:,i));
    end
end
%figure(2);
%plot(wl2,tot_ph_c);xlabel("Wavelength(nm)"); ylabel("Total Photon Count (excl. bg and dark current noise");
%% Daytime background
d_T = 2*dep/c;

spec=xlsread("C:\Users\tanlo\Documents\ASEN\Sem3 - Lidar Remote Sensing\HW\project\spectral_irrad.xlsx");
wl_s=spec(:,1)*1000;
s_irrad = spec(:,2);
s_irrad_interp= interp1(wl_s,s_irrad,wl2);
figure(3);
subplot(2,1,1);
plot(wl2,s_irrad_interp);xlabel("Wavelength (nm)"); ylabel("Spectral Irradiance W/m^2/ \mum")
nb = s_irrad_interp.*wl2*1E-6*(fov*0.001/2)^2*pi()*qde*r_A*d_T*(d_wl*1E-9)*R_trans/h/c/4/pi();
subplot(2,1,2);
plot(wl2,nb); xlabel("Wavelength (nm)"); ylabel("Daytime Background Photon Count");

% dark current noise
nd = csp*d_T;
%% Total Signal
[row, col] = size(tot_ph_c);
for m = 1:row
    Nt(m,:)=poissrnd(tot_ph_c(m,:))+poissrnd(nb)+nd;
end
figure(4);
plot(wl2,Nt(10,:)); hold on; plot(wl2,nb);  xlabel("Wavelength (nm)"); ylabel("Photon count at 1900m"); legend("Photon Count at 1900m", "Background photon count");
%% Sensitivity analysis
shotsnum=1;
figure(5);
for b= 1:length(peak_loc)
    snr(b,:)= sqrt(shotsnum)*Nt(b,:)./sqrt(Nt(b,:)+2*(nb+nd));

end
plot(wl2,snr); xlabel("Wavelength (nm)"); ylabel("SNR at night-time");

