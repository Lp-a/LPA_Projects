
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_April16_xy1_after488_in_mea\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_April3_xy1_DEER_MEA_after488_trail2\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_March12_DEER_after488_MEA50mM\';
path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_May14_xy3_NV1_DEER_after488_20mW_no_OD_dye_inPMMA_with_RF_in_MEA\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_March12_DEER_after488_MEA50mM\';
%path ='C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_May15_xy2_single_NV_after488_20mW_no_OD_DEER_in_MEA_with_RF_on\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\20250616 functionalized diamond\r6NV2 DEER Rabi 655 MHz -4p2 dBm after 488 overnight\';
%path ='C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_May15_xy2_single_NV_after488_20mW_no_OD_DEER_in_MEA_with_RF_on\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\20250617 functionalized diamond\rA1NV1 DEER Rabi 655 MHz -4p2 dBm after 488 no OD 1 min\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\20250617 functionalized diamond\rC1NV1 DEER Rabi 655 MHz -4p2 dBm after 488 OD 2 1 min\';
%path = 'C:\Users\lpa3\Box\Backlund lab shared\data\ExperimentEditor Data\2025_April1_xy1_NV1_DEER_after488\';


file = 'Iteration_336';%530, 'Iteration_171'

load([path,file])

reference = savedData.reference;

signal = savedData.signal;

contrast = (reference-signal)./reference;


gammae = 2*pi*2.8e10;

h = 6.626e-34;

hbar = h/2/pi;

dNV = 12e-9;



Tsweep_start = (1e-9)*savedData.loadedScan.starts(1); %seconds

Tsweep_end = (1e-9)*savedData.loadedScan.ends(1);

nsteps = savedData.loadedScan.nsteps;

Tsweep_vec = linspace(Tsweep_start,Tsweep_end,nsteps)';

twait_start = (1e-9)*savedData.loadedScan.starts(2);

tau = Tsweep_start + twait_start;



prefactorguess = 5e12;

T1guess = 100e-6;%100e-6

T2guess = 200e-6;%0.005e-6

t_piguess = Tsweep_vec(contrast==min(contrast));

C0guess = max(contrast);

Aguess = 1e6;

Deltaguess = 2*pi*1e6;%2*pi*1e6;or 0



paramsGaussian = lsqnonlin(@(params) params(1)*contrast_GaussianfromBloch(Tsweep_vec,tau,params(2)*1e-6,params(3)*1e-6,params(4)*1e-6,params(5)*1e11)-contrast,[C0guess,t_piguess*1e6,T1guess*1e6,T2guess*1e6,prefactorguess*1e-11], [min(contrast),0.025,0.001,0.00001,0], [0.3,10,1000,1000,1000]);



prefactorfit = paramsGaussian(5)*1e11;

sigmafit = prefactorfit/( (1e-14)*3*pi/64*(gammae^4)*(hbar^2)/dNV^4 );

T2fit = paramsGaussian(4);
T1fit = paramsGaussian(3);



sigma_pernm = sigmafit*1e-18;

T2_ns = T2fit*1e3;
T1_us = T1fit;



contrastfitGaussian = paramsGaussian(1)*contrast_GaussianfromBloch(Tsweep_vec,tau,paramsGaussian(2)*1e-6,paramsGaussian(3)*1e-6,paramsGaussian(4)*1e-6,paramsGaussian(5)*1e11);



paramsSingle = lsqnonlin(@(params) params(1)*contrast_single(Tsweep_vec,tau,params(2)*1e-6,params(3)*1e6,params(4)*2*pi*1e6)-contrast, [C0guess,t_piguess*1e6,Aguess*1e-6,Deltaguess/2/pi*1e-6],[min(contrast),0.025,0,0], [0.3,10,100,100]);



contrastfitSingle = paramsSingle(1)*contrast_single(Tsweep_vec,tau,paramsSingle(2)*1e-6,paramsSingle(3)*1e6,paramsSingle(4)*2*pi*1e6);



%contrastfitGaussian = C0guess*contrast_GaussianfromBloch(Tsweep_vec,tau,t_piguess,T1guess,T2guess,prefactorguess);

%contrastfitSingle = C0guess*contrast_single(Tsweep_vec,tau,t_piguess,Aguess,Deltaguess);



figure('color','w')

plot(Tsweep_vec*1e9,contrast,'linewidth',2)

hold on

plot(Tsweep_vec*1e9,contrastfitGaussian,'linewidth',2)

%plot(Tsweep_vec*1e9,contrastfitSingle,'linewidth',2)

set(gca,'fontsize',14)

legend('Data','fit to Gaussian ensemble','fit to single target spin','location','southeast')

xlabel('T_{sweep} (ns)')

ylabel('Contrast (a.u.)')

fprintf('--- Fit Results ---\n');
fprintf(' T1 : %.2f us\n', paramsGaussian(3));
fprintf('T2:  %.4f ns\n', paramsGaussian(4)*1e3);



function mycontrast = contrast_GaussianfromBloch(Tsweep,tau,t_pi,T1,T2,prefactor)



    Omega = pi/t_pi;

    Dmat = [0,0,0;...

        0,0,-Omega;...

        0,Omega,0];

    Gammamat = [-1/T2,0,0;...

        0,-1/T2,0;...

        0,0,-1/T1];

    Wmat = Dmat + Gammamat;

    invWmat = inv(Wmat);

    invGmat = inv(Gammamat);



    mycontrast = zeros(size(Tsweep));



    for ii = 1:length(Tsweep)

        % disp([num2str(ii) '/' num2str(length(Tsweep))])

    

        Ts = Tsweep(ii);

    

        expWt = expm(Wmat*Ts);

        expGt = expm(Gammamat*(tau-Ts));

    

        ChiI = ( expWt-eye(3) )*invWmat;

        ChiII = ( expGt-eye(3) )*invGmat*expWt;

        ChiIII = ChiI*expGt*expWt;

        ChiIV = ChiII*expGt*expWt;

    

        Chitot = ChiI + ChiII - ChiIII - ChiIV;

        Chisq = Chitot*(Chitot');

        meansqphase = Chisq(3,3);

    

        mycontrast(ii) = exp(-prefactor*meansqphase);

    

    end



end



function mycontrast = contrast_single(Tsweep,tau,t_pi,A,Delta)



    Omega = pi/t_pi;



    Omegatilde = sqrt(Omega^2 + Delta^2);

    Omegaprimetilde = sqrt(Omega^2 + (Delta-A)^2);

    

    C2 = cos(Omegaprimetilde*Tsweep/2);

    S2 = sin(Omegaprimetilde*Tsweep/2);

    

    C4 = cos(Omegatilde*Tsweep/2);

    S4 = sin(Omegatilde*Tsweep/2);



    mycontrast = (Omega^2)/Omegatilde/Omegaprimetilde*S2.*S4.*cos(A*(tau-Tsweep/2)) + cos(A*Tsweep/2).*( C2.*C4 + S2.*S4*Delta*(Delta-A)/Omegatilde/Omegaprimetilde ) +-sin(A*Tsweep/2).*( S2.*C4*(Delta-A)/Omegaprimetilde - C2.*S4*Delta/Omegatilde );





end
