clear
close all

%first plot our DEER Rabi for reference. Code in this section should be
%identical to calculate_ourDEERRabi_GaussianfromBlochEqns.m

omegaL = 2*pi*(10.7e6)*(225*1e-4);
curlyA = 0.5;

t_pi = 100e-9;
Omega = pi/t_pi;
tau = 900e-9;
T1 = 20e-6;%4*T2;
T2 = 0.005e-6;
curlyB = 1e12;
C0 = 0.3;

Dmat = [0,0,0;...
    0,0,-Omega;...
    0,Omega,0];
Gammamat = [-1/T2,0,0;...
    0,-1/T2,0;...
    0,0,-1/T1];
Wmat = Dmat + Gammamat;
invWmat = inv(Wmat);
invGmat = inv(Gammamat);

Tsweep = linspace(0,900e-9,51);
normalized_signal = zeros(size(Tsweep));

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

    normalized_signal(ii) = C0/2*exp(-curlyB*meansqphase/2 -curlyA/2);

end

figure
plot(Tsweep,normalized_signal)
figure
plot(Tsweep,normalized_signal-normalized_signal(1))
%%

Ts = t_pi;

expWt = expm(Wmat*Ts/2);
expGt = expm(Gammamat*(tau-Ts/2));

ChiI = ( expGt-eye(3) )*invGmat;
ChiII = ( expWt-eye(3) )*invWmat*expGt;
ChiIII = ( expWt-eye(3) )*invWmat*expWt*expGt;
ChiIV = ( expGt-eye(3) )*invGmat*expWt*expWt*expGt;

Chitot = ChiI + ChiII - ChiIII - ChiIV;
meansqphase12 = norm((Chitot')*[0;0;1])^2;

Tcorr_sweep = linspace(1e-9,100e-6,5000);
correlation_signal = zeros(size(Tcorr_sweep));


for nn = 1:length(correlation_signal)

    Tcorr = Tcorr_sweep(nn);

    V = expm(Gammamat*Tcorr)*expGt*expWt*expGt*expWt;
    Vprime = expm(Gammamat*(Tcorr-Ts)/2)*expWt*expm(Gammamat*(Tcorr-Ts)/2)*expGt*expWt*expGt*expWt;
    
    meansqphase34 = norm((V')*(Chitot')*[0;0;1])^2;
    meansqphase34prime = norm((Vprime')*(Chitot')*[0;0;1])^2;
    crossphase = [0,0,1]*Chitot*V*(Chitot')*[0;0;1];
    crossphaseprime = [0,0,1]*Chitot*Vprime*(Chitot')*[0;0;1];

    correlation_signal(nn) = -C0/4*exp(-curlyB*meansqphase12/2 - curlyA/2)*(...
        exp(-curlyB*meansqphase34prime/2 -curlyA/2)*exp(curlyB*crossphaseprime + curlyA*cos(omegaL*(Tcorr+2*tau)))+...
        - exp(-curlyB*meansqphase34prime/2 -curlyA/2)*exp(-curlyB*crossphaseprime - curlyA*cos(omegaL*(Tcorr+2*tau)))+...
        - exp(-curlyB*meansqphase34/2 -curlyA/2)*exp(curlyB*crossphase + curlyA*cos(omegaL*(Tcorr+2*tau))) +...
        exp(-curlyB*meansqphase34/2 -curlyA/2 )*exp(-curlyB*crossphase - curlyA*cos(omegaL*(Tcorr+2*tau))) );

end
%%
figure
plot(Tcorr_sweep*1e6,-correlation_signal)
ylim([-max(abs(normalized_signal-normalized_signal(1))) max(abs(normalized_signal-normalized_signal(1)))])