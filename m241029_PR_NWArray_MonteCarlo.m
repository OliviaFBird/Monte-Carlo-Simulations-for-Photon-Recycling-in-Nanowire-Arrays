tic %start tracking runtime
%% Key Assumptions
%Cylinder/Wire shaped geometry 
%Photons behave as particles and recycling can be predicted by solid angle
%Photon emission is equally likly to occur at any point on nanowire surface (infinite mobility within wire)

%User required inputs: independent variable range (Line 12), Experimental Absorbance and PL Data (Lines 39-43),
%nanowire geometry parameters (Lines 51-57), radiative rate coefficient (line 61),
%fraction of light expected to scatter radially (lines 64-65), refractive index of material (lines 67-71), 
%Laser excitation parameters (73-76), how many photons to simulate at once (and thus runtime, Line 133)

%% Setup
clear all
VaryParam = linspace(500,3300,20);%Set your desired independent variable
PLRatio = zeros(length(VaryParam),1); %Initialize tracker for the PL Ratio
ScattRatio = zeros(length(VaryParam),1); %Initialize tracker for Reflection ratio
PLTot = zeros(length(VaryParam),1); %Initialize tracker for Total PL 

% Initialize flags for tracking run progress
flag10 = false;
flag20 = false;
flag30 = false;
flag40 = false;
flag50 = false;
flag60 = false;
flag70 = false;
flag80 = false;
flag90 = false;
%% Count Different types of events
AbsE = 0; %Count Absorbtion events
TransE = 0; %Count Transmission events
ScatE = 0; %Count Scattering Events 
NRRec = 0; %Count non-radiative recombination events
RadRec = 0; %Count radiative recombination events 
Lost = 0; %Count light lost through radiative emission that doesn't hit 2nd wire
Impossible = 0; %initialize an impossible case tracker for troubleshooting

%Import Absorbance and PL Spectra as csv
%NOTE: Data must be in the same folder as this code to import properly

Abs = importdata (""); %Import Absorbance in absorbance units/OD
PL = importdata(""); %Import PL spectrum. Should be normalized
ExpData = importdata(""); %Import experimental IRF (optional)
AbsMin = 206; %Minimum wavelength of absorbance spectrum
AbsMax = 1111; %Maximum wavelength of absorbance spectrum
%NOTE: This code assumes that the range of your PL spectrum is less than or
%equal to the range of your absorbance spectrum. Interpolation may not work
%properly if this is not true

for j = 1:length(VaryParam)
    %Nanorod/wire parameters
    lengthNW = 2000; %length of nanowire in nm
    rad = 250; %radius of nanowire in nm 
    EEd = 800 ;
    d = VaryParam(j); %center-to-center interparticle distance in nm (must be greater than 2*rad)
    
    PLQY = 0.5; %Quantum yield of emission as fraction
    
    VNW = (pi*rad^2*lengthNW)/((10^7)^3); %Volume of NW in cm^3
    
    RadCoef = 1.2E-10; %Radiative Rate Coefficient in cm^3/s Source: https://www.ioffe.ru/SVA/NSM/Semicond/InP/electric.html
    kRad = RadCoef/VNW; % radiative recombination rate constant (in 1/s)
    kNRad = (kRad/PLQY)-kRad; %non-radiative recombination rate constant (in 1/s)
    SideScattExt = 0.1; %estimated fraction of scattered light to go to the side from excitation
    SideScatt = 0.1; %scattered light scattered sideways if the light came from another wire

    %RefIn = 4.3784; %refractive index of wires at excitation wavelength
    %ExtScatt = (abs((1-RefIn)/(1+RefIn)))^2; %estimated fraction of light scattered from Fresnel Equation
    ExtScatt = 0.41; %Fraction of light scattered experimentally
    RefInPL = 3.4138; %refractive index of wires at PL peak wavelength
    ExtScattPL = (abs((1-RefInPL)/(1+RefInPL)))^2; %fraction of PL scattered

    %Input information on your excitation beam 
    Exnm = 415; %wavelength of excitation beam (nm) 
    ExP = 0.3/1000 ; %Power of excitation laser (Watts)
    Ext = 1 ; %Duration of excitation (s)
    
    %% Static Calculations 
    h = 6.62607015*10^(-34); %plank's constant
    c = 299792458; %speed of light
    
    AbsY = Abs.data(:,2); %Assign 2nd column of absorbance data as your Y data (OD). If your data is plotted in a different column change the "2" to the appropriate column
    AbsYSide = AbsY*2*rad/lengthNW;
    AbsXnm = Abs.data(:,1); %Assign X axis (nm) data to first column. Again, change the "1" to the appropriate column if your data is not in the same form. 
    AbsX = h*c./(AbsXnm/10^9)*6.242E18; %Convert absorbance x-axis to eV
    PLXnm = PL.data(:,1); %Assign PL X-axis (nm) as 1st column
    PLX = h*c./(PLXnm/10^9)*6.242E18; %Convert PL x-axis to eV
    PLY = PL.data(:,2); %Assign PL Y-axis (counts) to 2nd column
    PLYNorm = normalize(PLY,"norm",1); %normalize PL y-axis
    PLint = trapz(PLYNorm); %Confirm normalization worked properly. PLint = 1
    
    
    % Interpolate the PL to match the x-values of the Abs dataset
    AbsMinev = h*c./(AbsMin/10^9)*6.242E18; %convert AbsMin into ev
    AbsMaxev = h*c./(AbsMax/10^9)*6.242E18; %convert AbsMax into ev
    PLX(end+1) = AbsMinev; %extend PL spectrum to match Abs spectrum (change the '6' to the maximum value of your absorbance spectrum)
    PLX(end+1) = AbsMaxev; %extend PL spectrum to match Abs spectrum (change the '6' to the maximum value of your absorbance spectrum)
    PLYNorm(end+1) = 0; %add 0 to end of PL spectrum out to end of Abs spectrum
    PLYNorm(end+1) = 0; %add 0 to end of PL spectrum out to beginning of Abs spectrum
    PLY_interp = normalize(interp1(PLX, PLYNorm, AbsX, 'linear', 'extrap'), "norm", 1); %Interpolate PL spectrum to Abs spectrum
    PLinterp_int = trapz(PLY_interp); %confirm still normalized after interpolation. Should = 1
    
    %Multiply PL spectrum by Absortance to get overlap 
    Overlap = PLY_interp.*(1-10.^(-AbsY));
    OverlapIntegral = trapz(Overlap); %integrate overlap spectrum. This is equal to the fraction of emitted photons which could be absorbed
    %% Optional: Plot PL, abs, and overlap to ensure all data was uploaded and normalized correctly
    %Plot Abs, PL, and overlap
    %plot(AbsX,AbsY),
    %hold on
    %plot (AbsX,PLY_interp)
    %plot(AbsX,Overlap)
    %xlim([1.1 2])
    %ylim([0, 1])
    %legend
    %hold off
    %%
    phoH=lengthNW/2; %Set an initial value for photon height, will vary each loop
    
    NWSA = 2*pi*rad*lengthNW+2*pi*rad^2;%Nanowire Surface Area
    SideSA = 2*pi*rad*lengthNW; %surface area of NW excluding top and bottom
    XAngle = 2*atan(rad/d); %What angle in the XY plane will result in PR to adjacent NWs
    ZAngle = pi-atan((d-rad)/(lengthNW-phoH))-atan((d-rad)/(phoH)); %What angle in Z will result in PR to adjacent NWs
    
    XChance = XAngle/(2*pi);%divide occluded angle by total angle to get probability 
    ZChance = ZAngle/pi; %divide occluded angle by total angle to get probability 
    
    %PLQY = kRad/(kRad + kNRad); %calculate PL quantum yield from rate constants
    
    EPho = h*c/(Exnm/10^9); %Calculate Energy of a photon in J
    NPhoTot = ExP*Ext/EPho; % Calculate Total number of photons in experiment
    OnePhoRate = (ExP/EPho)^(-1); %Calculate time between photons in s 
    
    PhoBunchSize = 1E7; %how many photon movements do you want to simulate at once? (Change this variable to control simulation runtime. Larger bunch = faster runtime but worse signal : noise)
    NPho = NPhoTot / PhoBunchSize; %total photon bunches to simulate 
    PhoRate = (EPho*PhoBunchSize)/ExP; %Calculate time between photon bunches in s 
    kPho = PhoRate^(-1);
    
    %find how much light is intially absorbed from the laser 
    ExtVal = interp1(AbsXnm,AbsXnm,Exnm,'nearest'); %Find value of Abs spectrum closest to excitation wavelength
    ExtInd = find(AbsXnm == ExtVal); %Find index of ExtVal
    ExtTrans = 10.^(-AbsY(ExtInd)); %Fraction of laser light transmitted
    ExtTransSide = 10.^(-AbsYSide(ExtInd)); %Fraction of laser light transmitted

    ExtAbs = 1-ExtTrans; %Fraction of laser light absorbed from top
    ExtAbsSide = 1-ExtTransSide;%fraction of laser light absorbed from side
    PLAbsSide = OverlapIntegral;%fraction of PL light absorbed from side
    
    %%
    RunArray = zeros(1,2); %Initialize Array to track photons at each wire
    DetectionArrayPL = zeros(1,2); %Initialize array to track all detected emission events
    DetectionArrayScatt = zeros(1,2); %Initialize array to track all detected scattering events 
    %%
    while NPho > 0 %Run until all excitations have been exhausted
        if RunArray(1,1) > 0 && RunArray(1,2) > 0 %Both wires are excited
            ktot = (kPho + kNRad*RunArray(1,1)*PhoBunchSize + kRad*RunArray(1,1)*PhoBunchSize + kNRad*RunArray(1,2)*PhoBunchSize + kRad*RunArray(1,2)*PhoBunchSize); %define all available pathways
            r = rand; %Draw random number 
            if r < (kPho/ktot) %new excitation 
                NPho = NPho - 1; %next excitation occurs
                r = rand;
                if r < ExtScatt %light scatters
                    r = rand;
                    ScatE = ScatE + 1;
                    if r < SideScattExt %light is scattered outward
                        r= rand;
                        if r < XChance %light hits wire 2
                            r = rand;
                            if r < ExtScatt % light is scattered off wire 2
                            ScatE = ScatE + 1;
                            r = rand;
                                if r < SideScatt %light is scattered outward
                                    Lost = Lost +1;
                                else % light is scattered up
                                    DetectionArrayScatt(1,2) = DetectionArrayScatt(1,2) + 1; %detect light
                                end
                            else %light interacts with wire 2
                                r = rand;
                                if r < ExtAbsSide %light is absorbed by wire 2
                                    RunArray(1,2) = RunArray(1,2) + 1;
                                    AbsE = AbsE + 1;
                                else %light is transmit 
                                    TransE = TransE +1;
                                end
                            end
                        else %light is lost 
                            Lost = Lost + 1;
                        end
                    else %light is scattered up
                        DetectionArrayScatt(1,1) = DetectionArrayScatt(1,1) + 1; 
                    end
                else %light doesn't scatter
                    r = rand;
                    if r < ExtAbs %light is absorbed
                        RunArray(1,1) = RunArray(1,1) + 1;
                        AbsE = AbsE + 1;
                    else %light is transmit
                        TransE = TransE + 1;
                    end
                end
            elseif r > kPho/ktot && r < (kPho + kRad*RunArray(1,1)*PhoBunchSize)/ktot %rad rec at wire 1
                r = rand;
                RadRec = RadRec + 1;
                RunArray(1,1) = RunArray(1,1) - 1;
                if r < (SideSA/NWSA)%Emission from Side
                    rx = rand; 
                    rz = rand;
                    if rx < XChance && rz < ZChance %emission hits wire 2
                        r = rand;
                        if r < ExtScattPL %light is scattered 
                            ScatE = ScatE + 1;
                            r = rand;
                            if r<SideScatt %light is scattered away
                                Lost = Lost + 1;
                            else %light is scattered up
                                DetectionArrayPL(1,2) = DetectionArrayPL(1,2) + 1;
                            end
                        else %light is not scattered 
                            r = rand;
                            if r < PLAbsSideegral %light is absorbed
                                RunArray(1,2) = RunArray (1,2) + 1;
                                AbsE = AbsE + 1;
                            else %light is transmitted
                                TransE = TransE + 1;
                            end
                        end
    
                    else %light is lost
                        Lost = Lost + 1;
                    end
                elseif r > (SideSA/NWSA) && r < (SideSA/NWSA+pi*rad^2) %Emission from Top
                    DetectionArrayPL(1,1) = DetectionArrayPL(1,1) + 1;
                else %Emission from bottom
                    Lost = Lost + 1;
                end
            elseif r > (kPho + kRad*RunArray(1,1)*PhoBunchSize)/ktot && r < (kPho + kRad*RunArray(1,1)*PhoBunchSize + kNRad*RunArray(1,1)*PhoBunchSize)/ktot % nonrad rec at wire 1
                RunArray(1,1) = RunArray (1,1) - 1;
                NRRec = NRRec + 1;
            elseif r > (kPho + kRad*RunArray(1,1)*PhoBunchSize + kNRad*RunArray(1,1)*PhoBunchSize)/ktot && r < (kPho + kRad*RunArray(1,1)*PhoBunchSize + kNRad*RunArray(1,1)*PhoBunchSize + kNRad*RunArray(1,2)*PhoBunchSize)/ktot % nonrad rec at wire 2
                RunArray(1,2) = RunArray (1,2) - 1;
                NRRec = NRRec + 1;
            else % rad rec at wire 2
                RunArray(1,2) = RunArray(1,2) - 1;
                RadRec = RadRec + 1;
                r = rand;
                if r < (1-SideSA/NWSA)/2 %light emits upward
                    DetectionArrayPL (1,2) = DetectionArrayPL(1,2) + 1;
                elseif r > (1-SideSA/NWSA)/2 && r < (1-SideSA/NWSA)  %light emits Down
                    Lost = Lost + 1;
                else %light emits out
                    rx = rand;
                    rz = rand;
                    if rx < XChance && rz < ZChance %emission hits wire 1
                        r = rand;
                        if r < ExtScattPL % light is scattered off wire 1
                            r = rand;
                            ScatE = ScatE + 1;
                            if r < SideScatt %light scatters away
                                Lost = Lost + 1;
                            else %light scatters up
                                DetectionArrayPL(1,1) = DetectionArrayPL(1,1) + 1;
                            end
                        else %light interacts with wire 1
                            r = rand;
                            if r < PLAbsSideegral %light is absorbed
                                RunArray(1,1) = RunArray(1,1) + 1;
                                AbsE = AbsE + 1;
                            else %light is transmitted 
                                TransE = TransE + 1;
                            end
                        end
                    else %light is lost
                        Lost = Lost + 1;
                    end
                end
            end 
        elseif RunArray(1,1) > 0 && RunArray(1,2) == 0 %Central wire is excited, neighboring wire is not
            ktot = (kPho + kNRad*RunArray(1,1)*PhoBunchSize + kRad*RunArray(1,1)*PhoBunchSize);
            r = rand;
            if r < (kPho/ktot) %get another excitation 
                NPho = NPho - 1; %next excitation occurs
                r = rand;
                if r < ExtScatt %light scatters
                    r = rand;
                    ScatE = ScatE + 1;
                    if r < SideScattExt %light is scattered outward
                        r= rand;
                        if r < XChance %light hits wire 2
                            r = rand;
                            if r < ExtScatt % light is scattered off wire 2
                            ScatE = ScatE + 1;
                            r = rand;
                                if r < SideScatt %light is scattered outward
                                    Lost = Lost +1;
                                else % light is scattered up
                                    DetectionArrayScatt(1,2) = DetectionArrayScatt(1,2) + 1; 
                                end
                            else %light interacts with wire 2
                                r = rand;
                                if r < ExtAbsSide %light is absorbed by wire 2
                                    RunArray(1,2) = RunArray(1,2) + 1;
                                    AbsE = AbsE + 1;
                                else %light is transmit 
                                    TransE = TransE +1;
                                end
                            end
                        else %light is lost 
                            Lost = Lost + 1;
                        end
                    else %light is scattered up
                        DetectionArrayScatt(1,1) = DetectionArrayScatt(1,1) + 1; 
                    end
                else %light doesn't scatter
                    r = rand;
                    if r < ExtAbs %light is absorbed
                        RunArray(1,1) = RunArray(1,1) + 1;
                        AbsE = AbsE + 1;
                    else %light is transmit
                        TransE = TransE + 1;
                    end
                end
                    
            elseif r > (kPho/ktot) && r < ((kPho+kRad*RunArray(1,1)*PhoBunchSize)/ktot) %radiative recombination at wire 1
                r = rand;
                RadRec = RadRec + 1;
                RunArray(1,1) = RunArray(1,1) - 1;
                if r < (SideSA/NWSA)%Emission from Side
                    rx = rand;
                    rz = rand;
                    if rx < XChance && rz < ZChance %emission hits wire 2
                        r = rand;
                        if r < ExtScattPL %light is scattered 
                            ScatE = ScatE + 1;
                            r = rand;
                            if r<SideScatt %light is scattered away
                                Lost = Lost + 1;
                            else %light is scattered up
                                DetectionArrayPL(1,2) = DetectionArrayPL(1,2) + 1;
                            end
                        else %light is not scattered 
                            r = rand;
                            if r < PLAbsSideegral %light is absorbed
                                RunArray(1,2) = RunArray (1,2) + 1;
                                AbsE = AbsE + 1;
                            else %light is transmitted
                                TransE = TransE + 1;
                            end
                        end
    
                    else %light is lost
                        Lost = Lost + 1;
                    end
                elseif r > (SideSA/NWSA) && r < (SideSA/NWSA+pi*rad^2) %Emission from Top
                    DetectionArrayPL(1,1) = DetectionArrayPL(1,1) + 1;
                else %Emission from bottom
                    Lost = Lost + 1;
    
                end
            else %nonradiative recombination
                NRRec = NRRec + 1;
                RunArray(1,1) = RunArray(1,1) - 1;
            end
           
        elseif RunArray(1,1) == 0 && RunArray(1,2) > 0 %Neighboring wire is excited, central wire is not
            ktot = (kPho + kNRad*RunArray(1,2)*PhoBunchSize + kRad*RunArray(1,2)*PhoBunchSize);
            r = rand;
            if r < (kPho/ktot) %get another excitation 
                NPho = NPho - 1; %next excitation occurs
                r = rand;
                if r < ExtScatt %light scatters
                    r = rand;
                    ScatE = ScatE + 1;
                    if r < SideScattExt %light is scattered outward
                        r= rand;
                        if r < XChance %light hits wire 2
                            r = rand;
                            if r < ExtScatt % light is scattered off wire 2
                            ScatE = ScatE + 1;
                            r = rand;
                                if r < SideScatt %light is scattered outward
                                    Lost = Lost +1;
                                else 
                                    DetectionArrayScatt(1,2) = DetectionArrayScatt(1,2) + 1; 
                                end
                            else %light interacts with wire 2
                                r = rand;
                                if r < ExtAbsSide %light is absorbed by wire 2
                                    RunArray(1,2) = RunArray(1,2) + 1;
                                    AbsE = AbsE + 1;
                                else %light is transmit 
                                    TransE = TransE +1;
                                end
                            end
                        else %light is lost 
                            Lost = Lost + 1;
                        end
                    else %light is scattered up
                        DetectionArrayScatt(1,1) = DetectionArrayScatt(1,1) + 1; 
                    end
                else %light doesn't scatter
                    r = rand;
                    if r < ExtAbs %light is absorbed
                        RunArray(1,1) = RunArray(1,1) + 1;
                        AbsE = AbsE + 1;
                    else %light is transmit
                        TransE = TransE + 1;
                    end
                end
            elseif r > (kPho/ktot) && r < ((kPho + kRad*RunArray(1,2)*PhoBunchSize)/ktot) %radiative recombination at wire 2
                RunArray(1,2) = RunArray(1,2) - 1;
                RadRec = RadRec + 1;
                r = rand;
                if r < (1-SideSA/NWSA)/2 %light emits upward
                    DetectionArrayPL (1,2) = DetectionArrayPL(1,2) + 1;
                elseif r > (1-SideSA/NWSA)/2 && r < (1-SideSA/NWSA)  %light emits Down
                    Lost = Lost + 1;
                else %light emits out
                    rx = rand;
                    rz = rand;
                    if rx < XChance && rz < ZChance %emission hits wire 1
                        r = rand;
                        if r < ExtScattPL % light is scattered off wire 1
                            r = rand;
                            ScatE = ScatE + 1;
                            if r < SideScatt %light scatters away
                                Lost = Lost + 1;
                            else %light scatters up
                                DetectionArrayPL(1,1) = DetectionArrayPL(1,1) + 1;
                            end
                        else %light interacts with wire 1
                            r = rand;
                            if r < PLAbsSideegral %light is absorbed
                                RunArray(1,1) = RunArray(1,1) + 1;
                                AbsE = AbsE + 1;
                            else %light is transmitted 
                                TransE = TransE + 1;
                            end
                        end
                    else %light is lost
                        Lost = Lost + 1;
                    end
                end
            else %nonradiative recombination 
                RunArray(1,2) = RunArray(1,2) - 1;
                NRRec = NRRec - 1;
            end
        elseif RunArray(1,1) == 0 && RunArray(1,2) == 0 % Neither wire is excited
            NPho = NPho - 1; %next excitation occurs
            r = rand;
            if r < ExtScatt %light scatters
                r = rand;
                ScatE = ScatE + 1;
                if r < SideScattExt %light is scattered outward
                    r= rand;
                    if r < XChance %light hits wire 2
                        r = rand;
                        if r < ExtScatt % light is scattered off wire 2
                        ScatE = ScatE + 1;
                        r = rand;
                            if r < SideScatt %light is scattered outward
                                Lost = Lost +1;
                            else %light is scattered up
                                DetectionArrayScatt(1,2) = DetectionArrayScatt(1,2) + 1; 
                            end
                        else %light interacts with wire 2
                            r = rand;
                            if r < ExtAbsSide %light is absorbed by wire 2
                                RunArray(1,2) = RunArray(1,2) + 1;
                                AbsE = AbsE + 1;
                            else %light is transmit 
                                TransE = TransE +1;
                            end
                        end
                    else %light is lost 
                        Lost = Lost + 1;
                    end
                else %light is scattered up
                    DetectionArrayScatt(1,1) = DetectionArrayScatt(1,1) + 1; 
                end
            else %light doesn't scatter
                r = rand;
                if r < ExtAbs %light is absorbed
                    RunArray(1,1) = RunArray(1,1) + 1;
                    AbsE = AbsE + 1;
                else %light is transmit
                    TransE = TransE + 1;
                end
            end
                
        else 
            Impossible = Impossible + 1;
        end
    end 
    %Calculate and display % complete
    if j > length(VaryParam)*0.9 && ~flag90
        disp('90 %')
        flag90 = true;
    elseif j > length(VaryParam)*0.8 && ~flag80
        disp('80 %')
        flag80 = true;
    elseif j > length(VaryParam)*0.7 && ~flag70
        disp('70 %')
        flag70 = true;
    elseif j > length(VaryParam)*0.6 && ~flag60
        disp('60 %')
        flag60 = true;
    elseif j > length(VaryParam)*0.5 && ~flag50
        disp('50 %')
        flag50 = true;
    elseif j > length(VaryParam)*0.4 && ~flag40
        disp('40 %')
        flag40 = true;
    elseif j > length(VaryParam)*0.3 && ~flag30
        disp('30 %')
        flag30 = true;
    elseif j > length(VaryParam)*0.2 && ~flag20
        disp('20 %')
        flag20 = true;
    elseif j > length(VaryParam)*0.1 && ~flag10
        disp('10 %')
        flag10 = true;
    end
    PLTot(j) = DetectionArrayPL(2)*PhoBunchSize; %Multiply detected bunches by bunch size to get total photons 
    PLRatio(j) = DetectionArrayPL(2)/DetectionArrayPL(1); %Calculate PL ratio by dividing detection array at neighbor by detection at central 
    ScattRatio(j) = DetectionArrayScatt(2)/DetectionArrayScatt(1); %Calculate scattering/reflection ratio
end
toc
%% Outputs
figure(1); %Create plot of PL and reflection ratio 
plot(VaryParam,PLRatio,'r')
hold on
plot(VaryParam,ScattRatio, 'b')
legend('PL Ratio','Reflection Ratio')
hold off
detected = sum(DetectionArrayPL) + sum(DetectionArrayScatt);
figure(2); %create pie chart of loss mechanisms
LossMechNames = ["Detected";"Transmitted";"Missed Wire";"Non-radiative Recombination"];
LossMechNum = [detected; TransE; Lost; NRRec];
tbl = table(LossMechNames,LossMechNum);
pie(LossMechNum)
lgd = legend(LossMechNames);
%%
OutputVary = VaryParam'; %transpose Variable parameter for easier export to graphing programs

