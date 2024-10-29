%% geometric calculations for photon recycling in nanowire arrays. 
%Users should input their independent variable range in Line 4 and input
%wire geometry in lines 9-11
Vary=linspace(50,1000,5000); %Input your independent variable over the desired range (e.g. radius)

GeometricProb = zeros(length(Vary),1); %Initialize output matrix to store values

for j = 1:length(Vary) %Run for all values of your variable parameter  
    rad = Vary(j); %Wire Radius
    h = 2000; %Wire Height
    EEd = 800;%Edge-to-Edge Distance
    d = EEd+2*rad; %Center-to-Center Distance
    
    SideChance = 2*pi*rad*h/(2*pi*rad^2+2*pi*rad*h); %Calculate the probability the photon is emitted from the side (not top or bottom)
    XAngle = 2*atan(rad/d); %Calculate the angle needed to hit the next wire in the XY plane 
    XChance = XAngle/(2*pi); %Calculate the probability of hitting the next wire in the XY plane
    
    n = round(h); %number of Z positions to calculate from
    ZArray = zeros(1,n); %initialize z matrix of all z probabilities
    for i = 1:n %Calculate the probability of being within the correct Z angle at different z heights
        z = i;
        ZAngle = pi-atan((d-rad)/(h-z))-atan((d-rad)/(z)); %calculate angle needed to hit neighboring wire in z plane
        ZChance = ZAngle/pi; %calculate probability of hitting neighboring wire in the z angle
        ZArray(1,i) = ZChance; %Add z probability to array
    end 
    
    ZChanceAvg = mean(ZArray); %Average the probability of being at the correct Z angle over all positions of z
    TotalProb = SideChance*XChance*ZChanceAvg; % calculate the total probability of absorbing a photon at the 2nd wire
    GeometricProb(j,1) = TotalProb; %Add probability of reabsorbtion to output matrix
end 
plot(Vary,GeometricProb) %plot reabsorbtion probability against independent variable