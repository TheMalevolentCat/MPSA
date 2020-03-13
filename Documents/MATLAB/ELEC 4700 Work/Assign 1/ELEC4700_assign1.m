close all
clc, clear

global c

c.m_0 = 9.10938215e-31; % electron mass (kg)
c.m_n = c.m_0 * 0.26; %effective e mass (kg)
c.q_0 = 1.60217653e-19;  % electron charge (q)
c.k = 1.38e-23; %boltzman constant in S.I. units
c.T = 300; %Tempurature in Kelvin

pNum = 100; %Number of particles
obsNum = 5; %Number of particles you want to observe the motion of
pVel = zeros(pNum, 2); %Velocity matrix col 1 is X velocity col 2 is y velocity
pPos = zeros(pNum, 2); %Position matrix col 1 is X pos col 2 is y pos
Vefcn = @(T) sqrt(3 .* c.k .* T ./ c.m_n); %function for avg e velocity given tempurature

timesteps = 1000;
dt = 0.75e-15; %sec
TVec = zeros(timesteps, 1);
Pscat = 1 - exp(-dt / 0.2e-12);

pdxy = zeros(pNum, 2);
BigX = zeros(pNum, timesteps);
BigY = zeros(pNum, timesteps);

MFPmat = zeros(pNum, timesteps);

%Histogram Creation
mu = 0;
sigma = sqrt(c.k * c.T / c.m_n);
bd = makedist('Normal', 'mu', mu, 'sigma', sigma); %distribution for velocity is stated as normal in the slides

sampNum = pNum * 2;
probVecX = zeros(sampNum, 1);
probVecY = zeros(sampNum, 1);

for i = 1:sampNum
    probVecX(i) = rand;
    probVecY(i) = rand;
end

SampXVecVel = zeros(sampNum, 1);
SampYVecVel = zeros(sampNum, 1);
SampTotalVec = zeros(sampNum, 1);

SampXVecVel = icdf(bd, probVecX);
SampYVecVel = icdf(bd, probVecY);
SampTotalVec = (SampXVecVel.^2 + SampYVecVel.^2).^(1/2);

figure(1)
histogram(SampTotalVec, 100);
xlabel('Velocity (m/s)');
ylabel('Number of Particles');

%Simulation Calculations
for n = 1: pNum
  pVel(n,1) = SampXVecVel(n);
  pVel(n,2) = SampYVecVel(n);
  
  pPos(n,1) = rand .* 200e-9; %apply randomness to position
  pPos(n,2) = rand .* 100e-9;
  
  while (pPos(n,1) <= 120e-9) & (pPos(n,1) >= 80e-9) & (pPos(n,2) >= 60e-9) | (pPos(n,1) <= 120e-9) & (pPos(n,1) >= 80e-9) &  (pPos(n,2) <= 40e-9)
      pPos(n,1) = rand .* 200e-9; %apply randomness to position
      pPos(n,2) = rand .* 100e-9;
  end
  
end

%Actual Simulation
for i = 1: timesteps
    
 pdxy = pVel .* dt;    
 pPos = pPos + pdxy;
 
 TotalVelSqr = pVel(:,1).^2 + pVel(:,2).^2;
 TAvg = mean(TotalVelSqr) * c.m_n / (2 * c.k);
 TVec(i) = TAvg;
 
 for n = 1:pNum %this is the inclusion of scattering, can be commented out to remove scattering
     if Pscat > rand
         pVel(n,1) = SampXVecVel(pNum + n);
         pVel(n,2) = SampYVecVel(pNum + n);
         MFPmat(n,i) = sqrt(pVel(n,1)^2 + pVel(n,2)^2);
     end
 end
 
 if i > 1 %this is the inclusion of the box, can be commented out to remove box
     for m = 1:pNum
        if (pPos(m,1) <= 120e-9) & (pPos(m,1) >= 80e-9) & (pPos(m,2) >= 60e-9) %upper box
           if (BigX(m, i - 1) < 80e-9) & (BigY(m, i - 1) > 60e-9) %collision with left side
               pVel(m,1) = -pVel(m,1);
           elseif (BigX(m, i - 1) > 120e-9) & (BigY(m, i - 1) > 60e-9) %collision with right side
               pVel(m,1) = -pVel(m,1);
           elseif (BigX(m, i - 1) < 120e-9) & (BigX(m, i - 1) > 80e-9) %collision with bottom side
               pVel(m,2) = -pVel(m,2);
           end
        elseif (pPos(m,1) <= 120e-9) & (pPos(m,1) >= 80e-9) &  (pPos(m,2) <= 40e-9) %lower box
           if (BigX(m, i - 1) < 80e-9) & (BigY(m, i - 1) < 40e-9) %collision with left side
               pVel(m,1) = -pVel(m,1);
           elseif (BigX(m, i - 1) > 120e-9) & (BigY(m, i - 1) < 40e-9) %collision with right side
               pVel(m,1) = -pVel(m,1);
           elseif (BigX(m, i - 1) < 120e-9) & (BigX(m, i - 1) > 80e-9) %collision with top side
               pVel(m,2) = -pVel(m,2);
           end
        end
     end
 end
 
 for m = 1:pNum %for loop put in because it was the only way I found the border behaviour works well
     if(((pPos(m,1)) >= (200e-9 - dt * Vefcn(300))) | ((pPos(m,1)) <= (- dt * Vefcn(300)))) %behaviour at x border
        pPos(m,1) = 0;
     end    

     if((pPos(m,2) >= (100e-9 - dt * Vefcn(300))) | (pPos(m,2) <= (dt * Vefcn(300)))) %behaviour at y border
        pVel(m,2) = pVel(m,2)* -1;
     end 
 end
 
 BigX(:,i) = pPos(:, 1);
 BigY(:,i) = pPos(:, 2);
  
end

%Tempurature Plot
figure(2)
plot(1:timesteps, TVec);
ylabel('Tempurature (K)');
xlabel('Timesteps (.75e-14 sec)');

%MFP and MFT calculations

MFPtime = zeros(pNum, timesteps);
totalScats = 0;

for j = 1 : pNum
    count = 0;
    for i = 1 : timesteps
        count = count + 1;
        if MFPmat(j, i) > 0
            totalScats = totalScats + 1;
            MFPtime(j,i) = count;
            count = 0;
        end
    end
end
MFPtime = MFPtime .* dt;
MFPdist = MFPtime .* MFPmat;
finalMFP = sum(MFPdist, 'all') / totalScats;
finalMFT = sum(MFPtime, 'all') / totalScats;

%{
%Trajectory Animation
for i = 1: timesteps %loop to plot what the trajectory looks like in time
    for j = 1: obsNum %no code to connect lines because it creates lines from where particle crosses x boundary
         
        figure(2)
        plot(BigX(j,i), BigY(j,i), 'r.');
        xlim([0 200e-9]);
        ylim([0 100e-9]);
        xlabel('Distance in Meters (M)');
        ylabel('Distance in Meters (M)');
        hold on;
        
        pause(0.001);
    end 
end 
%}

for j = 1: pNum %Trajectory plot
                        
    figure(3)
    plot(BigX(j,:), BigY(j,:), '.', 'color', [rand rand rand]);
    xlim([0 200e-9]);
    ylim([0 100e-9]);
    xlabel('X Distance in Meters (M)');
    ylabel('Y Distance in Meters (M)');
    hold on;

end 

figure(4) %Density Plot
hist3(pPos, 'Nbins', [20 10], 'CDataMode','auto','FaceColor','interp');
xlabel('X Distance in Meters (M)');
ylabel('Y Distance in Meters (M)');
zlabel('Number of electrons');
fprintf('The mean free path is : %d metres \n', finalMFP);
fprintf('The mean free time is : %d seconds \n', finalMFT);

 %3D Tempurature Plot

TempZVec = (pVel(:,1).^2 + pVel(:,2).^2) ./ (3*c.k) .* c.m_n;
x = linspace(0, 200e-9, 1000);
y = linspace(0, 100e-9, 500);
[Xi, Yi] = meshgrid(x, y);
X = pPos(:, 1);
Y = pPos(:, 2);
Z = TempZVec;
Zi = griddata(X, Y, Z, Xi, Yi);

figure(5)
surf(Xi,Yi,Zi);
shading interp
xlabel('X Distance in Meters (M)');
ylabel('Y Distance in Meters (M)');
zlabel('Tempurature of electrons');
