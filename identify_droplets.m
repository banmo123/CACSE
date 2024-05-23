clc;clear;close

%% Identify droplets
I=imread(".tif"); % the picture of droplets
subplot(3,4,[1,2,3,5,6,7,9,10,11]); 
imshow(I);
lower = 50; % lower threshold of the radius size(pixels)
higher = 100; % upper threshold of the radius size(pixels)
[centers,radii] = imfindcircles(I,[lower higher],"method","TwoStage",'ObjectPolarity','dark', ...
    "Sensitivity",0.9,"EdgeThreshold",0.1); % identify droplets
viscircles(centers,radii,'Edgecolor','b','LineWidth',0.1); % show the circles

%%  Export the diameter to Excel
BarScale = 100/155; %micron length corresponding to each pixel
Diameters = radii*2*BarScale;
writematrix(Diameters,".xls");

%%  Plot frequency distribution histograms
Ave_Dia = mean(Diameters); % average diameters of all the droplets identified above
Mu_Dia = std(Diameters); % standard deviation
cv = Mu_Dia/Ave_Dia; % variation coefficient
bins = (lower:1:higher)*BarScale*2;
H=histogram(Diameters,bins,'Normalization','probability'); %概率密度
text(200,0.16,['CV=',num2str(100*cv),'%'],"FontSize",10);

% Plot probability density function curves
hold on
bins=(lower:0.1:higher)*BarScale*2;
f = exp(-(bins-Ave_Dia).^2./(2*Mu_Dia^2))./(Mu_Dia*sqrt(2*pi));
plot(bins,f,'LineWidth',1.5)
set(gca,'FontName','Times New Roman','Fontsize',5)
xlabel("Droplet Diameter (μm)","FontSize",7);
ylabel("Frequency","FontSize",7);

