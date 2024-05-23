%% theoretical deduction
clc;clear;close;
filepath = "theroy.xlsx";
A = readmatrix(filepath,"Sheet","theroy"); % import the table
ac = []; % centrifugal acceleration
color = ["#FF8000";"#0080FF";"#33FF66";"#B34DCC"];% Colors are used to differentiate between different centrifugal accelerations

for i = 1:1:4 % four centrifugal accelerations
     ac = [ac;A((i-1)*7+1,2)];
end
ac = round(ac,2);

field = 4;% Column 4 is the gravity field and column 6 is the experimental data

channel = [20;30;40]; % microchannel diamension
gama = 6.02/1000; % interfacial tension（N/m）
figure(1) % resistance factor plotted against diameter
hold on
box on
d_all = []; % all diameters
rs_all = []; % all resistance factors
combinedlabel = [];
for i = 1:1:3 % three sizes of microchannels
    w = (10 + 10*i)/10^6;% microchannel diamensions
    h = w;
    v = A(59:85,8*(i-1)+1);% import V_in,normal_1:27
    % vacant rectangular microchannel 30:56,vacant circular microchannel 59:85
    d = A(59:85,8*(i-1)+field)/10^6;
    r = d/2;% radius
    rs = 2*gama.*(1/w+1/h-1./r)./v;% resistance factor
    if i == 1
        pattern = 'o';
    elseif i == 2
        pattern = 'square';
    else
        pattern = 'diamond';
    end
    for j = 1:1:4 % four centrifugal accelerations
        dotplot = scatter(d(7*(j-1)+2:7*(j-1)+6)*10^6,rs(7*(j-1)+2:7*(j-1)+6),12 ...
           ,'MarkerFaceColor','none','MarkerEdgeColor',color(j),'Marker',pattern);
        label = "w = h = " + num2str(w*10^6) + "μm, " + ...
            "a = " + num2str(ac(j))+ " m/s^2";
        combinedlabel = [combinedlabel;label'];
    end
    d_all = [d_all;d];
    rs_all = [rs_all;rs];
end

line = [];
a = [];% coefficient α
b = [];% coefficient β
rsquare = [];
for i = 1:1:4 % four centrifugal accelerations
    d = []; % diameter
    rs = []; % resistance factor
    for j = 1:1:3
        % combine the data of three microchannels
        d = [d;d_all(27*(j-1)+7*(i-1)+2:27*(j-1)+7*(i-1)+6)];
        rs = [rs;rs_all(27*(j-1)+7*(i-1)+2:27*(j-1)+7*(i-1)+6)];
    end
    d = 10^6*d;
    d(isnan(d)==1)=[];
    rs(isnan(rs)==1)=[];

    % exponential fit (μm)
    if field == 4 % simulation data
        lower = [130000,-0.05];% normal'Lower',[130000,-0.05],'Upper',[600000,-0.01]
        upper = [1000000,-0.01];%vacant rectangular microchannel [1000000,-0.01];vacant circular microchannel[1000000,-0.01]
    elseif field == 6 % experimental data
        lower = [130000,-0.02];%'Lower',[130000,-0.02],'Upper',[180000,-0.01]
        upper = [180000,-0.01];
    end
    
    expFitType = fittype('a * exp(b * x)','independent','x','coefficients',{'a','b'});
    [cfun,rsquare1] = fit(d, rs, expFitType,'Lower',lower,'Upper',upper, ...
        'StartPoint',[mean(d),mean(rs)]); 
    
    % Obtain the fitting parameters for four centrifugal accelerations
    a = [a;cfun.a]; 
    b = [b;cfun.b];
    rsquare = [rsquare;rsquare1.rsquare];

    xi = 90:0.1:max(d);
    yi = cfun(xi);
    line = [line;plot(xi, yi,'--','Color',color(i))];
    
end
hold off  
a = round(a,1);
b = round(b,3);
rsquare = round(rsquare,2);
label = "a = " + num2str(roundn(ac,-1)) + "m/s^2," + ...
    "R_f = " + num2str(a) + "e"+num2str(b)+"D" + ...
    ",R^2 = " + num2str(rsquare);
combinedlabel = [combinedlabel;label];
legend(combinedlabel)
xlabel("Diameter (μm)",'FontName','Arial',"FontSize",7)
ylabel("Resistance factor (Pa·s/m)",'FontName','Arial',"FontSize",7)

pbaspect([1 0.7 1])
% xlim([90 250])
set(gcf,'PaperUnits','centimeters')
set(gca,'FontName','Arial',"FontSize",5)

figure(2) % predictive results
box on
hold on
p = [];
dot = [];
combinedlabel = [];
for i = 1:1:4
    for j =1:1:3
        syms d_pred % the independent variable to be predicted(diameter)
        w = (10+10*j)/10^6;%通道宽度、高度(m)
        h = w; 
        v_fun(d_pred) = 2*gama.*(1/w+1/h-2./(d_pred/10^6))/(a(i).*exp((b(i)).*d_pred));
        d_want = linspace(50,450,100);  % predicted diameter range
        x = A(7*(i-1)+2:7*(i-1)+6,8*(j-1)+1);% import V_in
        y = A(7*(i-1)+2:7*(i-1)+6,8*(j-1)+field);% diameter
        if j == 1 % three microchannels
            pattern = 'o';
            group = '20μm';
            linestyle = ':';
        elseif j == 2
            pattern = 'square';
            group = '30μm';
            linestyle = '--';
        else
            pattern = 'diamond';
            group = '40μm';
            linestyle = '-.';
        end
        dot(i,j) = scatter(x,y,12,pattern,'MarkerFaceColor','none', 'MarkerEdgeColor',color(i),'DisplayName',group);
        p(i,j) = plot(v_fun(d_want),d_want,'LineStyle',linestyle,'Color',color(i));
        p1(i,j) = plot(v_fun(d_want),d_want,'LineStyle',linestyle,'Color',color(i));
       
    end
end
combinedLabel1 = num2str(channel) + " μm";
combinedLabel2 = "a = " + num2str(ac) + " m/s^2";
combinedLabel3 = num2str(channel) + " μm-predict";
combinedLabel = [combinedLabel1',combinedLabel2',combinedLabel3'];
lgd = legend([dot(1,:)';p(:,3);p1(1,:)'], combinedLabel);

hold off
xlabel("Velocity (m/s)",'FontName','Arial','FontSize',7)
ylabel("Diameter (μm)",'FontName','Arial','FontSize',7)
xlim([0,0.6*(field/2-1)])
pbaspect([1 0.7 1])
set(gca,'FontName','Arial',"FontSize",5)


%% Heat mapping of predictions
clc
k = 1000;% intervals
v_vecter = linspace(0,0.5,k); %the range of V_in
d_matrix = zeros(4,k);

for i = 1:1:4
    w = 40/10^6;
    h = w;
    v_fun2(d_pred) = 2*gama.*(1/w+1/h-2./(d_pred/10^6))/(a(i).*exp((b(i)).*d_pred));
    d_want = linspace(50,1000,10000);  % range of predicted diameters 
    v_want = double(v_fun2(d_want));
    for j = 1:1:length(v_vecter) 
        number = next2(v_want,v_vecter(j));% index of the closest of the two
        d_matrix(i,j) = d_want(number);
        f_matrix(i,j) = v_vecter(j)*w*h/(4/3*pi*(d_want(number)/10^6/2)^3);
    end

end

figure(1) % droplet daimeter prediction
h1 = heatmap(d_matrix,"GridVisible","off","CellLabelColor","none");
h1.FontName = "Arial";
h1.FontSize = 18;

figure(2) % generation frequency prediction
h2 = heatmap(f_matrix,"GridVisible","off","CellLabelColor","none");
h2.FontName = "Arial";
h2.FontSize = 18;

function output = next2(v,v_exp)
%Find the index that is closest to the desired V_in in the predicted V_in
    mine = abs(v(1)-v_exp);
    output = 1;
    for i = 1:1:length(v)
        if  mine > abs(v(i)-v_exp)
            mine = abs(v(i)-v_exp);
            output = i;
        end
    end
end