%% calculate the pressure, capillary pressure and flow rate in the centrifugal field
clear;clc;close
color = ["#FF8000";"#0080FF";"#33FF66";"#B34DCC"];

h = (20:20/300:40)/10^6; % microchannel height(m)
w = h; % microchannel width
l = 5/10^3; % microchannel length(m)
x = h/w; % the ratio of microchannel height to width
eta_w = 1.003/10^3; % viscosity of fluid（Pa*s）, here we used water as the fluid
y = 1-(192*x/power(pi,5)*(tanh(pi/2/x)+ ...
    tanh(3*pi/2/x)/power(3,5)+  ...
    tanh(5*pi/2/x)/power(5,5))); 
R = 12*eta_w*l./(h.^3.*w*y); % hydraulic resistance of a microchannel（Pa.s.m-3）
rpm = 2000:10:5000; % rotational speed（rpm）
[W,RPM] = meshgrid(w,rpm);
r1 = 17.5/1000;% the surface of fluid to the centrifugal center
r2 = 6.5/1000;% the height of fluid
omg = 2*pi*rpm/60; % angular velocity(rad)
rho_w = 998.2; % density of fluid（kg/m3）

P_w_first = 0.5*rho_w*(omg.^2)*((r1+r2)^2-r1^2); % pressure of fluid in the initial stage（Pa）
P_w_final = 0.5*rho_w*(omg.^2)*((r1+r2)^2-(r1+r2-l)^2); % pressure of fluid at the end stage (Pa)

figure(1) % the fluid pressure against rotational speed
plot(rpm,P_w_first,'Color',color(1))
hold on
plot(rpm,P_w_final,'Color',color(2))
plot(rpm,P_w_first - P_w_final,'Color',color(3))
area = fill([rpm fliplr(rpm)],[P_w_first fliplr(P_w_final)],'r','FaceColor',color(4),'edgealpha', '0', 'facealpha', '0.1');
set(gca,'FontName','Arial','FontSize',5)
ylabel("Pressure (Pa)",'FontName','Arial','FontSize',7)
xlabel("Rotational speed (rpm)",'FontName','Arial','FontSize',7)
legend("The beginning","The end","The decline")
pbaspect([1 0.7 1])
hold off

sita = 110; % contact angle of walls of microchannels
gama = 6.02/1000; % interfacial tension between the fluid and oil（N/m）
P_c = 2*cos(sita)*gama*(1./h + 1./w);% capillary pressure(Pa)
P_a = P_w_first + P_c; % fluid pressure + capillary pressure in the initial stage
P_a2 = P_w_final + P_c; % fluid pressure + capillary pressure at the end stage
speed = P_a./(R.*W.*W); % velocity in the microchannels in the initial stage
speed2 = P_a2./(R.*W.*W); % velocity in the microchannels at the end stage

figure(2) % the fluid velocity against rotational speed and microchannel diamensions
mesh(10^6*W,RPM,speed,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor','none')
hold on
mesh(10^6*W,RPM,speed2,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor','none')
mesh(10^6*W,RPM,speed-speed2,'FaceAlpha',0.5,'FaceColor','flat','EdgeColor','none')
test = speed-speed2;
set(gca,'FontName','Arial','FontSize',5)
cb = colorbar();
set(get(cb,'title'),'string','Velocity (m/s)','fontsize',5)
xl = xlabel("Microchannel Dimension (μm)",'FontName','Arial','FontSize',7);
yl = ylabel("Rotational speed (rpm)",'FontName','Arial','FontSize',7);
zlabel("Velocity (m/s)",'FontName','Arial','FontSize',7)
set(xl,'Rotation',17,'Position',[15 -1500 0.3]);   
set(yl,'Rotation',-25,'Position',[1 500 0.3]);   

%% keep w/h and a constant, get a series of V_in
r1 = 16.5/10^3;% the surface of fluid to the center
r2 = 7.5/10^3;% the height of fluid
w = 40/10^6;% the width of microchannal(um) (20/30/40)
h = 40/10^6;% the height of microchannel(um)

RPM = 2000;% rotate speed
omg = 2*pi*RPM/60;% rad

a_constant = omg^2*(r1+r2);% centrifugal acceleration


rho_w = 998.2; % the density of dispersed phase（kg/m3）
rho_o = 810.6; % the density of continuous phase（kg/m3）

gama = 6.02/1000;% interfacial tension of dispersed and continuous phase（N/m）
sita = 110;% contact angle

p_w = 0.5*rho_w*(omg^2)*((r1+r2)^2-r1^2);% the pressure of fluid
p_c = 2*cos(sita)*gama*(1/h + 1/w);% capillary pressure(Pa)
p_a = p_c + p_w;% final pressure

x = h/w; 
l = 5/1000; % length of microchannels（m）
eta_w = 1.003/1000; % viscosity of water（Pa*s）

y = 1-(192*x/power(pi,5)*(tanh(pi/2/x)+ ...
    tanh(3*pi/2/x)/power(3,5)+  ...
    tanh(5*pi/2/x)/power(5,5))); 
R = 12*eta_w*l/(power(h,3)*w*y); %resistance of microchannels（Pa.s.m-3）

v_min = p_a/(R*h*w);% minimal flow rate in the microchannels(V_in)


r2 = 27.5/1000; % highest liquid height
omg = (a_constant/(r1+r2))^0.5; 
p_w = 0.5*rho_w*(omg^2)*((r1+r2)^2-r1^2);%the pressure of fluid
p_a = p_c + p_w;% final pressure
v_max = p_a/(R*h*w);% maximum flow rate in the microchannels(V_in)


v = (v_min:(v_max-v_min)/4:v_max); % five flow rates between minimal and maximum flow rate
Ca = v*eta_w/gama; % capillary number
p_a = v*R*h*w;
p_w = p_a - p_c;
r2 = (p_w+(p_w.^2+(a_constant*rho_w*r1)^2).^0.5)/(a_constant*rho_w)-r1;
omg = (a_constant./(r1+r2)).^0.5;
RPM = 30.*omg/pi;
p1 = 0.5*rho_w.*(omg.^2).*((r1+r2-l).^2-r1^2);
array = [v',Ca',1000*r2',RPM'];