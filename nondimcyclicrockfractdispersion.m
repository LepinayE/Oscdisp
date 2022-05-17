%% NONDIM Cyclic Rock and Fracture model 2D with Dispersion

%{

y = b/2 * Y
t = b^2/(4*K) * \tau
u = U0* U 
x = (b^2 * U0)/(4*K)
and let \gamma =  b^2 * w / 4*K

Fracture model using forward FDM T(x,t)

Equation is T_t + Sin(\gammat) T_X = (K_disp /U0) T_XX - * df/dY(at Y=1)
which turns in to the FDM forward: when sin(\gammat)>0
T(i,k+1) = T(i,k ) - sin(\gamma*k*dt)*(dt/dX) (T(i,k)
-T(i-1,k))  - *(dt/dY)(f(j+1,k,i) - f(j,k,i)) +
(dt/dX^2)(K_disp/U0)(T(i+1,k) -2 T(i,k) + T(i-1,k))

when sin(wt)<0

T(i,k+1) = T(i,k ) - sin(\gamma*k*dt)*(dt/dX) (T(i+1,k)
-T(i,k))  - *(dt/dY)(f(j+1,k,i) - f(j,k,i)) +
(dt/dX^2)(K_disp/U0)(T(i+1,k) -2 T(i,k) + T(i-1,k))


% speed of fluid flow :
 u = U0sin(wt) so non dim u^ = sin(\gamma\tau)

Rock Model using forward FDM

Solving f_t = f_YY where f = f(y,t,x)
Using FDM s.t f(j, k+1,i) = f(j,k,i) + *dt/(dY)^2 (f(j+1,k,i) + f(i,j-1,k) -2 f(i,j,k) )

with 
IC
f(y,0,x) = 0
T(x,0) = 0

BCs
f(Y = 1,t,x) = T(x,t)

Injection: T(0,t) = 1 for 0< t< pi/(\gamma)
otherwise T_X (0,t) = 0 for pi/(\gamma)< t < 2 * pi/(\gamma)

Ends of channel held at initial temp
f(inf,t,X) = 0
f(Y,t,inf) = 0

T(inf,t) = 0


%}

close all;
clear;
clc;
clf;


% Set up 
K = 1*10^(-7); % Heat capacity Heat diffusivity could be up to -7 (molecular)
b = 1; %1 % height of fracture
K_disp = 3*10^(-5); % thermal dispersivity or could be much bigger
%K_disp = K_disp * 100;

U0 = 10^(-5); %1 ;% speed of fluid flow


%Vertical height of domain
LY = 10; %10 Length of rock cap
NY = 351; %10 Number of elements

Y_vec = linspace(0, LY, NY);  %Discretise Yspace starting from 0 and incrementing to Ly with length of vector Ny (LY/(NY -1))
dY = Y_vec(2) - Y_vec(1); %Concequence of linspace 

%Horizontal length of domain
LX = 10;  %length of fracture
NX = 301;  % number of elements %NX cannot be smaller than NY

X_vec = linspace(0, LX, NX);  %Discretise Xspace starting from 0 and incrementing to Lx with length of vector Nx
dX = X_vec(2) -X_vec(1);  


%Choosing time step so that d\tow < dX/(3), Courant condition for rock
%d\tow < ((dY)^2)/2 and Courant condition for fracture d\tow <(U0(dX)^2)/(2*K_disp)

if (dY^2)/ 4 <= dX/(3) && (dY^2)/ 4 <= (U0*(dX)^2)/(4*K_disp)

    dt = ((dY^2)/ 4) ;%Courant condition needs d\tow < ((dY)^2)/2

elseif (U0*(dX)^2)/(4*K_disp) <= (dY^2)/ 4 && (U0*(dX)^2)/(4*K_disp) <= dX/(3)

    dt  = (U0*(dX)^2)/(4*K_disp);
else
    dt = dX/(3);

end

spacet =1/dt; % = 1; % 1 year
t_vec = 0:dt:spacet *dt; %Discretise time
 

%Choose w such that pi/w is spacet *dt
ncyc = 1; % number of full cycles

gamma = 2*ncyc*pi/((length(t_vec) -1)*dt);  % (pi + 1)/ (20 *dt);
period  = ((length(t_vec) -1 )*dt) / ncyc; 

%Temperature matrix for rock f(y,t,x) and for fracture T(x,t)
f = zeros (length(Y_vec), length(t_vec), length(X_vec));
T = zeros (length(X_vec), length(t_vec));



%ICs and BCs

%{ 
Not needed since Matrix is of zeroes
% Fracture temp T, rock temp f(y,t,x)

% BC at end of domain held at temp = 0
T(end, :) = 0;
f(end, :,:) = 0;
f(:,:,end) = 0;

%Initial temp at t=0 is 0
T(:,1) = 0;
f(:, 1,:) = 0;

%}

% FDM Matrix & solver

%FDM Matrix for rock f

A0 = (1 - 2* dt/(dY)^2) * ones(length(Y_vec),1); 
Ap1 = dt/(dY)^2 *  ones(length(Y_vec)-1,1);
Am1 = dt/(dY)^2 *  ones(length(Y_vec)-1,1);

A = diag(A0) + diag(Ap1, 1) + diag (Am1, -1); %2D array

A(1,1) = 1;
A(end,end) = 1; % needed ?

A = repmat(A, 1, 1, length(X_vec)); % creates 3D array by repeatedly stacking A in the 3rd directions


% Determining f rock temperature and T fracture temperature sequentially
K_D = ((K_disp)/U0) * (dt/(dX)^2);

if K_D >= 1
   msg = 'Courant condition not satisfied for fracture';
     error(msg)
end

for k = 1 : length(t_vec) -1 
    
    if sin (gamma*k* dt) >= 0 

        %Injection BC at x = 0, y=0 for 0< t< pi/w
        lt = 1: length(t_vec); 
        T(1,lt) = 1; 
        
      
        %FDM Matrix for fracture T
        C0 = (1 - sin(gamma*k*dt)*(dt/dX) - 2*K_D ) * ones(length(X_vec),1) ;
        Cm1 = (sin(gamma*k*dt)*(dt/dX) + K_D) * ones(length(X_vec)-1,1);
        Cp1 = K_D* ones(length(X_vec)-1,1);

        C = diag(C0) + diag(Cm1,-1) + diag (Cp1,1);
        C(1 ,1) = 1 ;

    
    elseif sin (gamma *k *dt) < 0 

        %FDM Matrix for fracture T
        C0 = (1 + sin(gamma*k*dt)*(dt/dX) - 2*K_D ) * ones(length(X_vec),1);
        Cp1 = (-sin(gamma*k*dt)*(dt/dX) + K_D) * ones(length(X_vec)-1,1);
        Cm1 = K_D* ones(length(X_vec)-1,1);

        C = diag(C0) + diag(Cp1,1)+ diag(Cm1,-1);
    end

    

    %Fracture at time k 
    D = T(:, k);
    
    M1 = squeeze(f(2,k,:) - f(1,k,:)); % Heat flux term, heat loss to rock
    
    T2 = C*T(:,k) + (dt/dY)*(M1); % Find Temp at next time step

    
    T(2:length(X_vec) - 1, k+1) = T2(2:length(X_vec) - 1); 

       
    if sin(gamma*k*dt) <0 
    
        T(1, k+1) = T2(1);  
    
    end

   

    %Rock at time k 
     
    f2 = pagemtimes(A,f( :, k,:)) ;
    
    f( 2:length(Y_vec) - 1, k+1,:) = f2( 2:length(Y_vec) - 1, 1, :);

    % relate fracture temp to rock temp at boundary y = b/2
    
    f(1,k+1,:) = T(:,k+1);
   

    % relate fracture temp to rock temp at boundary y = b/2
    
    f(1,k,end) = T(end,k);

end

 
%Checking errors

fractsum = sum(T,2);
 if fractsum(end) > 0.0001
     msg = 'Fracture Temperature changes past domain size';
     error(msg)
 end

 
 rocksum = sum(f(:,:,2),2);
 if rocksum(end) > 0.0001
     msg = 'Rock Temperature changes past domain size';
     error(msg)
 end

%% Plotting

figure (1)

xpos = [5 31 61 91 181];
makel = (xpos -1 ) * dX;
makel = round(makel, 3);

subplot(4,1,1)
plot( X_vec,T(:,2),'*r', X_vec,T(:,floor(length(t_vec) /(ncyc *2))),'*b');
hold on 
plot(X_vec, T(:, floor( ((ncyc *2) -1)* length(t_vec)/(ncyc *2) )),'*k', X_vec,T(:,length(t_vec)),'*c');
hold on 

title('Heat Diffusion in Fracture T(x,t) at different times')
xlabel('x')
ylabel('Fracture Temperature')
legend ('Immediately after Injection', 'End of first Injection','End of Eighth Injection', 'End')


subplot(4,1,2)
plot( f(:,floor(length(t_vec) /(ncyc *2)),xpos(1)),Y_vec,'*r',f(:,floor(length(t_vec) /(ncyc *2)),xpos(2)),Y_vec,'*b',f(:,floor(length(t_vec) /(ncyc *2)),xpos(3)),Y_vec,'*k');
hold on 
plot(f(:,floor(length(t_vec) /(ncyc *2)),xpos(4)),Y_vec,'*c',f(:,floor(length(t_vec) /(ncyc *2)),xpos(5)),Y_vec,'*g'); %'linewidth',2);
hold on

title('Heat Diffusion in Rock at End of first Injection')
xlabel('Rock Temperature f(y,t,x)')
ylabel('y')
legend(sprintfc('x = %.3f', makel, false))


subplot(4,1,3)
plot( f(:,floor( ( (ncyc *2) -1)*length(t_vec) /(ncyc *2) ),xpos(1)),Y_vec,'*r',f(:,floor(((ncyc *2) -1)*length(t_vec)/(ncyc *2)),xpos(2)),Y_vec,'*b',f(:,floor(((ncyc *2) -1)*length(t_vec) /(ncyc *2)),xpos(3)),Y_vec,'*k');
hold on 
plot(f(:,floor(((ncyc *2) -1)*length(t_vec) /(ncyc *2)),xpos(4)),Y_vec,'*c',f(:,floor(((ncyc *2) -1)*length(t_vec) /(ncyc *2)),xpos(5)),Y_vec,'*g'); %'linewidth',2);
hold on

title('Heat Diffusion in Rock at End of Last Injection')
xlabel('Rock Temperature f(y,t,x)')
ylabel('y')
legend(sprintfc('x = %.3f', makel, false))

subplot(4,1,4)
plot( f(:,(length(t_vec) ),xpos(1)),Y_vec,'*r',f(:,(length(t_vec) ),xpos(2)),Y_vec,'*b',f(:,(length(t_vec)),xpos(3)),Y_vec,'*k');
hold on 
plot(f(:,(length(t_vec) ),xpos(4)),Y_vec,'*c',f(:,(length(t_vec) ),xpos(5)),Y_vec,'*g'); %'linewidth',2);
hold on

title('Heat Diffusion in Rock at End of time ')
xlabel('Rock Temperature f(y,t,x)')
ylabel('y')
legend(sprintfc('x = %.3f', makel, false))


%% Temp split between injection and extraction and per cycles
nlines = 6; %number of lines wanted on graph

kindex = 2 : floor((length(t_vec)-1) /(2*nlines)): length(t_vec)-1;
kindex2 = kindex *dt;

id = sin(gamma*kindex*dt)>=0 ;

%Fracture temp
%Injection
figure(2)

Tinj = T(:,kindex(id));

plot(X_vec, Tinj ,'linewidth',2)
drawnow
title('Heat in Fracture for different times during Injection')
legend(sprintfc('Time = %.2f yr', kindex2(id), false))
xlabel('x')
ylabel('Fracture Temperarture T(x,t)')
    
%Extraction 
figure(3)

Text = T(:,kindex(~id));

plot(Text,'linewidth',2)
drawnow
title('Heat in Fracture for different times during Extraction')
legend(sprintfc('Time = %.2f yr', kindex2(~id), false))
xlabel('x')
ylabel('Fracture Temperarture T(x,t)')


%Rock temp 
%Injection
figure(4)

plot(Y_vec,f(:,kindex(id),1),'linewidth',2);

drawnow
title('Heat in Rock for different times during Injection at x = 1')
xlabel('y')
ylabel('Rock Temperature')
legend(sprintfc('Time = %.2f yr', kindex2(id), false))

%Extraction
figure(5)

plot(Y_vec,f(:,kindex(~id),1),'linewidth',2);

title('Heat in Rock for different times during Extraction at x = 1')
xlabel('y')
ylabel('Rock Temperature')
legend(sprintfc('Time = %.2f yr', kindex2(~id), false))
%% Heat Ratio

lt = 1:length(t_vec);
id2 = sin(gamma*lt*dt)>=0 ; % id2 will be true 1 when we are injecting
lt2 = lt;


firstint = find(~id2,1) -1 ; % Last index value in lt2 when sin >0  during first cycle
lt2 = mod(lt2,firstint);

lt3 = find(~lt2); % finds the index position of the 0s in lt2 which are the end of the injection cycles

lt2(lt3) = firstint;

lt2(~id2) = 0; %lt2 = lt value when injecting and 0 otherwise

    
F1 = shiftdim(f,2); % permutes f s.t. F1 = f(x,y,t)
sumrock = squeeze(sum(F1,[1 2])); % sums all vales of F1 on each page -> Heat in Rock at each time step


sumfract = sum(T); %Heat in fracture at each time step in vector form

Heataq = sumfract(lt) *dX  + (sumrock(lt).') *dX *dY; % Total Heat in Aquifer at time step t

Heatin = sin(gamma*dt *lt2) .* lt2 *dt; % Heat only going in during injections Tempinj x Sin(gamma * dt) x dt and restarts at each injection
Heatin = cumsum(Heatin);

Heatratio = Heataq./Heatin;% Ratio of total heat in Aquifer vs heat injected in system (at each time step)
Heatratio(~id2) = 0;

R = Heataq./Heatin;
% R = diff(Heataq)./diff(Heatin); %Ratio of difference of heat in aquifer between each time step vs heat injected in system at each time step

% ratio plots
figure(6)
plot(Heatratio)
xlabel('Time step')
ylabel('Heat ratio')
title('Ratio of Heat in Aquifer to Heat injected in system during cycle at time step t')
xlim([5 length(t_vec)])


%% Conservation of Heat
lt = 1:length(t_vec);
id2 = sin(gamma*lt*dt)>=0 ; % id2 will be true 1 when we are injecting
lt2 = lt;


firstint = find(~id2,1) -1 ; % Last index value in lt2 when sin >0  during first cycle
lt2 = mod(lt2,firstint);

lt3 = find(~lt2); % finds the index position of the 0s in lt2 which are the end of the injection cycles

lt2(lt3) = firstint;

lt2(~id2) = 0; %lt2 = lt value when injecting and 0 otherwise

    
F1 = shiftdim(f,2); % permutes f s.t. F1 = f(x,y,t)
sumrock = squeeze(sum(F1,[1 2])); % sums all vales of F1 on each page -> Heat in Rock at each time step


sumfract = sum(T); %Heat in fracture at each time step in vector form

Heataq = sumfract(lt) *dX  + (sumrock(lt).') *dX *dY; % Total Heat in Aquifer at time step t

Heatin = sin(gamma*dt *lt2) .* lt2 *dt; % Heat only going in during injections Tempinj x Sin(gamma * dt) x dt and restarts at each injection

Heatratio = Heataq./Heatin;% Ratio of total heat in Aquifer vs heat injected in system (at each time step)
Heatratio(~id2) = 0;


R = diff(Heataq)./diff(Heatin); %Ratio of difference of heat in aquifer between each time step vs heat injected in system at each time step

% ratio plots
figure(6)
plot(Heatratio)
xlabel('Time step')
ylabel('Heat ratio')
title('Ratio of Heat in Aquifer to Heat injected in system during cycle at time step t')
xlim([5 length(t_vec)])

figure(7)
plot(R)
xlabel('Time step')
ylabel('Heat different ratio of each time step')
title(['Ratio of the difference of Heat in Aquifer between time t and time ' ...
    't-1 to the difference of Heat injected in system between time t and time t-1']) 


%Recovery temp
figure(8)
for n = 1:ncyc

    txt2 = ['Cycle number = ',num2str(n)];
        
    extrT1 = [];
    kmin = round(((2*n-1)*pi)/(gamma * dt));
    kmax1 = floor((n*2*pi)/(gamma * dt));

    extrT1 = T(2,kmin:kmax1);
   
    plot(extrT1, 'linewidth',2, 'DisplayName', txt2);
    hold on

    legend('Location','northeastoutside') 
    title('Extraction temperature of fracture water at injection point T(0,t)')
    xlabel('Time')
    ylabel('Extraction Temperarture ')

end

%Heat in Aquifer
figure(18)

Hhat = Heataq./(sqrt(lt)) ; 

plot((t_vec).^(1/2),Hhat)
title ('Total heat in system rescaled H/surd t','Interpreter','tex')
xlabel('\surd t','Interpreter','tex')
ylabel('H/\surd t','Interpreter','tex')


figure(19)

n = floor ( length(t_vec) /ncyc);

avgHeataq =movmean(Heataq,n);

plot(avgHeataq)
title('Total heat in system averaged over each cycle')
ylabel('Average total heat in system')
xlabel('time')

figure(25) %log-log plot of moving mean of total heat in Aquifer
tl2 =(n/2):(length(t_vec) - (n/2));
tl2 = log(tl2);
plot(tl2, avgHeataq( (n/2):(length(t_vec) - (n/2)) ))

xlabel('log(time)')
ylabel('log(averaged Total heat in Aquifer)')
title('Log-Log graph of Total Heat in Aquifer averaged over each cycle between t = 2000 and 13000')


%% Contour Plots of Fracture temp 
figure(9)

contour(t_vec,X_vec,T)
xlabel('Time')
ylabel('x-direction')
title('Contour of Fracture temperature')
ylim ([0,40])
colorbar

% Contour Plots of Rock temp at y = 10
figure(10)

faty = squeeze(f(10,:,:));
faty = shiftdim(faty,1);

contour(faty)
xlabel('Time')
ylabel('x-direction')
title('Contour of Rock temperature at y = 10')
ylim ([0,150])
colorbar

% Rescaled x/t^(1/2) Contour plot for fracture temp 
figure(13)
[t_g, x_g] = meshgrid(t_vec, X_vec);

for i = 1:length(t_vec);

    x_g2(:,i) = x_g(:,i) / sqrt(t_vec(i));

end
contour(t_g, x_g2, T)
ylim([0 2])
ylabel('x/\surd t','Interpreter','tex');
xlabel('time t')
title('Contour plot of Fracture Temperature T(x/\surd t , t)','Interpreter','tex');
colorbar

% Rescaled x/t^(1/2) Contour Plots of Rock temp at y = 10

figure(14)

contour(t_g, x_g2, faty)
ylim([0 2])
ylabel('x/\surd t','Interpreter','tex');
xlabel('time t')
title('Contour plot of Rock Temperature f(y = 10,t, x/\surd t)','Interpreter','tex');
colorbar


%%
% Changes in Fracture temp through the cycles

figure(11)
plot (t_vec, T(2,:),'r','LineWidth', 2)
hold on 
plot(t_vec, T(15,:),'b','LineWidth',2)
ylabel('Fracture temperature')
xlabel('time')
title ('Fracture temperature variations through the cycles at different positions')
legend ('x=2', 'x=15')

% Changes in Rock temp through the cycles at y =10

figure(12)
plot (t_vec, faty(2,:),'r','LineWidth', 2)
hold on
plot(t_vec, faty(15,:),'b','LineWidth', 2)
ylabel('Rock temperature')
xlabel('time')
title ('Rock temperature variations through the cycles at different positions at y=10')
legend ('x=2', 'x=15')



%% Center of Mass and mean temperature

% Fracture temp Center of mass
figure(15)

lx = 1:length(X_vec);

TiXi = lx * T; % Sum of TiXi for each time t

fracTbar = TiXi ./ sum(T); % Center of mass of fracture temp
plot(fracTbar)
xlabel('time')
ylabel('Center of Mass of Temperature in Fracture')
title('Center of Mass of Fracture Temperature through the cycles')

%Fracture temp mean temp
figure(16)

plot(mean(T))
xlabel('time')
ylabel('Mean Fracture Temperature')
title('Mean Fracture Temperature through the cycles')


figure(17)

n = floor ( length(t_vec) /ncyc);

avgfracTbar = movmean(fracTbar,n);
plot(avgfracTbar)

title('Average of the center of mass of the fracture temperature throught the cycles')
xlabel('time')
ylabel('Average of Center of mass of fracture temperature')

%log - log plot
%{
figure(20)
x1 = 1 :  length(t_vec);
x2 = log (x1);
y1 = log(avgfracTbar);

plot(x2,y1)
title('Log-Log plot of Average of the center of mass of the fracture temperature throught the cycles')
xlabel('Log of time')
ylabel('Log of Average of fracture temperature')

hold on 
x3 = 5:15;


y3 = 0.25 * x3 -0.15;
 plot (y3)
%}

%Rock temp, center of mass in y direction

figure(21)

ly = 1:length(Y_vec);
fiYi = squeeze(pagemtimes(ly, f)) ; %sum of fiyi

rockfbar = fiYi ./squeeze (sum(f,1)); %COM of rock temp
rockfbar = rockfbar .';

n = floor ( length(t_vec) /ncyc);
t1 = 1: n: length(t_vec);

plot(X_vec, rockfbar(:,t1)) 
drawnow
legend(sprintfc('Time = %d', t1, false))
xlabel('x')
ylabel('Center of Temperature in Rock ')
title('Horizontal Center of mass of rock temperature through the cycles')




%Rock temp, mean in y direction
figure (22)

meanfract = squeeze(mean(f)).';
t1 = 1: n: length(t_vec);

plot(X_vec, meanfract(:,t1))
drawnow
legend(sprintfc('Time = %d', t1, false))
xlabel('x')
ylabel('Horizontal mean Temperature in Rock ')
title('Mean of rock temperature through the cycles')

%Center of mass of rock temperature in 1st cycle
figure(23)

for t = 1: 100: length(t_vec)/ncyc ; 
    
    txt3 =  ['Time = ',num2str(t)];
    plot(X_vec, meanfract(:,t),'DisplayName', txt3) 
    hold on 
    
    xlim([0 35])
    ylim([0,0.2])
    drawnow 
  
    xlabel('x')
    ylabel('Horizontal Center of Temperature in Rock')
    ylim([0 0.05])
    title('Center of mass of rock temperature in 1st cycle')
    legend show
end
%%
%video of temperature changing in rock through the cycles
figure(24)
for i = 1:150

imshow(squeeze(f(:,100 * i,:)).')

colormap jet
caxis([0 1])
pause(0.1)

xlabel('x')
ylabel('y')
title('Temperature change in rock during cycles')
end

