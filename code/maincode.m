% experimental data
% collected from a vertical-take-off-landing (VTOL) aerial device test
% identification procedure for quasi-LPV of linear systems
% updated April 04, 2023  (13h)

clear all, close all, clc, format long, format compact,

disp(' .... computing code for robust linear system from VTOL (it may take some minutes) ...')


fid = fopen('listaQuasiLPV.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline)
    nome{count} = sprintf('%s',tline);
    tline = fgetl(fid);
    count = count+1;
end

fclose(fid);


%load Data_Identification_08-Mar-2023_14h00m.mat
%load Data_Identification_09-Mar-2023_14h34m.mat
%load Data_Identification_14-Feb-2023_14h24m.mat
%load Data_Identification_14-Mar-2023_12h20m.mat
%load Data_Identification_16-Mar-2023_16h51m.mat
%load Data_Identification_16-Mar-2023_17h14m.mat
%load Data_Identification_16-Mar-2023_17h29m.mat


for cx=1:max(size(nome))
    
    load(nome{cx});
    
    duty = ans.inputU(100:end);
    y1 = ans.angle(100:end);
    y2 = ans.velocity(100:end);
    
    x1real = y1';
    x2real = y2';
    %u = duty/10;  %procedure to set '0<u<1'
    
    %dn = duty;  % duty cycle normalized
    %u =10*(-0.0001402372 + 0.0009364425*dn + 0.0003874307*(dn.^2));
    u = -0.001051853 + 0.007023813*duty + 0.002905935*(duty.^2);
    
    
    Ts=TS;
    
    
    VarEPS = 1e-6;
    
    Nit = max(size(duty));
    Ts=TS;
    t=[0:Ts:Ts*Nit];
    
    P{1} = diag([VarEPS VarEPS 10^3 10^3  VarEPS 10^3]);
    Delta{1} = [1 ; 10*TS ; 0; 0; 0.3; 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % procedure of Kalman filter to identify
    % the parameters of the linear system
    % x(k+1) = A*x(k) + B*u(k) + cte
    % We write the system above in the form of
    % \theta(k+1) = \theta(k) + w(k)
    %  x(k) = m(k)*theta(k) + e(k)
    % The last two equations feed the KF
    for k=2:Nit
        m{k-1} = [x1real(k-1) x2real(k-1)  0   0  0  0;
            0  0  x1real(k-1) x2real(k-1) u(k-1) 1];
        K{k} = P{k-1}*m{k-1}'*inv( eye(2,2) + m{k-1}*P{k-1}*m{k-1}' ) ;
        Delta{k} = Delta{k-1} + K{k}*( [x1real(k) x2real(k)]' - m{k-1}*Delta{k-1} );
        P{k} = diag([VarEPS VarEPS 10^3 10^3  VarEPS 10^3])  + P{k-1} - K{k}*m{k-1}*P{k-1};
        P{k} = (P{k}+P{k}')/2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    heap_a11=[]; heap_a12=[];
    heap_a21=[]; heap_a22=[];
    heap_b21=[]; heap_b22=[];
    vecx1=[];vecx2=[];
    xsim{1}= [x1real(1) x2real(1)]';
    for k=1:Nit
        vecx1 = [vecx1 xsim{k}(1)];
        vecx2 = [vecx2 xsim{k}(2)];
        A = [Delta{k}(1:2,:)';
            Delta{k}(3:4,:)'];
        B = [0 0; Delta{k}(5:6,:)'];
        xsim{k+1} = A*xsim{k} + B*[u(k); 1];
        
        heap_a11 = [heap_a11  A(1,1)];
        heap_a12 = [heap_a12  A(1,2)];
        heap_a21 = [heap_a21  A(2,1)];
        heap_a22 = [heap_a22  A(2,2)];
        heap_b21 = [heap_b21 B(2,1)];
        heap_b22 = [heap_b22 B(2,2)];
    end
    
    
    for k=1:Nit
        A_d{k} = [heap_a11(k) heap_a12(k);
            heap_a21(k) heap_a22(k)];
        B_d{k} = [0; heap_b21(k)];
    end
    
    
    error = norm(x1real - vecx1);
    
    cx
    if (error>10)
        disp(nome{cx});
        %error
    end
end


figure(1)
subplot(2,1,1)
hold on
plot(Ts*[1:max(size(x1real))],x1real,'r-.','LineWidth',2);
plot(Ts*[1:max(size(x1real))],vecx1,'k','LineWidth',1)
hold off
grid
legend('real','estimated');
xlabel('seconds'),ylabel('x1')

subplot(2,1,2)
hold on
plot(Ts*[1:max(size(x2real))],x2real,'r-.','LineWidth',2);
plot(Ts*[1:max(size(x2real))],vecx2,'k','LineWidth',1)
hold off
grid
legend('real','estimated');
xlabel('seconds'),ylabel('x2')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)

subplot(4,2,1)
plot(Ts*[1:max(size( u ))], u ,'b','LineWidth',2);
ylabel('u'),
grid

subplot(4,2,2)
plot(Ts*[1:max(size( heap_a11 ))],heap_a11,'k','LineWidth',2);
grid, ylabel('a11'),

subplot(4,2,3)
plot(Ts*[1:max(size( heap_a12 ))],heap_a12,'k','LineWidth',2);
grid, ylabel('a12'),

subplot(4,2,4)
plot(Ts*[1:max(size( heap_a21 ))],heap_a21,'k','LineWidth',2);
grid, ylabel('a21'),

subplot(4,2,5)
plot(Ts*[1:max(size( heap_a22 ))],heap_a22,'k','LineWidth',2);
grid, ylabel('a22'),

subplot(4,2,6)
plot(Ts*[1:max(size( heap_b21 ))],heap_b21,'k','LineWidth',2);
grid, ylabel('b12'),

subplot(4,2,7)
plot(Ts*[1:max(size( heap_b22 ))],heap_b22,'k','LineWidth',2);
grid, ylabel('b22');


return

% 
% vecx=[-4:0.0001:4];
% vecy=[]; vecyH=[];
% k=2;
% for i=1:max(size(vecx))
%     vecy = [vecy sign(vecx(i))];
%     vecyH = [vecyH tanh(k*vecx(i))];
% end
% 
% figure(101)
% hold on
% plot(vecx,vecy,'r','LineWidth',4)
% plot(vecx,vecyH,'b','LineWidth',2)
% hold off

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % procedure to let the matrices belong to a polytope
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% points =[];
% for k=1:max(size(heap_a11))
%     vector{k} = [heap_a11(k) heap_a12(k) heap_a21(k) heap_a22(k)]; %heap_b1(k) heap_b2(k)];
%     points = [points; vector{k}];
% end
% indices = convhulln(points);
% vertices = unique(vec(indices));
% 
% 
% count =1;
% Bvec{1}=[heap_a11(1) heap_a12(1) heap_a21(1) heap_a22(1)];
% for k=1:max(size(heap_a11))
%     vector{k} = [heap_a11(k) heap_a12(k) heap_a21(k) heap_a22(k)]; %heap_b1(k) heap_b2(k)];
%     if (norm(Bvec{count}-vector{k})>0.001)
%         count = count+1;
%         Bvec{count} = vector{k};
%     end
% end
% 
% 
% points2 =[];
% for k=1:max(size(Bvec)) %max(size(heap_a11))
%     points2 = [points2; Bvec{k}];
% end
% indices2 = convhulln(points2);
% vertices2 = unique(vec(indices2));
% 
% 
% return
% 

u=[0:0.001:10];

figure(41)
vecF=[];
for k=1:max(size(u))
    vecF = [vecF  -0.001051853 + 0.007023813*u(k) + 0.002905935*u(k)^2];
end
plot(u,vecF)


