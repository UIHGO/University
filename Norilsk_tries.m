function [] = Norilsk_tries
  clear all
  global nl
  global m mm1 mp1 mp2  %m - число КО на [0,zL]   
  global ibc_i0 ibc_mp1
  global nt
  global iorder
  global rho_d c_d W_tot lambda_M lambda_T u_phase
  global eps D_e
  global t0 tk t dt it_per li
  global zL 
  global lambda c_rho
  global ae ae_ip ae_im be
  global n_lr    %n_lr - номер литологического слоя КО
  
  global he theta
  global z u
  global zw hz
  global zw_z
  global f_Ci0  f_Pi0
  global f_Cmp1 f_Pmp1
  global S_C S_P
  global alpha_w beta_w gamma_w
  global rkappa c_w c_ice q_E
  
% Начало программы main
  Data_preparation();
  for it=1:nt
      t=t0+it*dt;
      Material_properties_and_boundary_conditions();
      SLAU();
  end
  OUTPUT();
% !end program main

%///////////////////////////
function [] = Data_preparation()
  eps=1.e-20;   %eps=0.
  nl=3;
  li=12;
%Температура окружающей среды:
%Норильск
 %theta(1)=-26.8; theta(2)=-26.03; theta(3)=-20.81; %зима
 %theta(4)=-13.07; theta(5)=-4.26; theta(6)=7.25; %весна
 %theta(7)=14.45; theta(8)=11.01; theta(9)=4.1; %лето
 %theta(10)=-8.13; theta(11)=-20.96; theta(12)=-24.59; %осень
%Прошлое
 theta(1)=-28.3; theta(2)=-28; theta(3)=-22; %зима
 theta(4)=-13.8; theta(5)=-5.8; theta(6)=6.2; %весна
 theta(7)=14.9; theta(8)=10.07; theta(9)=3.8; %лето
 theta(10)=-9.3; theta(11)=-21.7; theta(12)=-25.5; %осень
%RCP2.6
 %theta(1)=-23.4; theta(2)=-22.88; theta(3)=-17.91; %зима
 %theta(4)=-10.42; theta(5)=-1.86; theta(6)=9.4; %весна
 %theta(7)=16.35; theta(8)=13.16; theta(9)=6.5; %лето
 %theta(10)=-5.48; theta(11)=-18.06; theta(12)=-21.44; %осень
%RCP8.5
 %theta(1)=-17.7; theta(2)=-17.4967; theta(3)=-12.8433; %зима
 %theta(4)=-5.67; theta(5)=2.57333; theta(6)=13.5167; %весна
 %theta(7)=20.15; theta(8)=17.2767; theta(9)=10.9333; %лето
 %theta(10)=-0.73; theta(11)=-12.9933; theta(12)=-16.0567; %осень
%Коэффициент теплообмена h:
 he(1)=1.032449; he(2)=0.832335;  he(3)=0.729155;
 he(4)=0.680178; he(5)=0.884148;  he(6)= 12.6000;
 he(7)=13.4;   he(8)=14.0000;  he(9)= 14.3;
 he(10)=2.041904; he(11)=1.061922;  he(12)= 0.939475;

  iorder=2;
  rkappa=79.4;
  c_w=1.006;
  c_ice=0.49;
  q_E=0.0;   % q_E=0.043;
      
  it=0;
  t0=0.0;
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  tk=438000; 
  nt=60000;
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  dt=(tk-t0)/nt; %dt=7.3 часа.
  it_per=100; %it_per*dt=100*7.3=730 - 100 шагов = 1 месяц
  
  zL=20.0;
  m=200;
  mp2=m+2;
  mp1=m+1;
  mm1=m-1;
     
  ibc_i0=2; %1;
  ibc_mp1=2; %1;

  t=t0; 
  f_Ci0=he(1)*theta(1);  f_Pi0=-he(1);  
  f_Cmp1=q_E; f_Pmp1=0.0;

%Построение сетки:
  GRID();
%Вычисление номера литологической разности КО (i):      
  number_lithologic();
  
%Ввод начальных условий:
    %Универсальные (средние за год на основе многолетних измерений)
  for i=1:2    u(i)=-2.5;      end
  for i=3:11    u(i)=-1.6;      end
  for i=12:129  u(i)=-2.1;      end
  for i=130:mp2   u(i)=-2.8;      end
    %Январь
  %for i=1:2    u(i)=-13;      end
  %for i=3:11    u(i)=-4.5;      end
  %for i=12:129    u(i)=-1.5;      end
  %for i=130:mp2  u(i)=-2.8;      end
    %Июль
  %for i=1:2   u(i)=14.8;      end
  %for i=3:11    u(i)=4;      end 
  %for i=12:129    u(i)=-2;      end
  %for i=130:mp2    u(i)=-2.8;      end
  
%Объемные источники тепла отсутствуют:  
  for i=1:m
      S_C(i)=0.0;      
      S_P(i)=0.0;
  end     

%Свойства грунтов:
  rho_d(1)=1540.0; %плотность
  c_d(1)=0.23;     %теплоемкость
  W_tot(1)=0.25;   %влажность
  lambda_M(1)=1.5; %теплопроводность мерзлого
  lambda_T(1)=1.3;%теплопроводность талого
  u_phase(1)=0.0;  %температура фазового перехода

  rho_d(2)=1500.0;
  c_d(2)=0.23;
  W_tot(2)=0.2;
  lambda_M(2)=1.13;
  lambda_T(2)=0.95;
  u_phase(2)=0.0;

  rho_d(3)=2900.0;
  c_d(3)=0.2;
  W_tot(3)=0.0;
  lambda_M(3)=1.12;
  lambda_T(3)=1.12;
  u_phase(3)=0.0;
  
% Вычисление коэффициентов функции Ww=Ww(u):  
  W_w_coef();
  
%  fprintf(['Press any key to continue...\n'])
%  pause
end % Data_preparation

%//////////////////
function [] = GRID()
%координаты границ литологических разрезов по оси z;
  hz_c=zL/m;
  zw_z(1)=0.0;
  zw_z(2)=90*hz_c; 
  zw_z(3)=130*hz_c;
  zw_z(4)=zL;

%построение сетки - zw(i) (координаты граней КО) во всей области    
% zw(1)=0.0;
  for i=1:m
      hz(i)=zL/m; %размеры КО
      zw(i)=(i-1)*hz(i);  
  end
  zw(mp1)=zL;
  
%построение сетки - z(i) (координаты узловых точек) во всей области
  z(1)=zw(1);
  for i=1:m
      z(i+1)=0.5*(zw(i)+zw(i+1));
  end
  z(mp2)=zw(mp1);
end  % end GRID
      

%////////////////////////////////
function [] = Material_properties_and_boundary_conditions() 
 l_t=fix(it/it_per)+1;
 while(l_t>li)
      l_t=l_t-li;
 end
 
 lambda=zeros(m);
 c_rho=zeros(m);
 
 for i=1:m
     nv = n_lr(i+1);
     if (u(i+1)<u_phase(nv)+0.0001)
        lambda(i)=lambda_M(nv);
        zn_W_w=W_w(nv,u(i),alpha_w,beta_w,gamma_w);             
        zn_Der_W_w=Derivative_W_w(nv,u(i),alpha_w,beta_w);       
        c_M=rho_d(nv)*(c_d(nv)+c_ice*(W_tot(nv)-zn_W_w)+c_w*zn_W_w+...
                       rkappa*zn_Der_W_w);
        cl=c_M;
     else 
        lambda(i)=lambda_T(nv);
        c_T=rho_d(nv)*(c_d(nv)+c_w*W_tot(nv));
        cl=c_T;
     end 
     zn_W_w_phase=W_w(nv,u_phase(nv),alpha_w,beta_w,gamma_w);
     c_rho(i)=cl+rkappa*rho_d(nv)*...
              (W_tot(nv)-zn_W_w_phase)*d_f(u(i)-u_phase(nv));
 end 
  
   f_Ci0=he(l_t).*theta(l_t);   
   f_Pi0=-he(l_t); 
end  %end Material_properties_and_boundary_conditions

%/////////////////////
function [] = SLAU() 

  Matrix_elements();
%!---------------------------
%!     Прогонка по i        !
%!---------------------------
  a=ones(m);
  a_m=ones(m);
  a_p=ones(m);
  b=ones(m);
  for i=1:m
      a(i)=ae(i);
      a_m(i)=ae_im(i);
      a_p(i)=ae_ip(i);
      b(i)=be(i);
  end
         
  xx=progonka(a,a_m,a_p,b);
  for i=1:m
      u(i+1)=xx(i);
  end
end % end SLAU

%////////////////
function [] = Matrix_elements()  
%Вычисление коэффициентов матрицы СЛАУ
	       
  beta=4.0/3.0;
  if(iorder == 1) 
     beta=1.0;
  end
  for i=1:m
      Vol=hz(i);
      ae_0=c_rho(i)*Vol/dt;
      be(i)=S_C(i)*Vol+ae_0*u(i+1);
      ae(i)=ae_0-S_P(i)*Vol;
  end
      
%!Вычисление коэффициентов ae_im(i) и ae_ip(i) во внутренних точках области:
      
%!Вдоль оси z:
  for i=1:mm1
      D_e=2.*lambda(i)*lambda(i+1)/...
      ((hz(i)*lambda(i+1)+hz(i+1)*lambda(i)+eps));      
      ae_ip(i)=D_e+eps;
      ae_im(i+1)=ae_ip(i);
  end
      
%!Вычисление коэффициентов ae_i в граничных точках области:

%!На границе i=0:
  D_e=lambda(1)/(0.5*hz(1))+eps;
  cl=beta*D_e;                  %!ae_im(1)=beta*D_e
  ae_ip_0=cl;                   %!ae_ip_0=ae_im(1)
  ae_im_0=(beta-1.)*ae_ip(1);
  ae_ip(1)=ae_ip(1)+ae_im_0;
  ae_im(1)=cl;                  %!ae_im(1)=ae_im(1)

  if(ibc_i0==1)
     u(1)=ua(z(1),t);
     be(1)=be(1)+ae_im(1)*u(1);
  else
     cl1=ae_ip_0-f_Pi0;
     cl=ae_im(1)/cl1;
     ae(1)=ae(1)-ae_ip_0*cl;   
     ae_ip(1)=ae_ip(1)-ae_im_0*cl;
     be(1)=be(1)+f_Ci0*cl;
  end                                                                
  ae(1)=ae(1)+ae_im(1);
  ae_im(1)=0.0;
         
%!На границе i=mp1:
  D_e=lambda(m)/(0.5*hz(m))+eps;
  cl=beta*D_e;                   %!ae_ip(m)=beta*D_e
  ae_im(mp1)=cl;                 %!ae_im(mp1)=ae_ip(m)
  ae_ip(mp1)=(beta-1.)*ae_im(m); 
  ae_im(m)=ae_im(m)+ae_ip(mp1);
  ae_ip(m)=cl;                   %!ae_ip(m)=ae_ip(m)         

  if(ibc_mp1==1)
	 u(mp2)=ua(z(mp2),t);
     be(m)=be(m)+ae_ip(m)*u(mp2);
  else
     cl1=ae_im(mp1)-f_Pmp1;
     cl=ae_ip(m)/cl1;
     ae(m)=ae(m)-ae_im(mp1)*cl;
     ae_im(m)=ae_im(m)-ae_ip(mp1)*cl;
     be(m)=be(m)+f_Cmp1*cl;
  end
  ae(m)=ae(m)+ae_ip(m);
  ae_ip(m)=0.0;

  for i=1:m
      ae(i)=ae(i)+ae_ip(i)+ae_im(i);
  end
end  % end Matrix_elements

%//////////////////////////////////////////
function [s] = progonka(aa,aa_m,aa_p,bb)
  al=ones(m); 
  bet=ones(m);
  al(1)=aa_p(1)/aa(1);
  bet(1)=bb(1)/aa(1);
  for kp=2:m
      rp=aa(kp)-aa_m(kp)*al(kp-1);
      al(kp)=aa_p(kp)/rp;
      bet(kp)=(bb(kp)+aa_m(kp)*bet(kp-1))/rp;
  end
  s(m)=bet(m);

  for kp=m-1:-1:1 
      s(kp)=al(kp)*s(kp+1)+bet(kp);
  end
end


%////////////////////////
function [s] = d_f(u)
  delta=2;
  s=0;
  if (abs(u) <= delta) 
      s=0.5/delta;
  end
end

%//////////////////////////////////
function [] = W_w_coef()
%РСН
  %u1(1)=-0.3;  Ww_u1(1)=0.14;
  %u2(1)=-1.;   Ww_u2(1)=0.12;
  %u3(1)=-10.;  Ww_u3(1)=0.08;
  
  %u1(2)=-0.3;  Ww_u1(2)=0.12;
  %u2(2)=-1.;   Ww_u2(2)=0.08;
  %u3(2)=-10.;  Ww_u3(2)=0.05;
%New
  u1(1)=-0.3;  Ww_u1(1)=0.1435;
  u2(1)=-1.;   Ww_u2(1)=0.1189;
  u3(1)=-10.;  Ww_u3(1)=0.082;
  
  u1(2)=-0.3;  Ww_u1(2)=0.09506;
  u2(2)=-1.;   Ww_u2(2)=0.078764;
  u3(2)=-10.;  Ww_u3(2)=0.05432;
  
  u1(3)=0.;  Ww_u1(3)=0.;
  u2(3)=0.;  Ww_u2(3)=0.;
  u3(3)=0.;  Ww_u3(3)=0.;

  
  alpha_w=zeros(1,nl);
   beta_w=zeros(1,nl);
  gamma_w=zeros(1,nl);
  
  for i=1:nl
      w1 = abs(Ww_u1(i)-Ww_u2(i));
      w2 = abs(Ww_u1(i)-Ww_u3(i));
      
      if ((w1<0.0001)||(w2<0.0001))
         alpha_w(i) = 0;
         beta_w(i) = 1.0;
         gamma_w(i) = Ww_u1(i);
      else
      beta_w(i)=(Ww_u1(i)*u2(i)*u1(i)-Ww_u2(i)*u2(i)*u1(i)+Ww_u3(i)*u1(i)*u3(i)-...
                 Ww_u3(i)*u2(i)*u3(i)-Ww_u1(i)*u1(i)*u3(i)+Ww_u2(i)*u2(i)*u3(i))/...
                 (Ww_u1(i)*u2(i)+Ww_u2(i)*u3(i)+Ww_u3(i)*u1(i)-Ww_u1(i)*u3(i)-...
                 Ww_u2(i)*u1(i)-Ww_u3(i)*u2(i)); 
      gamma_w(i)=(Ww_u1(i)*beta_w(i)+Ww_u2(i)*u2(i)-Ww_u2(i)*beta_w(i)-...
                 Ww_u1(i)*u1(i))/(u2(i)-u1(i));
      alpha_w(i)=(Ww_u1(i)-gamma_w(i))*(beta_w(i)-u1(i));
   end
  end

% %Вывод на экран коэффициентов функции Ww=Ww(u):
%   for i=1:nl
%       fprintf('lithological = %14.0f \n',i);    
%       fprintf('alpha_w(i) = %16.8f \n',alpha_w(i));
%       fprintf('beta_w(i) = %17.8f \n',beta_w(i));
%       fprintf('gamma_w(i) = %16.8f \n',gamma_w(i));
%   end 
end

%///////////////////////////////////////////////
function [ww] = W_w(nlv,util,alphal,betal,gammal)
  ww=alphal(nlv)/(betal(nlv)-util)+gammal(nlv);
end

%////////////////////////////////////////////////////
function [ww] = Derivative_W_w(nlv,util,alphal,betal)
  ww=alphal(nlv)/((betal(nlv) - util)*(betal(nlv) - util));
end

%///////////////////////////////
function [] = number_lithologic()
 for i=2:m+1
     if ((z(i)>zw_z(1)) && (z(i)<zw_z(2))) n_lr(i) = 1;
     end 
     if ((z(i)>zw_z(2)) && (z(i)<zw_z(3))) n_lr(i) = 2;
     end
     if (z(i)>zw_z(3)) n_lr(i)=3;
     end
 end
end

%/////////////////////
function [] = OUTPUT()
 for i=2:m+1
     im1=i-1;
     xg(im1)=z(i);
     ug(im1)=u(i);
    fprintf(1,'    z(%d)=%4.2f, u(%d)=%8.4f\n',im1,xg(im1),im1,ug(im1));
 end 
 hold on  
%figure;
 plot(xg,ug,'-k');  % r c g b m k y
 xlabel('z');
 ylabel('u');
 grid;  
end %end OUTPUT
end