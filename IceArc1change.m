function [] = IceArc1change
  clear all
  global nl
  global m mm1 mp1 mp2  %m - число КО на [0,zL]   
  global ibc_i0 ibc_mp1
  global nt
  global iorder
  global rho_dT rho_w rho_water c_dwaterT c_dT lambda_MT lambda_MTrho lambda_TT lambda_TTemp u_phase
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
  global rkappa c_w q_E
  
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
  nl=1; % число литологических слоев
  li=12;
%Температура окружающей среды:
%Северный полюс
 %theta(1)=-30.6735294117647; theta(2)=-31.11875; theta(3)=-31.37142857;
 %theta(4)=-23.26; theta(5)=-10.21142857; theta(6)=-0.562857143; 
 %theta(7)=1.877142857; theta(8)=0.067567568; theta(9)=-8.661764706; 
 %theta(10)=-18.44242424; theta(11)=-25.37142857; theta(12)=-28.725; 
%RCP2.6
 %theta(1)=-27.6735294117647; theta(2)=-28.11875; theta(3)=-28.37142857;
 %theta(4)=-20.26; theta(5)=-7.21142857; theta(6)=2.437142857; 
 %theta(7)=4.877142857; theta(8)=3.067567568; theta(9)=-5.661764706; 
 %theta(10)=-15.44242424; theta(11)=-22.37142857; theta(12)=-25.725;
%RCP8.5
 %theta(1)=-20.6735294117647; theta(2)=-21.11875; theta(3)=-21.37142857;
 %theta(4)=-13.26; theta(5)=-0.21142857; theta(6)=9.437142857; 
 %theta(7)=11.877142857; theta(8)=10.067567568; theta(9)=1.338235294; 
 %theta(10)=-8.44242424; theta(11)=-15.37142857; theta(12)=-18.725;
%1.5 C
 theta(1)=-26.17352941; theta(2)=-26.61875; theta(3)=-26.87142857;
 theta(4)=-18.76; theta(5)=-5.711428571; theta(6)=3.937142857; 
 theta(7)=6.377142857; theta(8)=4.567567568; theta(9)=-4.161764706; 
 theta(10)=-13.94242424; theta(11)=-20.87142857; theta(12)=-24.225;
%Коэффициент теплообмена h:
 he(1)=1.458890709; he(2)=1.45918942;  he(3)=1.461857699;
 he(4)=1.464902885; he(5)=1.460594677;  he(6)=1.459117811;
 he(7)=1.455429996; he(8)=1.451234234;  he(9)=1.451735837;
 he(10)=1.453491003; he(11)=1.458843801;  he(12)=1.463016893;
%------------------------------------------------------------ 
%Земля Франца-Иосифа
 %theta(1)=-23.0148; theta(2)=-23.06; theta(3)=-23.1567;
 %theta(4)=-18.3339; theta(5)=-8.985; theta(6)=-1.38644; 
 %theta(7)=0.767797; theta(8)=0.177586; theta(9)=-2.6386; 
 %theta(10)=-10.3196; theta(11)=-16.4293; theta(12)=-21.0293; 
%RCP2.6
 %theta(1)=-19.0148; theta(2)=-19.06; theta(3)=-19.1567;
 %theta(4)=-14.3339; theta(5)=-4.985; theta(6)=2.61356; 
 %theta(7)=4.767797; theta(8)=4.177586; theta(9)=1.3614; 
 %theta(10)=-6.3196; theta(11)=-12.4293; theta(12)=-17.0293; 
%RCP8.5
 %theta(1)=-11.0148; theta(2)=-11.06; theta(3)=-11.1567;
 %theta(4)=-6.3339; theta(5)=3.015; theta(6)=10.61356; 
 %theta(7)=12.767797; theta(8)=12.177586; theta(9)=9.3614; 
 %theta(10)=1.6804; theta(11)=-4.4293; theta(12)=-9.0293;
 %Земля Франца-Иосифа Коэффициент теплообмена h:
 %he(1)=1.389296508; he(2)=1.104079959;  he(3)=0.951550648;
 %he(4)=0.860422844; he(5)=0.806737701;  he(6)=0.935798623;
 %he(7)=3.351470977; he(8)=10.30505755;  he(9)=7.792177559;
 %he(10)=3.956285514; he(11)=2.557027996;  he(12)=1.880398672;
 
  iorder=2;
  
  %удельная теплота фазового перехода
  rkappa=79.4; %
  %78.63 для -2 градусов
  %79.4 для 0 градусов
 
  q_E=0.0;   % q_E=0.043;
      
  it=0;
  t0=0.0;
 %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 %50 лет
  %tk=429970; nt=58900;
  %tk=430700; nt=59000;
  %tk=431430; nt=59100;
  %tk=432160; nt=59200;
  %tk=432890; nt=59300;
  %tk=433620; nt=59400;
  %tk=434350; nt=59500;
  %tk=435080; nt=59600;
  %tk=435810; nt=59700;
  %tk=436540; nt=59800;
  %tk=437270; nt=59900;
  tk=438000; nt=60000;
 %Январь
  %tk=730; nt=100;
  %tk=35770; nt=4900;
  %tk=79570; nt=10900;
  %tk=167170; nt=22900;
  %tk=298570; nt=40900;
  %tk=429970; nt=58900;
 %Июль
  %tk=5110; nt=700;
  %tk=40150; nt=5500;
  %tk=83950; nt=11500;
  %tk=171550; nt=23500;
  %tk=302950; nt=41500;
  %tk=434350; nt=59500;
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
    %начальная температура
  for i=1:mp2  u(i)=-20;      end
    
  
%Объемные источники тепла отсутствуют:  
  for i=1:m
      S_C(i)=0.0;      
      S_P(i)=0.0;
  end     

%Свойства льда и воды
        %плотность
  rho_w=1000;
        %теплоемкость
  c_w=1.007; %теперь как с_dM
        %теплопроводность
  lambda_TT=0.489;%теплопроводность воды
  
  u_phase(1)=0.0;  %температура фазового перехода, для морской воды -1.8

end % Data_preparation

%//////////////////
function [] = GRID()
%координаты границ литологических разрезов по оси z;
  %hz_c=zL/m;
  zw_z(1)=0.0;
  zw_z(2)=zL; %всего 1 слой

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
     rho_dT(i)=916.8/(1+0.000158*u(i));
     c_dT(i)=(2.114+0.007787*u(i))*0.2388;
     lambda_MT(i)=2.24*(1-0.0048*u(i)); 
     lambda_MTrho(i)=(lambda_MT(i)-0.0057*(916.8-rho_dT(i)))*0.86; 
     %-------------------------------
     rho_water(i)=995.7/(0.984+0.000483*u(i));
     c_dwaterT(i)=(4194-1.15*u(i)-0.015*u(i)*u(i))*0.2388; %для температуры от 10 до 100 градусов
     lambda_TTemp(i)=0.553*(1+0.003*u(i))*0.86; %10-100
     if (u(i+1)<u_phase(nv)+0.0001)
        lambda(i)=lambda_MTrho(i);       
        c_M=rho_dT(i)*c_dT(i);
        cl=c_M;
     else 
        if (u(i)<10+0.0001)
            lambda(i)=lambda_TT;
            c_T=rho_w*c_w;
            cl=c_T;
        else  
            lambda(i)=lambda_TTemp(i);
            c_T=rho_water(i)*c_dwaterT(i);
            cl=c_T;
        end
     end 
     if (u(i+1)<u_phase(nv)+0.0001)
        c_rho(i)=cl+rkappa*rho_dT(i)*d_f(u(i)-u_phase(nv));
     else
        if (u(i)<10+0.0001)
            c_rho(i)=cl+rkappa*rho_w*d_f(u(i)-u_phase(nv));
        else  
            c_rho(i)=cl+rkappa*rho_water(i)*d_f(u(i)-u_phase(nv));
        end
     end
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

%///////////////////////////////
function [] = number_lithologic()
 for i=2:m+1
     if ((z(i)>zw_z(1)) && (z(i)<zw_z(2))) n_lr(i) = 1;
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
    hF = gcf;
    hF.Position = [100 200 1280 720];
 plot(xg,ug,'-r');  % r c g b m k y
 xlabel('z');
 ylabel('u');
 grid;  
 end %end OUTPUT
end