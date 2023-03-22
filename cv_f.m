function [] = cv_f
  clear all
  global m mm1 mp1 mp2  %m - число КО на [0,xL]   
  global nt
  global kord
  global c_rho conduction_heat_k rho c q_E
  global eps D_e
  global t0 tk t dt li it_per
  global xL powerX
  global rLam lambda 
  global ae ae_ip ae_im be
  global h theta
  global x u
  global xw hx
  global f_Ci0 f_Cmp1
  global f_Pi0 f_Pmp1
  global S_C S_P 

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
%Температура окружающей среды:
    %РСН87-67
%theta(1)=-24.6; theta(2)=-23.2; theta(3)=-19.2; %зима 
%theta(4)=-9.5; theta(5)=3.7; theta(6)=13.1; %весна
%theta(7)=18.0; theta(8)=12.4; theta(9)=4.6; %лето
%theta(10)=-4.6; theta(11)=-16.4; theta(12)=-22.4; %осень
    %RCP2.6
%theta(1)=-21.2; theta(2)=-20.05; theta(3)=-16.3; %зима
%theta(4)=-6.85; theta(5)=6.1; theta(6)=15.25; %весна
%theta(7)=19.9; theta(8)=14.55; theta(9)=7; %лето
%theta(10)=-1.95; theta(11)=-13.5; theta(12)=-19.25; %осень
    %RCP8.5
%theta(1)=-15.5; theta(2)=-14.67; theta(3)=-11.23; %зима
%theta(4)=-2.1; theta(5)=10.53; theta(6)=19.37; %весна
%theta(7)=23.7; theta(8)=18.67; theta(9)=11.43; %лето
%theta(10)=2.8; theta(11)=-8.43; theta(12)=-13.87; %осень
    %Якутск
%theta(1)=-39; theta(2)=-34.5; theta(3)=-20.3; %зима 
%theta(4)=-4.6; theta(5)=7.5; theta(6)=16.2; %весна
%theta(7)=19.3; theta(8)=15.2; theta(9)=5.9; %лето
%theta(10)=-7.8; theta(11)=-27.6; theta(12)=-37.7; %осень
    %RCP2.6
%theta(1)=-35.6; theta(2)=-31.35; theta(3)=-17.4; %зима
%theta(4)=-1.95; theta(5)=9.9; theta(6)=18.35; %весна
%theta(7)=21.2; theta(8)=17.35; theta(9)=8.3; %лето
%theta(10)=-5.15; theta(11)=-24.7; theta(12)=-34.55; %осень
    %RCP8.5
theta(1)=-29.9; theta(2)=-25.97; theta(3)=-12.33; %зима
theta(4)=2.8; theta(5)=14.33; theta(6)=22.47; %весна
theta(7)=25; theta(8)=21.47; theta(9)=12.73; %лето
theta(10)=-0.4; theta(11)=-19.63; theta(12)=-29.17; %осень
%Коэффициент теплообмена h:
h(1)=0.3742; h(2)=0.3595;  h(3)=0.3471;
h(4)=0.3243; h(5)=0.5610;  h(6)= 12.6000;
h(7)=13.4;   h(8)=14.0000;  h(9)= 14.3;
h(10)=0.7275; h(11)=0.4600;  h(12)= 0.4204;

  eps=1.e-20; 
  kord=2;
  q_E=0.0;                  %!0.043    плотность потока из недр
  %РСН67-87          
  %conduction_heat_k=1.55;  %коэффициент теплопроводности мерзлого грунта (у нас пески)
  %rho=1390;                %плотность сухого грунта (cредняя для всех слоев) !1390
  %c=0.22;                  %удельная теплоемкость сухого грунта !
  %Якутск
  conduction_heat_k=1.3;  %коэффициент теплопроводности мерзлого грунта (у нас пески)
  rho=1800;                %плотность сухого грунта (cредняя для всех слоев) !1390
  c=0.25;                  %удельная теплоемкость сухого грунта !
  
  c_rho=c*rho;             %объемная  теплоемкость элемента 
  
  t0=0.0;                  %время начала
  tk=4378540;         %момент времени
  nt=599800; 
  xL=30.0;          %размер расчетной области
  m=300;            %число КО

  dt=(tk-t0)/nt;    %dt=0.002
  it_per=100;  
  li=12;            %месяцев

  powerX=1.0;       %параметр неравномерности сетки по оси х;
  mp2=m+2;          %характерные границчные точки для вывода
  mp1=m+1;
  mm1=m-1;

  for i=1:mp2
        u(i)=-1.0;  %нач температура
  end
  GRID();
end % Data_preparation

%//////////////////
function [] = GRID() 
  xw(1)=0.0;
  %построение сетки - xw(i) (координаты граней КО) во всей области
  xw(mp1)=xL;
    for i=2:m
        cl=(i-1)/m; 
        xw(i)=xL*cl^powerX;
     end
  %построение сетки - x(i) (координаты узловых точек) во всей области
  x(1)=xw(1);
  for i=1:m
      x(i+1)=0.5*(xw(i)+xw(i+1));
  end
  x(mp2)=xw(mp1);
  for i=1:m
      hx(i)=xw(i+1)-xw(i);  %размеры КО 
  end
end  % end GRID
      
%////////////////////////////////
function [] = Material_properties_and_boundary_conditions() 
  for i=1:m
      lambda(i)=conduction_heat_k; %коэффициент теплоемкости
      rLam(i)=c_rho;               %объемная  теплоемкость элемента 
  end
  for i=1:m               %s - (средняя)мощность внутренних источников тепла
      S_C(i)=0.0;         %нет внутренних источников
      S_P(i)=0.0;         % т.к. не зависит от температуры
  end
  l_t=fix(it/it_per)+1;
  while(l_t>li)
     l_t=l_t-li;
  end
    f_Ci0=h(l_t).*theta(l_t);%значение fc на границах c воздухом (конвективный темплообмен)   
    f_Pi0=-h(l_t);          %значения fp на границах с воздухом
    f_Cmp1=q_E;             %значение fc на границах c недрами (известное значение плотности теплового потока
    f_Pmp1=0.0;             %значения fp на границах с недрами
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
  beta=4.0/3.0;%коэффициент в выражениях для аппроксимации на границах (равен 1 и 4/3,
%соответственно, для первого и второго порядков аппроксимации);
  for i=1:m
      Vol=hx(i);
      ae_0=rLam(i)*Vol/dt;          %a маленькое р индекс 0
      be(i)=S_C(i)*Vol+ae_0*u(i+1); %b(i)=Di
      ae(i)=ae_0-S_P(i)*Vol;        %a маленькое р, но почему-то без двух слагаемых
  end  
%!Вычисление коэффициентов ae_im(i) и ae_ip(i) во внутренних точках области:    
%!Вдоль оси x:
  for i=1:mm1
      D_e=2.*lambda(i)*lambda(i+1)/...
      ((hx(i)*lambda(i+1)+hx(i+1)*lambda(i)+eps));   
      ae_ip(i)=D_e+eps; %аЕ
      ae_im(i+1)=ae_ip(i); %aW
  end     
%!Вычисление коэффициентов ae_i в граничных точках области:
%!На границе i=0:
  D_e=lambda(1)/(0.5*hx(1))+eps;
  cl=beta*D_e;                  
  ae_ip_0=cl;                   
  ae_im_0=(beta-1.)*ae_ip(1);
  ae_ip(1)=ae_ip(1)+ae_im_0;
  ae_im(1)=cl;                  
 %задана плотность теплового потока
  cl1=ae_ip_0-f_Pi0;
  cl=ae_im(1)/cl1;     
  ae(1)=ae(1)-ae_ip_0*cl;   
  ae_ip(1)=ae_ip(1)-ae_im_0*cl;
  be(1)=be(1)+f_Ci0*cl;                                                           
  ae(1)=ae(1)+ae_im(1);
  ae_im(1)=0.0;
         
%!На границе i=mp1:
  D_e=lambda(m)/(0.5*hx(m))+eps;
  cl=beta*D_e;               
  ae_im(mp1)=cl;                
  ae_ip(mp1)=(beta-1.)*ae_im(m); 
  ae_im(m)=ae_im(m)+ae_ip(mp1);
  ae_ip(m)=cl;                           
  cl1=ae_im(mp1)-f_Pmp1;
  cl=ae_ip(m)/cl1;          
  ae(m)=ae(m)-ae_im(mp1)*cl;
  ae_im(m)=ae_im(m)-ae_ip(mp1)*cl;
  be(m)=be(m)+f_Cmp1*cl;
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

%////////////////
function [] = OUTPUT()
  for i=2:m+1          %отбрасываем первую и последнюю точку
      im1=i-1;
      xg(im1)=x(i);
      ug(im1)=u(i);
      fprintf(1,'    x(%d)=%4.2f, u(%d)=%8.4f\n'...
      ,im1,xg(im1),im1,ug(im1));
  end
 hold on 
 %set('Color', 'w');
 plot(xg,ug,'-m'); % r c g b m k y
 xlabel('z');
 ylabel('u');
 grid on; 
 end %end OUTPUT

end