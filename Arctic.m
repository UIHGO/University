clear;
close all;
clc;

if (~exist('g','var') || ~isa(g,'GCSAL.GCSAL'))
    % The gcsal.h5 and gcsal.h5.info.mat files are available for download from the
    % website and should be placed in the h5_data directory.

    %%%%%%%%%%%%%% CHANGE THESE %%%%%%%%%%%%%%%%%%
    % Set this to wherever you put the gcsal.h5 file and gcsal.h5.info.mat
    % files downloaded from dropbox
    h5_dir = 'C:\Users\Владимир\Desktop\Dourse\GCSAL-master\h5_data';

    
    % Directory to code. The folder +GCSAL which contains this file should be
    % in this directory
    codebase_dir = 'C:\Users\Владимир\Desktop\Dourse';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Full path to .mat file with h5 info
    h5_file = fullfile(h5_dir, 'gcsal.h5');
    h5_mat_file = [h5_file '.info.mat'];

    % Set up Matlab path
    addpath(genpath(codebase_dir))

    %% Load GCSAL object from .mat file

    % This requires about 6 gb of RAM just to hold all of the header
    % information in memory

    % Normally you should load the GCSAL object from the .mat file but if it
    % doesn't exist on your path you can use the .h5 file. After using the .h5
    % file a .mat file will be created automatically for subsequent use
    if ~exist(h5_mat_file, 'file')
        g = GCSAL.GCSAL(h5_file);
    else
        g = GCSAL.GCSAL(h5_mat_file);
        g.h5_fname = h5_file;
    end
end

% Introduction: Printout the list of variables that are in the header data
% and the entries data.
% header data is for a single balloon launch - things like time, date, and location
% entries data is the measurements of the baloon - wind speed, pressure, etc.
g.defs.header.params
g.defs.entries.params

% For more details look at a single parameter
% g.defs.entries.params.wspd
% g.defs.entries.params.wdir
% g.defs.entries.params.gph
% g.defs.entries.params.temp
% g.defs.entries.params.press
%%


%% Stations_search
 %Антарктика -80 -70 -30 -25 1977
  stations7 =  g.station_search('LatLong', [72 74 50 55]);
  %%новая земля 72 74 50 55 1987
  %%исландия 78 80 11 13 1995 шпицберген
  
  yearall=1977; % начальный год
  gphvys=0.875; %На высоте
  
%%Температура, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries2 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys]); 
%считаем количество лет за которое записывались данные 
years = entries2.year(end)-entries2.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys], yearall+1}); 
for j=2:(2017-yearall) 
entries_years(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries3(j)=0;l=0; 
for i=1:length(entries_years(j).temp) 
if entries_years(j).temp(i) <=0 || entries_years(j).temp(i) >=0 
temp=entries_years(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries3(j) = entries3(j) + temp; 
end 
entries3(j)=entries3(j)/(length(entries_years(j).temp)-l); 
end 
%строим график средних температур за 30 лет 
%Прямая 
p=polyfit(year,entries3,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
f1=figure('Position',[50 50 800 800]);
set(f1, 'Color', 'w');
plot1=plot(year,entries3,x1,y1); 
ylabel('Температура, C'); xlabel(' '); 
set(plot1(1),'linewidth',1.1,'Marker','*');
set(plot1(2),'color','r','linewidth',1.5,'linestyle','--');
grid on; 
hold on;

%%Температура на высоте 5, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries21 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys+4.125]); 
%считаем количество лет за которое записывались данные 
years = entries21.year(end)-entries21.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years1(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+1}); 
for j=2:(2017-yearall) 
entries_years1(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries31(j)=0;l=0; 
for i=1:length(entries_years1(j).temp) 
if entries_years1(j).temp(i) <=0 || entries_years1(j).temp(i) >=0 
temp=entries_years1(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries31(j) = entries31(j) + temp; 
end 
entries31(j)=entries31(j)/(length(entries_years1(j).temp)-l); 
end 
%строим график средних температур за 40 лет 
%Прямая 
p=polyfit(year,entries31,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
plot1=plot(year,entries31,x1,y1); 
ylabel('Температура, C'); xlabel(' '); 
set(plot1(1),'color','k','linewidth',1.1,'Marker','^'); 
set(plot1(2),'color','m','linewidth',1.5,'linestyle',':');
grid on; 
hold on; 

%%Температура на высоте 10, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries22 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys+9.125]); 
%считаем количество лет за которое записывались данные 
years = entries22.year(end)-entries22.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years2(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+1}); 
for j=2:(2017-yearall) 
entries_years2(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries32(j)=0;l=0; 
for i=1:length(entries_years2(j).temp) 
if entries_years2(j).temp(i) <=0 || entries_years2(j).temp(i) >=0 
temp=entries_years2(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries32(j) = entries32(j) + temp; 
end 
entries32(j)=entries32(j)/(length(entries_years2(j).temp)-l); 
end 
%строим график средних температур за 40 лет 
%Прямая 
p=polyfit(year,entries32,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
plot1=plot(year,entries32,x1,y1); 
set(plot1(1),'linewidth',1.1,'Marker','o'); 
set(plot1(2),'color','c','linewidth',1.5,'linestyle','-.');
%legend('Поверхность','линейный тренд','5 км','линейный тренд','10 км','линейный тренд','location','northeast'); legend('boxon');
ylabel('Температура, C'); xlabel(' '); 
%X=['Изменение температуры вблизи Атлантического побережья Антарктиды за ', num2str(2017-yearall), 'лет']
%title(X); 
grid on; 
hold on; 

%print('-dpng','-r800','Температура 1-10');

%%Температура на высоте 15, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries24 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys+14.125]); 
%считаем количество лет за которое записывались данные 
years = entries24.year(end)-entries24.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years4(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+1}); 
for j=2:(2017-yearall) 
entries_years4(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries34(j)=0;l=0; 
for i=1:length(entries_years4(j).temp) 
if entries_years4(j).temp(i) <=0 || entries_years4(j).temp(i) >=0 
temp=entries_years4(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries34(j) = entries34(j) + temp; 
end 
entries34(j)=entries34(j)/(length(entries_years4(j).temp)-l); 
end 
%строим график средних температур за 40 лет 
%Прямая 
p=polyfit(year,entries34,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
%f11=figure('Position',[50 20 800 800]);
%set(f11, 'Color', 'w');
plot1=plot(year,entries34,x1,y1); 
set(plot1(1),'linewidth',1.1,'Marker','diamond');
set(plot1(2),'color','r','linewidth',1.5,'linestyle','--');
ylabel('Температура, C'); xlabel(' '); 
grid on; 
hold on; 

%%Температура на высоте 20, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries23 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys+19.125]); 
%считаем количество лет за которое записывались данные 
years = entries23.year(end)-entries23.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years3(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+1}); 
for j=2:(2017-yearall) 
entries_years3(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries33(j)=0;l=0; 
for i=1:length(entries_years3(j).temp) 
if entries_years3(j).temp(i) <=0 || entries_years3(j).temp(i) >=0 
temp=entries_years3(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries33(j) = entries33(j) + temp; 
end 
entries33(j)=entries33(j)/(length(entries_years3(j).temp)-l); 
end 
%строим график средних температур за 40 лет 
%Прямая 
p=polyfit(year,entries33,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
plot1=plot(year,entries33,x1,y1); 
ylabel('Температура, C'); xlabel(' '); 
set(plot1(1),'color','k','linewidth',1.1,'Marker','>'); 
set(plot1(2),'color','m','linewidth',1.5,'linestyle',':');
%X=['Изменение температуры вблизи Атлантического побережья Антарктиды за ', num2str(2017-yearall) ,' лет '] 
%title(X); 
grid on; 
hold on; 

%%Температура на высоте 25, C%% 
%находим показания температуры по этой станции на высоте в пределах до x км 
entries25 = g.query(stations7, {'temp','year'}, 'gph', [0 gphvys+24.125]); 
%считаем количество лет за которое записывались данные 
years = entries25.year(end)-entries25.year(1); 
%делаем массив просто из номеров годов нужно для графика в конце 
year(1)=0; 
for i=1:(2017-yearall) 
year(i)=yearall+i; 
end 
%заводим массив данных по каждому году чтобы каждая переменная отвечала 
%за 1 год 
entries_years5(1)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+1}); 
for j=2:(2017-yearall) 
entries_years5(j)=g.query(stations7, {'temp','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+j}); 
end 
%заводим цикл который будет складывать все температуры за год исключая 
%переменные NaN (l-считает количество NaN) 
for j=1:(2017-yearall) 
entries35(j)=0;l=0; 
for i=1:length(entries_years5(j).temp) 
if entries_years5(j).temp(i) <=0 || entries_years5(j).temp(i) >=0 
temp=entries_years5(j).temp(i); 
else 
temp=0; l=l+1; 
end 
entries35(j) = entries35(j) + temp; 
end 
entries35(j)=entries35(j)/(length(entries_years5(j).temp)-l); 
end 
%строим график средних температур за 40 лет 
%Прямая 
p=polyfit(year,entries35,1); 
x1=linspace(year(1),year(end)+7,100); 
y1=polyval(p,x1); 
plot1=plot(year,entries35,x1,y1); 
ylabel('Температура, C'); xlabel(' '); 
set(plot1(1),'linewidth',1.1,'Marker','pentagram');
set(plot1(2),'color','c','linewidth',1.5,'linestyle','-.');
legend('Поверхность','линейный тренд','5 км','линейный тренд','10 км','линейный тренд','15 км','линейный тренд','20 км','линейный тренд','25 км','линейный тренд','location','EastOutside'); legend('boxon');
X=['Изменение температуры на о.Южный за ', num2str(2017-yearall) ,' лет '] 
title(X); 
grid on; 
hold on;

% print('-dpng','-r800','Температура 15-25');

%......................................................................%


  %%ДАВЛЕНИЕ 1 км%%
  entries2P = g.query(stations7, {'press','year'}, 'gph', [0 gphvys]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P.year(length(entries2P.year))-entries2P.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P(j)=0;l=0;
   for i=1:length(entries_yearsP(j).press)
      if entries_yearsP(j).press(i) <=0 || entries_yearsP(j).press(i) >=0
          tempp=entries_yearsP(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P(j) = entries32P(j) + tempp;
   end
   entries32P(j)=entries32P(j)/(length(entries_yearsP(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  f2=figure('Position',[50 50 800 800]); 
  set(f2, 'Color', 'w');
  plot2=plot(yearP,entries32P,x1,y1);
  set(plot2(1),'linewidth',1.3,'Marker','*');
  set(plot2(2),'color','r','linewidth',1.5,'linestyle','--');
  ylabel('Давление, Па'); xlabel(' '); 
  grid on;
  hold on;
  
  
    %%ДАВЛЕНИЕ 5 км%%
  entries2P2 = g.query(stations7, {'press','year'}, 'gph', [0 gphvys+4.125]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P2.year(length(entries2P2.year))-entries2P2.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP2(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP2(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P2(j)=0;l=0;
   for i=1:length(entries_yearsP2(j).press)
      if entries_yearsP2(j).press(i) <=0 || entries_yearsP2(j).press(i) >=0
          tempp=entries_yearsP2(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P2(j) = entries32P2(j) + tempp;
   end
   entries32P2(j)=entries32P2(j)/(length(entries_yearsP2(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P2,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  plot2=plot(yearP,entries32P2,x1,y1);
  set(plot2(1),'color','k','linewidth',1.3,'Marker','^');
  set(plot2(2),'color','m','linewidth',1.5,'linestyle',':');
  ylabel('Давление, Па'); xlabel(' ');
  %legend('Поверхность','линейный тренд','5 км','линейный тренд','location','northeast'); legend('boxon');
  %X=['Изменение давления вблизи Атлантического побережья Антарктиды за ', num2str(2017-yearall) ,' лет ']
  %title(X);
  grid on;
  hold on;
    print('-dpng','-r800','Давление 1 - 5');
  
    %%ДАВЛЕНИЕ 10 км %%
  entries2P1 = g.query(stations7, {'press','year'}, 'gph', [0 gphvys+9.125]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P1.year(length(entries2P1.year))-entries2P1.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP1(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP1(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P1(j)=0;l=0;
   for i=1:length(entries_yearsP1(j).press)
      if entries_yearsP1(j).press(i) <=0 || entries_yearsP1(j).press(i) >=0
          tempp=entries_yearsP1(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P1(j) = entries32P1(j) + tempp;
   end
   entries32P1(j)=entries32P1(j)/(length(entries_yearsP1(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P1,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  %f2=figure('Position',[115 115 800 1000]); 
  %set(f2, 'Color', 'w');
  plot2=plot(yearP,entries32P1,x1,y1);
  set(plot2(1),'linewidth',1.3,'Marker','o');
  set(plot2(2),'color','r','linewidth',1.5,'linestyle','-.');
  ylabel('Давление, Па'); xlabel(' '); 
  grid on;
  hold on;
   
     %%ДАВЛЕНИЕ 15 км %%
  entries2P3 = g.query(stations7, {'press','year'}, 'gph', [0 gphvys+14.125]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P3.year(length(entries2P3.year))-entries2P3.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP3(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP3(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P3(j)=0;l=0;
   for i=1:length(entries_yearsP3(j).press)
      if entries_yearsP3(j).press(i) <=0 || entries_yearsP3(j).press(i) >=0
          tempp=entries_yearsP3(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P3(j) = entries32P3(j) + tempp;
   end
   entries32P3(j)=entries32P3(j)/(length(entries_yearsP3(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P3,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  plot2=plot(yearP,entries32P3,x1,y1);
  set(plot2(1),'linewidth',1.3,'Marker','diamond');
  set(plot2(2),'color','r','linewidth',1.5,'linestyle','--');
  ylabel('Давление, Па'); xlabel(' '); 
  %legend('10 км','линейный тренд','15 км','линейный тренд','location','northeast'); legend('boxon');
  %X=['Изменение давления вблизи Атлантического побережья Антарктиды за ', num2str(2017-yearall) ,' лет ']
  %title(X);
  grid on;
  hold on;
   print('-dpng','-r800','Давление 10 - 15');
   
   
       %%ДАВЛЕНИЕ 20 км %%
  entries2P4 = g.query(stations7, {'press','year'}, 'gph', [0 gphvys+19.125]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P4.year(length(entries2P4.year))-entries2P4.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP4(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP4(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P4(j)=0;l=0;
   for i=1:length(entries_yearsP4(j).press)
      if entries_yearsP4(j).press(i) <=0 || entries_yearsP4(j).press(i) >=0
          tempp=entries_yearsP4(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P4(j) = entries32P4(j) + tempp;
   end
   entries32P4(j)=entries32P4(j)/(length(entries_yearsP4(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P4,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  %f2=figure('Position',[115 115 800 1000]); 
  %set(f2, 'Color', 'w');
  plot2=plot(yearP,entries32P4,x1,y1);
  set(plot2(1),'color','k','linewidth',1.3,'Marker','>');
  set(plot2(2),'color','m','linewidth',1.5,'linestyle',':');
  ylabel('Давление, Па'); xlabel(' '); 
  grid on;
  hold on;
   
     %%ДАВЛЕНИЕ 25 км %%
  entries2P5 = g.query(stations7, {'press','year'}, 'gph', [0 gphvys+24.125]);
  %считаем количество лет за которое записывались данные
  yearsP = entries2P5.year(length(entries2P5.year))-entries2P5.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearP(1)=0;
  for i=1:(2017-yearall)
      yearP(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsP5(1)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsP5(j)=g.query(stations7, {'press','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32P5(j)=0;l=0;
   for i=1:length(entries_yearsP5(j).press)
      if entries_yearsP5(j).press(i) <=0 || entries_yearsP5(j).press(i) >=0
          tempp=entries_yearsP5(j).press(i);
      else 
          tempp=0; l=l+1;
      end
       entries32P5(j) = entries32P5(j) + tempp;
   end
   entries32P5(j)=entries32P5(j)/(length(entries_yearsP5(j).press)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pP=polyfit(yearP,entries32P5,1);
  x1=linspace(yearP(1),yearP(end)+7,100);
  y1=polyval(pP,x1);
  plot2=plot(yearP,entries32P5,x1,y1);
  set(plot2(1),'linewidth',1.3,'Marker','pentagram');
  set(plot2(2),'color','c','linewidth',1.5,'linestyle','-.');
  ylabel('Давление, Па'); xlabel(' '); 
  legend('Поверхность','линейный тренд','5 км','линейный тренд','10 км','линейный тренд','15 км','линейный тренд','20 км','линейный тренд','25 км','линейный тренд','location','EastOutside'); legend('boxon');
  X=['Изменение давления на о.Южный за ', num2str(2017-yearall) ,' лет ']
  title(X);
  grid on;
  hold on;
   print('-dpng','-r800','Давление 20 - 25');
   
   
   %......................................................................%
   
   
   
  
   
  %%СКОРОСТЬ ВЕТРА%%
  entries2W = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys]);
  %считаем количество лет за которое записывались данные
  yearsW = entries2W.year(length(entries2W.year))-entries2W.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearW(1)=0;
  for i=1:(2017-yearall)
      yearW(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsW(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsW(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32W(j)=0;l=0;
   for i=1:length(entries_yearsW(j).wspd)
      if entries_yearsW(j).wspd(i) <=0 || entries_yearsW(j).wspd(i) >=0
          tempW=entries_yearsW(j).wspd(i);
      else 
          tempW=0; l=l+1;
      end
       entries32W(j) = entries32W(j) + tempW;
   end
   entries32W(j)=entries32W(j)/(length(entries_yearsW(j).wspd)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pW=polyfit(yearW,entries32W,1);
  x1=linspace(yearW(1),yearW(end)+7,100);
  y1=polyval(pW,x1);
  f3=figure('Position',[50 50 800 800]); 
  set(f3, 'Color', 'w');
  %subplot(2,1,1);
  plot3=plot(yearW,entries32W,x1,y1);
  ylabel('Скорость ветра, м/с'); xlabel(' '); 
  set(plot3(1),'linewidth',1.3,'Marker','*');
  set(plot3(2),'color','r','linewidth',1.5,'linestyle','--');
  legend('График','Линейная','location','southeast'); legend('boxon');
  X=['Изменение скорости ветра на о.Южный за ', num2str(2017-yearall) ,' лет']
  title(X); 
  grid on;
  hold on;
 
   %%СКОРОСТЬ ВЕТРА на 5%%
  entries2WW = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys+4.125]);
  %считаем количество лет за которое записывались данные
  yearsWW = entries2WW.year(length(entries2WW.year))-entries2WW.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearWW(1)=0;
  for i=1:(2017-yearall)
      yearWW(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsWW(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsWW(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+4.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32WW(j)=0;l=0;
   for i=1:length(entries_yearsWW(j).wspd)
      if entries_yearsWW(j).wspd(i) <=0 || entries_yearsWW(j).wspd(i) >=0
          tempWW=entries_yearsWW(j).wspd(i);
      else 
          tempWW=0; l=l+1;
      end
       entries32WW(j) = entries32WW(j) + tempWW;
   end
   entries32WW(j)=entries32WW(j)/(length(entries_yearsWW(j).wspd)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pWW=polyfit(yearWW,entries32WW,1);
  x1=linspace(yearWW(1),yearWW(end)+7,100);
  y1=polyval(pWW,x1);
  %subplot(2,1,2);  
  plot3=plot(yearWW,entries32WW,x1,y1);
  ylabel('Скорость ветра, м/с'); xlabel(' '); 
  set(plot3(1),'color','k','linewidth',1.3,'Marker','^');
  set(plot3(2),'color','m','linewidth',1.5,'linestyle',':');
  %legend('График','Линейная','location','southeast'); legend('boxon');
  %title(X); 
  grid on;
  hold on;
  %print('-dpng','-r800','Скорость ветра');
  
  
   %%СКОРОСТЬ ВЕТРА 10%%
  entries2WW1 = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys+9.125]);
  %считаем количество лет за которое записывались данные
  yearsWW = entries2WW1.year(length(entries2WW1.year))-entries2WW1.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearWW(1)=0;
  for i=1:(2017-yearall)
      yearWW(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsWW1(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsWW1(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+9.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32WW1(j)=0;l=0;
   for i=1:length(entries_yearsWW1(j).wspd)
      if entries_yearsWW1(j).wspd(i) <=0 || entries_yearsWW1(j).wspd(i) >=0
          tempWW=entries_yearsWW1(j).wspd(i);
      else 
          tempWW=0; l=l+1;
      end
       entries32WW1(j) = entries32WW1(j) + tempWW;
   end
   entries32WW1(j)=entries32WW1(j)/(length(entries_yearsWW1(j).wspd)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pWW=polyfit(yearWW,entries32WW1,1);
  x1=linspace(yearWW(1),yearWW(end)+7,100);
  y1=polyval(pWW,x1);
  plot3=plot(yearWW,entries32WW1,x1,y1);
  ylabel('Скорость ветра, м/с'); xlabel(' '); 
  set(plot3(1),'linewidth',1.3,'Marker','o');
  set(plot3(2),'color','c','linewidth',1.5,'linestyle','-.'); 
  grid on;
  hold on; 
  
       %%СКОРОСТЬ ВЕТРА 15%%
  entries2WW4 = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys+14.125]);
  %считаем количество лет за которое записывались данные
  yearsWW = entries2WW4.year(length(entries2WW4.year))-entries2WW4.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearWW(1)=0;
  for i=1:(2017-yearall)
      yearWW(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsWW4(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsWW4(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+14.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32WW4(j)=0;l=0;
   for i=1:length(entries_yearsWW4(j).wspd)
      if entries_yearsWW4(j).wspd(i) <=0 || entries_yearsWW4(j).wspd(i) >=0
          tempWW=entries_yearsWW4(j).wspd(i);
      else 
          tempWW=0; l=l+1;
      end
       entries32WW4(j) = entries32WW4(j) + tempWW;
   end
   entries32WW4(j)=entries32WW4(j)/(length(entries_yearsWW4(j).wspd)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pWW=polyfit(yearWW,entries32WW4,1);
  x1=linspace(yearWW(1),yearWW(end)+7,100);
  y1=polyval(pWW,x1);
  plot3=plot(yearWW,entries32WW4,x1,y1);
  ylabel('Скорость ветра, м/с'); xlabel(' '); 
  set(plot3(1),'linewidth',1.3,'Marker','diamond');
  set(plot3(2),'color','r','linewidth',1.5,'linestyle','--'); 
  grid on;
  hold on;
  
  %%СКОРОСТЬ ВЕТРА 20%%
  %entries2WWW5 = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys+19.125]);
  %%считаем количество лет за которое записывались данные
  %yearsWW = entries2WWW5.year(length(entries2WWW5.year))-entries2WWW5.year(1);
  %%делаем массив просто из номеров годов нужно для графика в конце
  %%yearWW(1)=0;
  %for i=1:(2017-yearall)
  %    yearWW(i)=yearall+i;
  %end
  %%заводим массив данных по каждому году чтобы каждая переменная отвечала
  %%за 1 год
  %ntries_yearsWWW5(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+1});
  %for j=2:(2017-yearall)
  %entries_yearsWWW5(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+19.125], yearall+j});
  %end
  %%заводим цикл который будет складывать все температуры за год исключая
  %%переменные NaN (l-считает количество NaN)
  %for j=1:(2017-yearall)
  %    entries32WWW5(j)=0;l=0;
  % for i=1:length(entries_yearsWWW5(j).wspd)
  %    if entries_yearsWWW5(j).wspd(i) <=0 || entries_yearsWWW5(j).wspd(i) >=0
  %        tempWW=entries_yearsWWW5(j).wspd(i);
  %    else 
  %        tempWW=0; l=l+1;
  %    end
  %     entries32WWW5(j) = entries32WWW5(j) + tempWW;
  % end
  % entries32WWW5(j)=entries32WWW5(j)/(length(entries_yearsWWW5(j).wspd)-l);
  %end
  %строим график средних температур за 30 лет 
  %Прямая
  %pWW=polyfit(yearWW,entries32WWW5,1);
  %x1=linspace(yearWW(1),yearWW(end)+5,100);
  %y1=polyval(pWW,x1);
  %plot3=plot(yearWW,entries32WWW5,x1,y1);
  %ylabel('Скорость ветра, м/с'); xlabel(' '); 
  %set(plot3(1),'linewidth',1.3,'Marker','>','Color','k');
  %set(plot3(2),'color','m','linewidth',1.5,'linestyle','-.'); 
  %grid on;
  %hold on;
  
   %%СКОРОСТЬ ВЕТРА 25%%
  entries2WW6 = g.query(stations7, {'wspd','year'}, 'gph', [0 gphvys+24.125]);
  %считаем количество лет за которое записывались данные
  yearsWW = entries2WW6.year(length(entries2WW6.year))-entries2WW6.year(1);
  %делаем массив просто из номеров годов нужно для графика в конце
  yearWW(1)=0;
  for i=1:(2017-yearall)
      yearWW(i)=yearall+i;
  end
  %заводим массив данных по каждому году чтобы каждая переменная отвечала
  %за 1 год
  entries_yearsWW6(1)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+1});
  for j=2:(2017-yearall)
  entries_yearsWW6(j)=g.query(stations7, {'wspd','year'}, {'gph','year'}, {[0 gphvys+24.125], yearall+j});
  end
  %заводим цикл который будет складывать все температуры за год исключая
  %переменные NaN (l-считает количество NaN)
  for j=1:(2017-yearall)
      entries32WW6(j)=0;l=0;
   for i=1:length(entries_yearsWW6(j).wspd)
      if entries_yearsWW6(j).wspd(i) <=0 || entries_yearsWW6(j).wspd(i) >=0
          tempWW=entries_yearsWW6(j).wspd(i);
      else 
          tempWW=0; l=l+1;
      end
       entries32WW6(j) = entries32WW6(j) + tempWW;
   end
   entries32WW6(j)=entries32WW6(j)/(length(entries_yearsWW6(j).wspd)-l);
  end
  %строим график средних температур за 30 лет 
  %Прямая
  pWW=polyfit(yearWW,entries32WW6,1);
  x1=linspace(yearWW(1),yearWW(end)+7,100);
  y1=polyval(pWW,x1);
  plot3=plot(yearWW,entries32WW6,x1,y1);
  ylabel('Скорость ветра, м/с'); xlabel(' '); 
  set(plot3(1),'linewidth',1.3,'Marker','pentagram','Color','k');
  set(plot3(2),'color','m','linewidth',1.5,'linestyle',':'); 
  grid on;
  hold on;
  legend('Поверхность','линейный тренд','5 км','линейный тренд','10 км','линейный тренд','15 км','линейный тренд','25 км','линейный тренд','location','EastOutside'); legend('boxon');
  
