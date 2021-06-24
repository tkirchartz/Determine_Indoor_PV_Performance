% programm to determine eta with QE, different LED Spectra and given JV LED
% measurements for FF, Voc
%all results can be found in the structure results

% You can use the matlab code 'EfficiencyCalcuation_Luebke' to calculate effiencies
% for indoor photovoltaics as stated in our paper 'Comparing and Quantifying Indoor
% Performance of Indoor Organic Photovoltaics'.
% First use the dataTemplate origin file to match your data to the template and save it
% as txt files, so matlab should work fine.
% The program is designed to use multiple LEDs and samples. For every sample eta is
% calculated for a range of illuminances and for each LED. When all calcuations are
% done you can use the code convertResults2Mat, to save the findings in txt files for
% origin again (ImportResultsLEDToOrigin.oif is the origin filter you can use).
%We hope everything works fine!
%%
close all
clear all
clc

q=1.6e-19;                       %[As], elementary charge
c=299792458;                     %[m/s], speed of light
h=6.626070040e-34;               %[J*s], Planck's constangt



%find beginning of path for every used computer
%for example: ->C:\Users\luebke\scieboLaptop\ changes
p=mfilename('fullpath');
stelle=strfind(p,'\');
computer=p(1:stelle(3));
Path_def = [fileparts(p) '\'];        %% set the default directory for data
clear stelle

%if this doesnt work type in path an comment section above
%Path_def='C:\Users\luebke\Documents\....\'   %


%%%%%%%% things to set %%%%%

%choose range in lux, either you can just choose one number e.g. 200, to
%calculate the performance at an illuminance of 200 lux, or you chose
%multiple illuminances to calculate the performance for a whole range of
%illuminances, for example [200 500 1000]. Note that low light JV
%measurements must be done close to the desired lux levels
range=[200 500 1000]';        % [lux], chosen range should be close to JscVocFF data from JV measurement to ensure good interpolation, make sure range is a one collum vector!
% range=logspace(log10(200),log10(1E5),110)'    %example for a whole range of illuminancnes from 200 to 1E5 lux in logarithmic scale

%normalization: you can do the efficiency calculation for multiple LED
%spectra. To compare the results and focus on the dependence of the efficiency on different
%LED spectra, it can be handy to normalize the efficiency to a certain LED.
%If you don't want to normalize set normEta=[]; type in as string as indicated in the structure
%results.data.comment eg 'LED 2700 K This work';
normEta=[]; %chose one light source to norm eta and Pout or set to [],

%only compare: to check wether some data in literature credibly reports a
%certain efficiency, you can set onlyCompare=1; then instead of calculating
%the efficiency of a sample for the whole set of given LEDs, only the
%efficiency under the LED of the original paper (the ones you measured the
%JV low light measurements with) is taken into account. This gets a more
%overseeable results structure
onlyCompare=0;%set to one,so only data of orignally used LED is saved (to validate method/check if paper and me gets same results)

%%%%%%% get paths of files %%%%%%%%%%%%%

[FileNamesLED, Path] = uigetfile([Path_def '.txt'] , 'Multiselect', 'on','Select spectral irradiance data of LED spectra');  %

if iscell(FileNamesLED)==0                     %% create a new cell array if just one txt file is chosen (then FileNames is string instead of cell)
    helps=FileNamesLED;
    clear FileNamesLED
    FileNamesLED{1}=helps;
    clear helps
end

[FileNamesQE, ~] = uigetfile([Path_def '.txt'],'Select all Quantum Efficiency data which where used in experiment', 'Multiselect', 'on');
clear data_Si
if iscell(FileNamesQE)==0                     %% create a new cell array if just one txt file is chosen (then FileNames is string instead of cell)
    helps=FileNamesQE;
    clear FileNamesQE
    FileNamesQE{1}=helps;
    clear helps
end

[FileNamesJVLED, ~] = uigetfile([Path_def '*.txt'],'Select JV LED results for Voc(Jsc), FF(Jsc) information', 'Multiselect', 'on');
if iscell(FileNamesJVLED)==0                     %% create a new cell array if just one txt file is chosen (then FileNames is string instead of cell)
    helps=FileNamesJVLED;
    clear FileNamesJVLED
    FileNamesJVLED{1}=helps;
    clear helps
end

%% load files

%%%V(lambda)
B = load([Path 'V(lambda).csv']);   %get data from luminous effiency function
wV=B(:,1);           %get wavelength of V(lambda) table [nm]
V=B(:,2);           %get luminous effiency funtction (from o to 1 dimensionless)
clear B

%%AM 1.5
dataAM=importdata([Path 'astmg173(W_m2nm).xls']);          %get AM1.5 spectrum for checking JscAM
dataAM.w=dataAM.data(1:end,1);              %[nm]; energy array for AM1.5
dataAM.spIrra=dataAM.data(1:end,3);         %[W/m²/nm]; spectral irradiance AM1.5

%%  % load LED data (from different txt files)
%the code is written so that LED data can be in one text file with multiple collumns as 
% w[nm] spectral_Irradiance 1 [W/m²/nm] w2[nm] spectral_Irradiance 2[W/m²/nm]
%or in different text files (please use the the origin dataTemplate to ensure your input
% data has the right appearance) 
%the results are saved in the structure dataLED
n=0;
for ii=1:length(FileNamesLED)                               %loops through different LED text files
    clear A
    
    %get data
    if contains(FileNamesLED{ii}, '.txt')==1
        [A,~,~] = importdata([Path,FileNamesLED{ii}]);      %imports data from LED txt file as structre A
        for jj=1:size(A.data,2)/2                           %loops through different collumns of one text file
            dataLED(n+jj).w=A.data(:,2.*jj-1);                                      %import wavelength in [nm]
            dataLED(n+jj).w=dataLED(n+jj).w(~isnan(dataLED(n+jj).w));               %delete NaNs in wavelength vector
            dataLED(n+jj).spIrra=A.data(:,2.*jj);                                   %imports LED spectral irradiance [W/m²/nm]
            dataLED(n+jj).spIrra=dataLED(n+jj).spIrra(~isnan(dataLED(n+jj).spIrra));%delete NaNs in spectral irradiance vector
            dataLED(n+jj).comment=A.textdata(2.*jj);                                %imports name of LED
        end
        n=length(dataLED);                                   %current length of dataLED structure
    end
    
end

clear A comment 

%plot spectral irradiance data
colors=lines(length(dataLED)+1);
figure(1)
g=gcf;
g.WindowState='maximized';

for ii=1:length(dataLED)
    plot(dataLED(ii).w,dataLED(ii).spIrra./max(dataLED(ii).spIrra),'-','Color',colors(ii,:),'LineWidth',2)
    hold on
end
legend([dataLED.comment])
xlabel('wavelength [nm]')
ylabel('normalized spectral irradiance a.u.')
xlim([370 900])                             %set limits of axis
ylim([0 1])

%for saving figures
newdir=[Path 'Figures\'];           % set directory of Figures

if exist(newdir,'dir')==0           %check if new dir exist
    mkdir(newdir)                   %if directory ...\Path\Figures does not exist, create it
end
saveas(g,[newdir 'LEDSpectra.png']) %save LED spectra as png file
%% interpolare all LED data to same axis
%the results are saved in the structure dataLED

w=linspace(300,1000,(1000-300).*5+1)';     %define axis of wavelength from 300 to 1000 nm in 0.2 steps

for ii=1:length(dataLED)
    
    %interpolate data to same w axis
    
    %problem for double existing values (e.g 0.5 W/m²/nm at 520 and 630 nm), interp1 gets problem
    % the trick is to shift with the x and y axis with a small enough increment cumulatively , here 1E-10:
    uniquex=dataLED(ii).w+(linspace(0,1,length(dataLED(ii).w)).*1E-10)';            %create x axis shifted with an increasing increment
    uniquey=dataLED(ii).spIrra+(linspace(0,1,length(dataLED(ii).spIrra)).*1E-10)';  %create y axis shifted with an increasing increment
    %interpolation of LED data to the w axis
    dataLED(ii).spIrra_int=interp1(uniquex,uniquey,w,'linear'); %interpolation
    dataLED(ii).spIrra_int(isnan(dataLED(ii).spIrra_int))=zeros(sum(isnan(dataLED(ii).spIrra_int)),1);%set Irradiance to zero out of boundaries
    
    %check interpolation with a figure
    figure(2)
    g=gcf;
    plot(dataLED(ii).w,dataLED(ii).spIrra,'bo')
    hold on
    plot(w,dataLED(ii).spIrra_int,'r.')
    
    
    
    clear uniquex uniquey %clear variables for next interpolation
end


title('plot for checking interpolation of LED data');
legend('data','interpolation');
xlabel('wavelength nm')
ylabel('spectral irradiance [W/m²/nm]')
xlim([370 900])
saveas(g,[newdir 'LEDSpectra_interpolated.png'])

%% get QE data
%the code is written so that Qe data can be in one text file with multiple collumns as 
% w[nm] Qe_1 [%] w2[nm] Qe_2[%] ...... 
%or in different text files (please use the the origin dataTemplate to ensure your input
% data has the right appearance) 
%the results are saved in the structure dataQE

%create new directory/folder in  main path if not existent yet
newdir_txt=[Path 'TextFiles\'];
%check if new dir exist
if exist(newdir,'dir')==0
    mkdir(newdir)
end

figure(3)
g=gcf;

%load QE
dataQE=struct;      %preallocation structre for QE measurements
n=0;
for ii=1:length(FileNamesQE)                              %loops different QE files
    data_imp=importdata([Path FileNamesQE{ii}],'\t');      %load quantum efficiency of txt data as structure
    sampleNumber=length(data_imp.textdata)/2;              %get number of samples of one txt file
    for jj=1:sampleNumber                                 %loops different collumns of one Qe text file
        
        dataQE(jj+n).sample=data_imp.colheaders(1,jj.*2);                     %sample name as indicated in the text file
        dataQE(jj+n).w=data_imp.data(:,jj*2-1);                               %[nm], wavelength of Qe measuremnt
        dataQE(jj+n).w=dataQE(jj+n).w(~isnan(dataQE(jj+n).w));               %delete NaNs in wavelength vector
        dataQE(jj+n).QE=data_imp.data(:,jj*2)./100;                           %QE data dimensionless []
        dataQE(jj+n).QE=dataQE(jj+n).QE(~isnan(dataQE(jj+n).QE));            %delete NaNs in Qe vector
        
        uniquex=dataQE(jj+n).w+(linspace(0,1,length(dataQE(jj+n).w)).*1E-10)';  %unique vectors for interpolation, see comment above as the same is done for the LED data
        uniquey=dataQE(jj+n).QE+(linspace(0,1,length(dataQE(jj+n).QE)).*1E-10)';
        %interpolation QE to axis of LED measurement
        dataQE(jj+n).QE_int2=interp1(uniquex,uniquey,w,'linear');               %interpolation of QE to the w axis
        dataQE(jj+n).QE_int2(isnan(dataQE(jj+n).QE_int2))=zeros(sum(isnan(dataQE(jj+n).QE_int2)),1);%set QE to zero out of boundaries
        
        %%find inflection point of QE
        y=dataQE(jj+n).QE_int2;     %QE to calculate with
        nozero=find(y) ;            %exclude zero elements (out of boundary of QE) to exclude pole in gradient
        y=y(nozero);                %QE axis without zero elements         
        wcal=w(nozero);             %w axis without zero elements
        
        [maxi,idxmax]=max(y);      %find maximum of QE and limit vectors to right side of maximum to exclude other inflection points
        ycal=y(idxmax:end); wcal=wcal(idxmax:end);  %QE axis right to maximum and w axis right to maximum
        
        grady=gradient(ycal,wcal); %gradient
        [~,idx]=min(grady);   %find index of minimum of gradient
        
        dataQE(jj+n).winfl=wcal(idx);   %save inflection point in dataQE sturcture
        dataQE(jj+n).idxinfl=find(w==dataQE(jj+n).winfl); %save index of inflection point with original length of axis
        
        %plot for checking interpolation
        plot(dataQE(jj+n).w,dataQE(jj+n).QE,'bo')                                                  %plot original QE data
        hold on
        plot(w,dataQE(jj+n).QE_int2,'r.')                                                          %plot interpolated QE data 
        plot(dataQE(jj+n).winfl,dataQE(jj+n).QE_int2(dataQE(jj+n).idxinfl),'r*','MarkerSize',15)   %plot inflektion point
%         plot(wcal,grady,'g.')                                                                    %plot gradient (to check)
        xlabel('wavelength [nm]');
        ylabel('EQE');
        legend('data','interpolation','inflection point');
        clear uniquex uniquey y ycal wcal grady idxmax maxi idx
    
                %incorperate winfl in EQE txt data (for importing txt file
                %to origin later), only if function 'writecell' is existent
                %(from matlab 2019a on)
                if exist('writecell','file')==2     %check wether the function 'writecell' exsists in your current matlab version
                %create matrix of w, QE, normalzied QE and winflection
                data=[dataQE(jj+n).w dataQE(jj+n).QE.*100 dataQE(jj+n).QE./max(dataQE(jj+n).QE) dataQE(jj+n).winfl.*ones(length(dataQE(jj+n).w),1) ];  
                %create header
                header={'wavelength', 'Quantum efficiency', 'normalized QE', 'inflection point';...
                        'nm' ,'[%]','[ ]', 'nm';...
                      dataQE(jj+n).sample{:}, dataQE(jj+n).sample{:}, dataQE(jj+n).sample{:}, dataQE(jj+n).sample{:} };
                cell=[header; num2cell(data)]; %create cell with header an data
        
                name=strrep(dataQE(jj+n).sample{:},':','_');        %change ':' to '_' in sample name so matlab can save
                writecell(cell ,[newdir_txt 'EQEwinfl_' name '.txt'],'Delimiter',',');  %save the cell/data in TextFiles Folder
                end
    end
    n=length(dataQE);
    clear data_imp data header cell name sampleNumber   
end
xlim([300 1000])
title('plot for checking interpolation of QE')
savefig(g, [newdir 'EQE.fig']);
saveas(g, [newdir 'EQE.png']);

%% interpolate QE to solar spectrum and checking JscAM1.5
%the results are saved in the structure dataQE
for jj=1:length(dataQE)
    
    uniquex=dataQE(jj).w+(linspace(0,1,length(dataQE(jj).w)).*1E-10)';
    uniquey=dataQE(jj).QE+(linspace(0,1,length(dataQE(jj).QE)).*1E-10)';
    %interpolation QE to axis of AM1.5
    dataQE(jj).QE_AM=interp1(uniquex,uniquey,dataAM.w,'linear');
    dataQE(jj).QE_AM(isnan(dataQE(jj).QE_AM))=zeros(sum(isnan(dataQE(jj).QE_AM)),1);%set QE to zero out of boundaries

    dataQE(jj).JscAM= q./h./c*1e-10.*trapz(dataAM.w,dataQE(jj).QE_AM.*dataAM.spIrra.*dataAM.w);         %calculate JscAM1.5 via Jsc=h*c/q*int(Qe*spectralIrradiance*lambda)dlambda

    clear uniquex uniquey
end

%% plot LED Spectra and EQE
colors=lines(length(dataLED)+1);
% close all
for jj=1:length(dataQE)
    figure(3+jj)
    g=gcf;
    g.WindowState='maximized';
    plot(w,dataQE(jj).QE_int2,'k-','LineWidth',2)
    hold on
    for ii=1:length(dataLED)
        plot(dataLED(ii).w,dataLED(ii).spIrra./max(dataLED(ii).spIrra),'-','Color',colors(ii,:),'LineWidth',2)
    end
    legend(['EQE' [dataLED.comment] ])
    title(dataQE(jj).sample)
    xlabel('wavelength nm')
    ylabel('normalized spectral irradiance [a.u] and quantum efficiency [ ].')
    xlim([300 950])
    ylim([0 1])
    
    saveas(g,[newdir 'EQEwithLEDSpectra_' num2str(jj) '.png'])
end


%% calculation of normlazied input power
%normalization of LED data to one or multiple illuminance level specified
%in 'range' as described in section 2.1 of the paper
%the results of the normalization are saved in the structure dataLED

%get Vint for data axis
Vint=interp1(wV,V,w);     %interpolates V to the axis of spectral irradiance
Vint(isnan(Vint))=0;      %sets NaNs to 0

for ii=1:length(dataLED)  %    %loop through different LEDs
    
    %spectral illuminace [Lux/nm]
    dataLED(ii).spIllu_int=dataLED(ii).spIrra_int.*683.*Vint;                   %spectral Illuminance [Lux/nm](spectral_irradiance*683[lm/W]*luminous_Efficiency)
    
    %calculate integrated illuminace and irradiance
    dataLED(ii).integratedIrra=cumtrapz(w,dataLED(ii).spIrra_int);              %[lum/nm] cumulatively
    dataLED(ii).integratedIllu=cumtrapz(w,dataLED(ii).spIllu_int);              %[W/m²/nm] cumulatively
    %calculate Illuminance/Irradiance for reference spectra as they are (without factor)
    dataLED(ii).Pin_Ref=trapz(w,dataLED(ii).spIrra_int);                        %in [W/m²], input power of the reference spectrum,     Pin=integral(spectral_reference_irradiance over lambda)
    dataLED(ii).Illu_Ref=trapz(w,dataLED(ii).spIllu_int);                       %[Lux]    , illuminance of reference spectrum, Illuminance=integral(spectral_reference_illuminance over lambda)
    
    %define Lux scala
    dataLED(ii).Illu_f=range;         %[Lux] range specified in the beginning
    %calculate factor of spectra, so that given Lux is reached
    dataLED(ii).f=dataLED(ii).Illu_f./dataLED(ii).Illu_Ref;                     %[]
    dataLED(ii).Pin=dataLED(ii).f.*dataLED(ii).Pin_Ref.*100 ;                  %µW/cm², Pin=Irradiance_reference*Faktor(Lux)
    % dataLED(ii).Pin=trapz(w,dataLED(ii).f.*dataLED(ii).spIrra_int);           %W/m², alternative calculation
    % dataLED(ii).Illu=trapz(w,dataLED(ii).f.*dataLED(ii).spIllu_int);          %test, should be the defined illumination
    
end


%% load Jsc Voc FF data, which are measured at various intensities with one arbitraty light source
%to double check the illuminances of luxmeters can be typed in the txt
%file, as well as already calcuated PCE (e.g. if data is gathered from
%other papers); to make sure your input data has the right appearance
%please check the dataTemplate origin file
%the data is sorted and saved in a dataJVLED structure

dataJVLED=struct;
nn=0;    %length of structure
for iii=1:length(FileNamesJVLED)
    [A,~,~] = importdata([Path,FileNamesJVLED{iii}],'\t');          %loads in VocFFJsc data from Path
    materials=unique(A.textdata(3:end,1),'stable');                 %find materials names from file
    LEDs=unique(A.textdata(3:end,2),'stable');                  	%find LEDs names from file
    names={'Illu','Jsc','Jsccal','Voc','FF','Pin','Pout','eta'}';   %names of variables
    factor=[1 1 1 1 100 1 1 1]';                                    %factor for changing units; all currents loaded from file as in µA/cm² and all powers in µW/cm²
    
    
    for jj=1:length(materials)             %loops sample
        for kk=1:length(LEDs)              %loops LED
            trues=contains(A.textdata(3:end,1),materials(jj)) & contains(A.textdata(3:end,2),LEDs(kk));  %find indices of right material and LED
            if sum(trues)>0                %only write if combination of material and LED is existent
                dataJVLED(kk+nn).sample=materials(jj);
                dataJVLED(kk+nn).LED=LEDs(kk);
                for ii=1:length(names)      %dynamic writing of data into structure
                    dataJVLED(kk+nn).(names{ii})=A.data(trues,ii).*factor(ii);
                end
            end
        end
        nn=length(dataJVLED);  
    end
end
%plot JscVocFF
figure(111)
g=gcf;
g.WindowState='maximized';
colors1=lines(length(dataJVLED));
for ii=1:length(dataJVLED)
    subplot(1,2,1)
    loglog(dataJVLED(ii).Jsc,dataJVLED(ii).Voc,'o-','Color',colors1(ii,:))
    hold on
    
    subplot(1,2,2)
    loglog(dataJVLED(ii).Jsc,dataJVLED(ii).FF,'x-','Color',colors1(ii,:))
    hold on
end

subplot(1,2,1)
legend([dataJVLED.sample],'Location','SouthEast')
xlabel('Jsc[µA/cm²]');  ylabel('Voc [V]')
subplot(1,2,2)
legend([dataJVLED.sample],'Location','SouthEast')
xlabel('Jsc[µA/cm²]');  ylabel('FF [%]')
   
savefig(g, [newdir 'JscVocFF.fig']);
saveas(g, [newdir 'JscVocFF.png']);

clear nn LEDs materials factor names trues
%% Calulation of Jsc,Voc,FF (Pout) at defined lux level
%all results will be saved in structure results which will be sorted this way:
%the first 'level' will contain one entry for every sample, where
%information on the 'original' measurement, meaning the JV measurment under
%which the solar cell was acutally measured. Here you can find also the
%JscVocFF data. 
%the second 'level' of every sample can be found in results.data (doubleClick); if you
%e.g calculate the efficiency for two different LEDs you will find two rows
%here with all information about the used LED and the calculated
%Jsc_f,Voc_f,FF_f, eta_f... at all lux levels you specified in range.

results=struct;
dels=[];
for jj=1:length(dataJVLED)    %loops through samples of dataJVLED (every sample which was given read in by JscVocFF text files)
    results(jj).sample=dataJVLED(jj).sample{:};                          %save sample name in structure
    results(jj).LED=dataJVLED(jj).LED{:};                                %save LED name in structure, the one under which the sample was actually measured
    results(jj).data=dataLED;                                               %copy LED data from dataLED to the results structure 
    % find sample in QE
    results(jj).idxQE=find(strcmp([dataQE.sample],results(jj).sample));     %find index of current sample in dataQE structure
    results(jj).Jsc=dataJVLED(jj).Jsc;                                   %save measured JscVocFF in results structure
    results(jj).Voc=dataJVLED(jj).Voc;                                   %save measured JscVocFF in results structure
    results(jj).FF=dataJVLED(jj).FF;                                     %save measured JscVocFF in results structure
    
    %calculation
    if isempty(results(jj).idxQE)==0        %only calculate if Qe data of current sample exist in dataQE, idxQE is empty if QE measurement not found
        
        results(jj).winfl=dataQE(results(jj).idxQE).winfl;      %save inflection point in results structure
        
        for ii=1:length(results(jj).data)           %loop through different LEDs, here everything wil be saven in 'second level' of results.data
            % Calculate Jsc @ f, if you calculate for 3 illuminancnes [200
            % 500 1000] you will get vectors with 3 entries in Jsc_f
            results(jj).data(ii).Jsc_f=  q./h./c*1e-10.*trapz(w,dataQE(results(jj).idxQE).QE_int2.*dataLED(ii).f'.*dataLED(ii).spIrra_int.*w)'.*1e3;   %in µA/cm², Jsc=h*c/q*integral(QE*spectral_irradiance(Lux)*lamnbda  over lambda)
            results(jj).data(ii).integratedJsc_f=  q./h./c*1e-10.*cumtrapz(w,dataQE(results(jj).idxQE).QE_int2.*dataLED(ii).f'.*dataLED(ii).spIrra_int.*w)'.*1e3;   %in µA/cm²  %integrated current: for matrices/multiple illuminances: in rows f is constant; in collumns wavelength is constant
            
            results(jj).data(ii).Voc_f=interp1(results(jj).Jsc,results(jj).Voc,results(jj).data(ii).Jsc_f,'pchip');  %[V], interpolation of Voc data to the calculated Jsc from QE measurments
            results(jj).data(ii).FF_f =interp1(results(jj).Jsc,results(jj).FF, results(jj).data(ii).Jsc_f,'pchip');  %[%], interpolation of FF data to the calculated Jsc from QE measurments
            
            %calculate eta
            results(jj).data(ii).Pout_f=results(jj).data(ii).Jsc_f.*results(jj).data(ii).Voc_f.*results(jj).data(ii).FF_f.*0.01;            %[µW/cm²]                         %in W/m²
            results(jj).data(ii).eta_f=results(jj).data(ii).Pout_f./(results(jj).data(ii).Pin).*100;     %in [%]
            
        end  %ii, different LEDs
        
        %normalize Pout and eta
        if isempty(normEta)==0          %if you want to normalzie Pout and eta
            idxnorm=find(contains([results(jj).data.comment],normEta));         %find index of LED Pout and eta should be normalzied to, eg specified in results(1).data(1).comment
            for ii=1:length(results(jj).data)           %loop through different LEDs
                results(jj).data(ii).Pout_norm=results(jj).data(ii).Pout_f./results(jj).data(idxnorm).Pout_f;  %[%], normalization of Pout
                results(jj).data(ii).eta_norm=results(jj).data(ii).eta_f./results(jj).data(idxnorm).eta_f;  %[%], normalization of eta
            end  %ii, different LEDs
        else                            %if you dont want to normalize Pout and eta
            for ii=1:length(results(jj).data)           %loop through different LEDs
                results(jj).data(ii).Pout_norm=ones(size(range));
                results(jj).data(ii).eta_norm=ones(size(range));
            end  %ii, different LEDs
        end
        if onlyCompare==1            %only show LEDs which data has original JscVocFF (for comparison with original data)
            idxLED=find(contains([results(jj).data.comment],results(jj).LED)); %find index of LED in results.data where its the same as in results
            results(jj).data=results(jj).data(idxLED);
        end
    else
        dels=[dels jj];
        disp(['EQE for samples in dataJVLED of index ' num2str(dels) ' not found. The sample in dataJVLED for which EQE was not found is ignored']);
        
    end
    
    
    
end  %end FilenamesQE/samples
results(dels)=[]; %delete entries with no EQE
clear del idxLED idxnorm

%% plotting Jsc(Lux),Voc(Lux),FF(Lux),Pout(Lux),Pin(Lux),eta(Lux)
% plot for illuminance dependent data (makes more sense if you
%calculate for multiple illuminances)

markers='osx+dv^<>ph';                                              %define markers
varnamesY={'Jsc_f','Voc_f','FF_f','Pout_f','Pin','eta_f'};          %define variable names for dynamic plotting
ylabels={'J_{sc} [µA/cm²]' 'V_{oc} [V]' 'Fillfactor [%]'...         %define y labels
    'Pout [µW/m²]' 'Pin [µW/m²]' 'Efficiency [%]' };

for jj=1:length(results)        %loops samples
    for ii=1:length(dataLED)         %loop through different LEDs
        
        figure(111+jj)
        g=gcf;
        g.WindowState='maximized';      %full window mode
        for tt=1:length(varnamesY)      %loops through variables defined in varnamesY (Jsc,Voc,FF,Pout,Pin,eta)
            figure(111+jj);
            subplot(2,3,tt)             %creates subplots
            loglog(results(jj).data(ii).Illu_f,results(jj).data(ii).(varnamesY{tt}), 'Color', colors(ii,:),'MarkerSize',8, 'Marker',markers(ii) ,'LineStyle','none');  %dynamic loglog ploting of data
            hold on
            ylabel(ylabels{tt});
            xlabel('Illuminance [Lux]');
            xlim([range(1)*0.9 range(end)*1.1])  %set limits
            
            if tt==1                            %plot legend in first subplot
                l=legend([results(1).data.comment])%
                l.Position=[0.02 0.1 0.05 0.8]; %[x y l h ], set position and size of legend
                title(results(jj).sample)       %title of figure
            end   
        end  %end tt, varnames  
    end  %ii, different LEDs

%save figures
name=strrep(results(jj).sample,':','_');        %change ':' to '_' in sample name so matlab can save
saveas(g, [newdir 'results_f(Lux)_' name  '.png']);
saveas(g, [newdir 'results_f(Lux)_' name  '.fig']) 
end %end jj, samples

%% plot Jsc(LED),Voc(LED),FF(LED),Pout(LED),Pin(LED),eta(LED)
% show performance dependent of different LED color temperature

%shorten LED names for plot
xlabels={};
for ii=1:length(dataLED)
    if length(dataLED(ii).comment{:})>15
    xlabels=[xlabels dataLED(ii).comment{1}(1:15)];
    else
    xlabels=[xlabels dataLED(ii).comment{1}];  
    end
end
markers='osx+dv^<>ph';                                           %define markes
varnamesY={'Jsc_f','Voc_f','FF_f','Pout_f','Pin','eta_f'};       %define variable names for dynamic plotting
ylabels={'J_{sc} [µA/cm²]' 'V_{oc} [V]' 'Fillfactor [%]'...      %define y labels
    'Pout [µW/m²]' 'Pin [µW/m²]' 'Efficiency [%]' };
idx=1;      %choose index of illuminance you want to plot (if e.g first illuminance in range is 200 lux, for idx=1 all data for 200 lux is plotted)
axlim=[0 length(dataLED)+1];                                     %setting limit of x axis

for jj=1:length(results)        %loops samples
        
        figure(211)
        g=gcf;
        g.WindowState='maximized';      %full window mode
        for tt=1:length(varnamesY)      %loops through variables
            for ii=1:length(dataLED)         %loop through different LEDs
            subplot(2,3,tt)                  %creates subplots
            plot(ii,results(jj).data(ii).(varnamesY{tt})(idx), 'Color', colors(jj,:),'MarkerSize',8, 'Marker',markers(jj) ,'LineStyle','none');  %dynamic plotting of data
            hold on
            end  %end ii, LED
            
            ylabel(ylabels{tt});                                      %set y axis labels
            ax=gca;                                                   %defines current ax
            ax.XTick=1:length(dataLED) ;                              %set used ticks to 1 2 3 4 5 6
            ax.XTickLabelRotation=90;                                 %rotates label 90°
            ax.XTickLabel= xlabels;                                   %set xlabels
            ax.XLim=axlim;                                            %set axis limits
            if tt==6            %set different axlimits for better text in eta figure
            ax.XLim=axlim+0.5;     
            end  
        title(['@' num2str(results(1).data(1).Illu_f(1)) ' Lux'])   
        end  %tt, varnames
        
       subplot(2,3,6)
       text(ii+0.1*ii,results(jj).data(ii).eta_f(idx),results(jj).sample,'Color', colors(jj,:))       %show sample name as textbox next to efficiency 
end

%save figure
saveas(g, [newdir 'results_f(LED).png']);
saveas(g, [newdir 'results_f(LED).fig']) 

%% %saving data from results structure in text data to import in origin
%the function convertResults2mat gives different options to save the
%results structre to txt files, which can then be imported to origin if
%neccessary. For that the matlab function 'writecell' is needed, which is
%only supported for matlab versions >2019a
%all data in results.data (Jsc_f, Voc_f, FF_f eta_f......) for every sample
%and LED as well as the normalzied LED spectra can be saved
%you can specifiy some things that could be handy

%if you only want to save the performance for the best performance LED set
%only best=1; if you want to save the normalzied spectral illuminances/irradiances set
onlyBest=0; %set to one, so only data of best effiency is saved

%saveSpectra=1, then spectral illuminances/irradiances and integrated Jsc are saved;
%if you dont save any spectra (saveSpectra=0) then saveIllu and idxJsc are not important
%- you can define at which lux levels the spectra should be saved by setting saveIllu to the desired indicess of range
%  whole illuminances: saveIllu=1:length(range)  or just some saveIllu=[1 2 3]
%  depending on how many illuminance levels were calculated, saving all data might take a while
%- you can define from which sample the integrated Jsc should be saved by setting idxJsc to the index of results structure;

saveIllu=1;  %set indices for illuminances which should be saved, 
saveSpectra=1;             %set to 1 for saving LED spectra as txt data for loading into origin
idxJsc=1;                  %chose index of results, for which wavelength dependent data should be saved

if exist('writecell')==2
convertResults2mat(results,Path,saveIllu,w,saveSpectra,idxJsc,onlyBest)
end

% close all
% save([Path 'workspace'])





