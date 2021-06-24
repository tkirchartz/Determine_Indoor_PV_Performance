%convert JVAll and JVAlldark 2 Matlab table
% varin=JVAll;
%07.01.2021 version für eingelesene Paper daten zB CuiAdvMat
function convertResults2mat(results,drt,idx,w,Spectra,idxJsc,onlyBest)

%new directory fot txt files
newdir=[drt 'TextFiles\'];
                %check if new dir exist
                if exist(newdir,'dir')==0
                mkdir(newdir)
                end
 
%%%% write Voc FF Jsc messpunkte                
sampleSize=length(results);
varnamesLight={'Jsc', 'Voc', 'FF'};

    for kk=1:sampleSize         %loops samples
    
    headerLight={'\i(J)\-(sc)'	'\i(V)\-(oc)'	'\i(FF)';
             'mA\g(×)cm\+(-2)'  'V'	     '%';
             results(kk).sample results(kk).sample results(kk).sample};
    
    %write cell with real measured Jsc Voc FF of LED     
    lightSize=length(results(kk).Jsc);
    light=NaN(lightSize,3);
        for ii=1:length(varnamesLight)
           light(:,ii)=results(kk).(varnamesLight{ii}); 
        end
        
    cellLight=[headerLight; num2cell(light)];
    %problem writecell doesnt accept ':' in path to save
    if contains(results(kk).sample,':')
        split=regexp(results(kk).sample,':','split');
        name{kk,1}=[split{1} '_' split{2}];
    else
        name{kk,1}=results(kk).sample;
    end
    %save
    writecell(cellLight,[newdir 'JVLEDresults_' name{kk} '.txt'],'Delimiter',',');   %,
    end
    
%%%% write results f
varnamesResults=fields(results(1).data);
varnamesResults=varnamesResults([10:13 15:20]);

l=length(results(1).data);

for kk=1:sampleSize  %loops entry in results (sample)
          matResults=NaN(length(idx),length(varnamesResults)+3);
    if onlyBest==0
        for ll=1:l   %loopsdifferent LEDs
             headerResults={'Material' 'LED' 'wInflection' 'Illuminance \i(E)\-(v)' 'f' 'Irradiance \i(E)\-(e)' 'short circuit current \i(J)\-(sc)' 'open circuit voltage \i(V)\-(oc)' 'fill factor \i(FF)' 'output power density \i(P)\-(out)' 'efficiency \g(\i(h))' 'normalized output power density \i(P)\-(out,norm)' 'normalized efficiency \g(\i(h))\-(,norm)';...		
                            ' '         ' '   'nm'   'Lux' '_' 'µW\g(×)m\+(-2)' 'µA\g(×)cm\+(-2)' 'V' '%' 'µW\g(×)m\+(-2)' '%' '[]' '[]' ;	...	
                            results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample;...
                            results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:}  results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} };
     
            for ii=1:length(varnamesResults)
                matResults(:,ii+3)=results(kk).data(ll).(varnamesResults{ii})(idx);    
            end
            matResults(:,3)=round(results(kk).winfl,1);
            %write material and LED in first two rows
            cellResults=num2cell(matResults);
            cellResults(:,1)={results(kk).sample};
            cellResults(:,2)=results(kk).data(ll).comment;
            
            %concentrate header and mat
            cellResults=[headerResults; cellResults];
            %save
            writecell(cellResults,[newdir 'results_' results(kk).data(ll).comment{:} '_' name{kk} '.txt'],'Delimiter',',');   %,
  
         
        end
    else 
         [~,ll]=max([results(kk).data.eta_f]);
%         for ll=1:l   %loopsdifferent LEDs
             headerResults={'Material' 'LED' 'wInflection' 'Illuminance \i(E)\-(v)' 'f' 'Irradiance \i(E)\-(e)' 'short circuit current \i(J)\-(sc)' 'open circuit voltage \i(V)\-(oc)' 'fill factor \i(FF)' 'output power density \i(P)\-(out)' 'efficiency \g(\i(h))' 'normalized output power density \i(P)\-(out,norm)' 'normalized efficiency \g(\i(h))\-(,norm)';...		
                            ' '         ' '   'nm'   'Lux' '_' 'µW\g(×)m\+(-2)' 'µA\g(×)cm\+(-2)' 'V' '%' 'µW\g(×)m\+(-2)' '%' '[]' '[]' ;	...	
                            results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample results(kk).sample;...
                            results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:}  results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} results(kk).data(ll).comment{:} };
     
            for ii=1:length(varnamesResults)
                matResults(:,ii+3)=results(kk).data(ll).(varnamesResults{ii})(idx);    
            end
            matResults(:,3)=round(results(kk).winfl,1);
            %write material and LED in first two rows
            cellResults=num2cell(matResults);
            cellResults(:,1)={results(kk).sample};
            cellResults(:,2)=results(kk).data(ll).comment;
            
            %concentrate header and mat
            cellResults=[headerResults; cellResults];
            %save
            writecell(cellResults,[newdir 'results_' results(kk).data(ll).comment{:} '_' name{kk} '.txt'],'Delimiter',',');   %,
  
         
%         end
    end
        
       
        
end
 
%write wavelength dependent data
varnamesSpectra=fields(results(1).data);
varnamesSpectra=varnamesSpectra([4:7 14]);

if Spectra==1
    for kk=1:length(idx)          %loops different illuminance levels f
        for ll=1:l   %loopsdifferent LEDs
            headerSpectra={'wavelength \g(l)' 'spectral Irradiance\i( E)\-(e \g(l))'	'spectral Illuminance\i( E)\-(v \g(l))'	'integrated Irradiance'	'integrated Illuminance'	'integrated short circuit current density \i(J)\-(sc)';...
                               'nm'                 'W\g(×)m\+(-2)\g(×)nm\+(-1)'            	'Lux\g(×)nm\+(-1)'                  	'W\g(×)m\+(-2)'     	'Lux'                   'mA\g(×)cm\+(-2)';...
                            results(idxJsc).data(ll).comment{:} results(idxJsc).data(ll).comment{:} results(idxJsc).data(ll).comment{:} results(idxJsc).data(ll).comment{:} results(idxJsc).data(ll).comment{:} results(idxJsc).data(ll).comment{:};...
                            results(idxJsc).data(ll).f(idx(kk)) results(idxJsc).data(ll).f(idx(kk)) results(idxJsc).data(ll).f(idx(kk)) results(idxJsc).data(ll).f(idx(kk)) results(idxJsc).data(ll).f(idx(kk)) results(idxJsc).data(ll).f(idx(kk))...
            };
           matSpectra=NaN(length(w),length(varnamesSpectra));
                for ii=1:length(varnamesSpectra)
                    if ii~=5   %case for 1:4, spIrri, spIllu,integratedIrri and Integrated Illu
                         matSpectra(:,ii)=results(idxJsc).data(ll).(varnamesSpectra{ii}).*results(idxJsc).data(ll).f(idx(kk));    
                    else %case for integrated Jsc
                         matSpectra(:,ii)=results(idxJsc).data(ll).(varnamesSpectra{ii})(idx(kk),:)'; 

                    end
                end
            matSpectra=[w matSpectra];
            cellSpectra= [headerSpectra; num2cell(matSpectra)];

            writecell(cellSpectra,[newdir 'Spectra_' name{idxJsc,1} '_' results(idxJsc).data(ll).comment{:} '_' num2str(round(results(idxJsc).data(ll).Illu_f(idx(kk)))) 'Lux.txt'],'Delimiter',',');   %,



        end
    end
end

        
        
end