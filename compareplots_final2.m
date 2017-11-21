function compare_all_estimators()


clear all;
close all;
rng(2013)
%% Main variables
Totalpts=3000;
n=100;            %ambient dimension;
cycles=1 ;     % average performance over #cycles
percs=[80];  % percentage of outliers
k=1; %inlier subspace dimension
TM = 10;   % Trace of artificial covariance matrices
methodnames={'1','1R1','1R2','2','2R','3','3R'};
methodnames2={'SGD','R-SGD1','R-SGD2','Inc','R-Inc','MD','Robust - MD'};
numberofmethods=length(methodnames);

if(~isdeployed)
  cd(fileparts(which(mfilename))); %change working directory to script directory
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         | CV1 |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   cov = | CV2 |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         | CV3 |
CV1=[1:n];
CV1=CV1/sum(CV1);
cov=CV1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% loop over covariance matrices
for i=[1:size(cov,1)]
    %% loop over nout/nin percentages
    for perc=percs
        nout = round(Totalpts*perc/100);
        nin = Totalpts - nout;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        CV=diag(cov(i,:));
        

        INDsize=100;
        P = zeros(numberofmethods,Totalpts);  % The ith row of P store the principle angles as a function of iterations for the ith method
        pR2=zeros(1,INDsize);

        
        disp('Computing estimators...')
        %% loop over cycles
        for j=1:cycles
 
            [data,SBSP] = datagen2(n,k,nin,nout,CV,TM);
            [initvecs, initvals ] = pca_initialize_random_orthogonal(k,n);
            
            IND=0;
            if any(strcmp(methodnames,'1'))
                
                [p_angles,I,U1]=first(n,k,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles; 
            end
            disp('.')
            
            p=1;
            
            if any(strcmp(methodnames,'1R1'))
                [p_angles,I,U1]=firstR(n,k,data,initvecs,SBSP,p,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
                
            end
            disp('.')
            if any(strcmp(methodnames,'1R2'))
                [p_anglesR2,IR2,U1]=firstR2(n,k,data,initvecs,INDsize,SBSP,p,0);
                pR2=pR2+p_anglesR2;
                IND=IND+1;
                P(IND,:)=P(IND,:)+zeros(1,Totalpts)*4;
            end
            
            disp('.')
            if any(strcmp(methodnames,'2'))
                [p_angles,I,U]=second(n,k,data,initvecs,initvals,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'2R'))
                [p_angles,I,U,S]=secondR(n,k,data,initvecs,initvals,SBSP,0);
                IND=IND+1
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'3'))
                [p_angles,I,bvecs,bvals]=three(n,k,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'3R'))
                [p_angles,I,bvecs]=threeR(n,k,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
                
            end
            disp('.')
           

        end
    
        P=P/cycles;
        pR2 = pR2/cycles;

        filename=sprintf(['results/CV%d_Perc%d_Cyc%d.csv'], i, perc, j )
        csvwrite(filename,P);
        filename=sprintf(['results/PR2_CV%d_Perc%d_Cyc%d.csv'], i, perc, j )
        csvwrite(filename,pR2);
        
    
    end
    
    
end

disp('Plotting')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Plot

 s={'-dk','-hc','-^m','-og','-sb','-vg','-pr'};
 orange =  [ 0.9100 0.4100 0.1700];
 colrs={'k','c','m',orange,'b','g','r'};
 
 
 t = 1:numel(I);
 every = 100;

for CV_index = 1:size(cov,1)
    
    
    %% Loop over inlier/outlier percentages

    for perc_index = 1:length( percs )
        fig=figure;clf;
      
        perc = percs( perc_index )

        filename=sprintf(['results/CV%d_Perc%d_Cyc%d.csv'], CV_index, perc, cycles )
        title_name = sprintf(['outlier percentage = %d%% averaged over %i cycles'],  perc, cycles);
        data = csvread(filename);
        filename=sprintf(['results/PR2_CV%d_Perc%d_Cyc%d.csv'], CV_index, perc, cycles )
        if exist(filename, 'file')
            PR2 = csvread(filename);
        end
        
        data=data';
        for i=1:numberofmethods
            hold on;   
            if strcmp(methodnames{i},'1R2')             
                plot(IR2',PR2',s{i},'MarkerSize',12,'MarkerFaceColor',colrs{i});
            else              
                plot(t([1,100:100:end])',data([1,100:100:end],i),s{i},'MarkerSize',12,'MarkerFaceColor',colrs{i});
            end
            
            hold;
            
        end        
        
        
        legend(methodnames2,'Location','Northeast')
        xlabel( 'Iterations' );
        ylabel( 'Normalized subspace angle' );
        title(title_name,'FontWeight','bold');
    end
    
   
end
            
            
            