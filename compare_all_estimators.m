function compare_all_estimators(D,d,Totalpts,perc_outliers,num_iters,seed)
% compare_all_estimators(D,d,Totalpts,perc_outliers,num_iters,seed)
%
% Generate figures like 1-3 in paper
%
% D - ambient dimension
% d - subspace dimension
% Totalpts - total number of iterations/points
% perc_outliers - percentage of points which are outliers (can be vector)
% num_iters   - number of times to run experiment
% seed        - choose seed for reproducibility
%
%
%Example: compare_all_estimators(100,10,3000,45,1,2013)


rng(seed)
%% Main variables
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
CV1=[1:D];
CV1=CV1/sum(CV1);
cov=CV1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% loop over covariance matrices
for i=[1:size(cov,1)]
    %% loop over nout/nin percentages
    for perc=perc_outliers
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
        for j=1:num_iters
 
            [data,SBSP] = datagen2(D,d,nin,nout,CV,TM);
            [initvecs, initvals ] = pca_initialize_random_orthogonal(d,D);
            
            IND=0;
            if any(strcmp(methodnames,'1'))
                
                [p_angles,I,U1]=first(D,d,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles; 
            end
            disp('.')
            
            p=1;
            
            if any(strcmp(methodnames,'1R1'))
                [p_angles,I,U1]=firstR(D,d,data,initvecs,SBSP,p,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
                
            end
            disp('.')
            if any(strcmp(methodnames,'1R2'))
                [p_anglesR2,IR2,U1]=firstR2(D,d,data,initvecs,INDsize,SBSP,p,0);
                pR2=pR2+p_anglesR2;
                IND=IND+1;
                P(IND,:)=P(IND,:)+zeros(1,Totalpts)*4;
            end
            
            disp('.')
            if any(strcmp(methodnames,'2'))
                [p_angles,I,U]=second(D,d,data,initvecs,initvals,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'2R'))
                [p_angles,I,U,S]=secondR(D,d,data,initvecs,initvals,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'3'))
                [p_angles,I,bvecs,bvals]=three(D,d,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
            end
            disp('.')
            if any(strcmp(methodnames,'3R'))
                [p_angles,I,bvecs]=threeR(D,d,data,initvecs,SBSP,0);
                IND=IND+1;
                P(IND,:)=P(IND,:)+p_angles;
                
            end
            disp('.')
           

        end
    
        P=P/num_iters;
        pR2 = pR2/num_iters;

        filename=sprintf(['results/CV%d_Perc%d_Cyc%d.csv'], i, perc, j );
        csvwrite(filename,P);
        filename=sprintf(['results/PR2_CV%d_Perc%d_Cyc%d.csv'], i, perc, j );
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

    for perc_index = 1:length( perc_outliers )
        fig=figure;clf;
      
        perc = perc_outliers( perc_index );

        filename=sprintf(['results/CV%d_Perc%d_Cyc%d.csv'], CV_index, perc, num_iters );
        title_name = sprintf(['outlier percentage = %d%%, averaged over %i cycle(s)'],  perc, num_iters);
        data = csvread(filename);
        filename=sprintf(['results/PR2_CV%d_Perc%d_Cyc%d.csv'], CV_index, perc, num_iters );
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
            
            
            