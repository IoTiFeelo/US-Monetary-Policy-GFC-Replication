function StateSpaceP_new=EM_DHDFM_SS_block(X,modelSpec,StateSpaceP_old)

% EM algorithm for unrestricted dynamic factor model
%
% X_{t}=\Lambda_{0} F_{t} + ... + \Lambda_{h} F_{t-h} + v_{t};  
% F_{t}=A1 F_{t-1} + ... + Ap F_{t-p} + u_{t};
% v_{t}=D1 v_{t-1} + ... + Ds v_{t-s} + e_{t};
% F_{t}=[f_{1t}',...,f_{rt}']';
% u_{t}~iid(0,\Omega); e~iid(0,\Psi); (\Psi diagonal)
% X=[Nx1]; \Lambda_{j}=[Nxr]; F=[rx1]; e=[Nx1]; u=[rx1]; Ai=[rxr]; E(Fe')=0
% active restriction on quarterly loadings (nR=#of factor lags involved)
%
% SSF:
% X_{t}=C S_{t} + e_{t}
% S_{t}=A S_{t-1} + u_{t}
% e_{t}~N(0,R); u_{t}~N(0,Q); P_{t}=var(S_{t});
% S_{t}=[F_{t},F_{t-1},...,F_{t-max(p-1,nR)},... [r x (max(max(p,h),nR)+1)]
%        vM_{t},vM_{t-1},...,vM_{t-s+1},... monthly idio [nMxs]
%        vQ_{t},vQ_{t-1},...,vQ_{t-max(s-1,nR)}]'; quarterly idio [nQ x max(s,nR+1)]
% nS=nSf+nSiM+nSiQ; %number of states
% nSf=sum(r.*(max(max(p,h),nR)+1)); nSiM=nM*s; nSiQ=nQ*max(s,nR+1); 
% X=[Nx1]; C=[NxnS]; S=[nSx1]; e=[Nx1]; R=[NxN]; 
% A=[nSxnS]; u=[nSx1]; Q=[nSxnS]; 

% INPUT:
% X=[TxN] matrix of data
% modelSpec=structure
% r=number of factors
% p=number of lags in factor VAR
% h=number of factor lags in observation eqn
% s=number of lags in idio dynamic

% OUTPUT: structure containing parameters estimates (SSF)
% S=smoothed states
% P=smoothed states variance
% Xhat=smoothed data
% C,R,A,Q=system matrices


% built on Banbura & Modugno (2010)
% miranda 2014 smirandaagrippino@london.edu
%--------------------------------------------------------------------------



% initialize
modelSpec.s     =1;
modelSpec.h     =zeros(size(X,2),1);


T    =size(X,1); 
mX   =nanmean(X); 
vX   =nanstd(X); 

xNaN =bsxfun(@minus,X,mX); 
xNaN =bsxfun(@rdivide,xNaN,vX);


if ~isempty(StateSpaceP_old)
    
    SystemMatrices_init.A   =StateSpaceP_old.A; 
    SystemMatrices_init.C   =StateSpaceP_old.C;
    SystemMatrices_init.Q   =StateSpaceP_old.Q; 
    SystemMatrices_init.R   =StateSpaceP_old.R; 
    SystemMatrices_init.S   =StateSpaceP_old.S; 
    SystemMatrices_init.S   =StateSpaceP_old.P;
    
else
    
    SystemMatrices_init     =initializeEM(xNaN,modelSpec);
end

%--------------------------------------------------------------------------
%EM algorithm

% set parameters
nmaxit          =modelSpec.max_iter; 
dLLthreshold    =1E-3;  %for likelihood convergence

it=0; 

%loglikelihood
llconverged     =0; 
loglklhd_old    =-inf; 
loglklhd_store  =NaN(nmaxit,1);

%remove leading and ending nans for the estimation
optNaN.method   =3;                
optNaN.k        =3;

y_est   =remNaNs_spline(xNaN,optNaN);


%EM
while it<=nmaxit && llconverged==0

    [SystemMatrices_end,loglklhd]=EMalgorithm(y_est,modelSpec,SystemMatrices_init);
    

    %update matrices
    SystemMatrices_init.S =SystemMatrices_end.S;
    SystemMatrices_init.P =SystemMatrices_end.P;
    
    SystemMatrices_init.C =SystemMatrices_end.C;
    SystemMatrices_init.R =SystemMatrices_end.R;
    SystemMatrices_init.A =SystemMatrices_end.A;
    SystemMatrices_init.Q =SystemMatrices_end.Q;
    
    
        
    %check likelihood convergence
    llconverged =EMconverged(loglklhd,loglklhd_old,dLLthreshold,1,it);

    loglklhd_old        =loglklhd; 
    loglklhd_store(it+1)=loglklhd;
    
    it=it+1;
    
    if it==nmaxit && ~llconverged
        
        fprintf(1,'EM DID NOT CONVERGE: maximum number of iterations reached!\n');
        fprintf(1,'Iterations allowed: %i; Tolerance Level %1.2e \n',nmaxit,dLLthreshold);
    end
    
end

if llconverged
    
    fprintf(1,'EM converged at iteration %i of %i\n',it,nmaxit);
    fprintf(1,'Tolerance Level %1.2e \n',dLLthreshold);
end


% last run of kalman filter & smoother 
%@ iteration j states are estimated using system matrices @ iteration (j-1))
KalmanFilterOutput      =kalman_filter(xNaN,SystemMatrices_end);
KalmanSmootherOutput    =kalman_smoother(KalmanFilterOutput);

%load results
StateSpaceP_new.S       =KalmanSmootherOutput.S_smooth(:,1); %S_{0|0}
StateSpaceP_new.P       =SystemMatrices_end.P; %P_{0|0}

StateSpaceP_new.C       =SystemMatrices_end.C;
StateSpaceP_new.R       =SystemMatrices_end.R;
StateSpaceP_new.A       =SystemMatrices_end.A;
StateSpaceP_new.Q       =SystemMatrices_end.Q;


StateSpaceP_new.States  =KalmanSmootherOutput.S_smooth(:,2:end)'; %states
StateSpaceP_new.Xhat    =((KalmanSmootherOutput.S_smooth(:,2:end)'*...
                         SystemMatrices_end.C').*kron(vX,ones(T,1)))+kron(mX,ones(T,1)); %estimated X

StateSpaceP_new.mX      =mX; 
StateSpaceP_new.vX      =vX; 





% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %
%                              child fcts                                 %
% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %




function [SystemMatrices_end,loglklhd]=EMalgorithm(x,modelSpec,SystemMatrices_init)


[T,N]   =size(x);


%unpack model specification
r       =modelSpec.r;               %number of factors
p       =modelSpec.p;               %number of lags in factor VAR
s       =modelSpec.s;               %number of lags in idio dynamic

%
blocks  =modelSpec.blocks;          %block restrictions



nM      =modelSpec.nM;              %number of monthly variables
nB      =size(blocks,2);            %number of blocks


% *     *     *     *   E X P E C T A T I O N   *     *     *     *

KalmanFilterOutput          =kalman_filter(x,SystemMatrices_init);
KalmanSmootherOutput        =kalman_smoother(KalmanFilterOutput);


%smoothed states and their variance 
S_smooth                    =KalmanSmootherOutput.S_smooth;
P_smooth                    =KalmanSmootherOutput.P_smooth;
PP_smooth                   =KalmanSmootherOutput.PP_smooth;


%log likelihood
loglklhd                    =KalmanFilterOutput.loglklhd;



% *     *     *     *   M A X I M I Z A T I O N    *     *     *     *
%[Shumway&Stoffer(2000)]

nSf     =sum(r.*p); 
nSiM    =nM*s; 
nS      =nSf+nSiM; %number of states


%coefficients matrices
C_end   =zeros(N,nS); 
R_end   =zeros(N,1); 
A_end   =zeros(nS,nS); 
Q_end   =zeros(nS,nS);


P       =zeros(nS,nS); 
Pl      =zeros(nS,nS); 
PPl     =zeros(nS,nS);

for t=1:T
    
    P   =P+P_smooth{t+1};       %sum(var(S_{t})
    Pl  =Pl+P_smooth{t};        %sum(var(S_{t-1})
    PPl =PPl+PP_smooth{t+1};    %sum(cov(S_{t}S_{t-1})
end


% *~~~*~~~*~~~*~~~*~~~*~~~* STATE EQUATION *~~~*~~~*~~~*~~~*~~~*~~~*
% *~~~* FACTORS *~~~*

nSfi=cumsum([1 r.*p]);

for b=1:nB

    rb      =r(b); %# of factors per block
    
    
    % sufficient statistics
    E_FF    =S_smooth(nSfi(b):nSfi(b)+rb*p-1,2:end)*S_smooth(nSfi(b):nSfi(b)+rb*p-1,2:end)'+...
             P(nSfi(b):nSfi(b)+rb*p-1,nSfi(b):nSfi(b)+rb*p-1); %E(F_{t}F_{t}')

    E_FlFl  =S_smooth(nSfi(b):nSfi(b)+rb*p-1,1:end-1)*S_smooth(nSfi(b):nSfi(b)+rb*p-1,1:end-1)'+...
             Pl(nSfi(b):nSfi(b)+rb*p-1,nSfi(b):nSfi(b)+rb*p-1); %E(F_{t-1}F_{t-1}')

    E_FFl   =S_smooth(nSfi(b):nSfi(b)+rb*p-1,2:end)*S_smooth(nSfi(b):nSfi(b)+rb*p-1,1:end-1)'+...
             PPl(nSfi(b):nSfi(b)+rb*p-1,nSfi(b):nSfi(b)+rb*p-1); %E(F_{t}F_{t-1}')
    
         
    %
    tempA                       =E_FFl/E_FlFl; 
    blockA                      =zeros(rb*p,rb*p);
    blockA(1:rb,1:rb*p)         =tempA(1:rb,1:rb*p); 
    blockA(rb+1:end,1:end-rb)   =eye(rb*(p-1));
    
    A_end(nSfi(b):nSfi(b+1)-1,nSfi(b):nSfi(b+1)-1) = blockA;
    
    
    %
    tempQ                       =(1/T)*(E_FF-tempA*E_FFl'); 
    blockQ                      =zeros(rb*p,rb*p);
    blockQ(1:rb,1:rb)           =tempQ(1:rb,1:rb);
    
    Q_end(nSfi(b):nSfi(b+1)-1,nSfi(b):nSfi(b+1)-1) = blockQ;
    
end


% *~~~* IDIOSYNCRATIC *~~~*
%idioM
for i=1:nM
    
    % sufficient statistics
    E_eMeM   =S_smooth(nSf+1+(i-1)*s:nSf+i*s,2:end)*S_smooth(nSf+1+(i-1)*s:nSf+i*s,2:end)'+...
              P(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s); %E(eM_{t}eM_{t}')

    E_eMleMl =S_smooth(nSf+1+(i-1)*s:nSf+i*s,1:end-1)*S_smooth(nSf+1+(i-1)*s:nSf+i*s,1:end-1)'+...
              Pl(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s); %E(eM_{t-1}eM_{t-1}')
    
    E_eMeMl  =S_smooth(nSf+1+(i-1)*s:nSf+i*s,2:end)*S_smooth(nSf+1+(i-1)*s:nSf+i*s,1:end-1)'+...
              PPl(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s); %E(eM_{t}eM_{t-1}')
    
          
    %
    tempA                   =E_eMeMl/E_eMleMl; 
    blockA                  =zeros(s,s); 
    blockA(1,:)             =tempA(1,:); 
    blockA(2:end,1:end-1)   =eye(s-1);
    
    A_end(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s) = blockA;
   
        
    %
    tempQ                   =(1/T)*(E_eMeM-tempA*E_eMeMl'); 
    blockQ                  =zeros(s,s);
    blockQ(1,1)             =tempQ(1,1); 
    
    Q_end(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s) = blockQ;
    
end




% *~~~*~~~*~~~*~~~*~~~* MEASUREMENT EQUATION *~~~*~~~*~~~*~~~*~~~*
% define superblocks as combination of factors and group variables
% according to which one they load on distinguishing bwn M and Q

superB      =unique(blocks,'rows'); 
nSB         =size(superB,1);

%factors and their lags
MloadPerSB  =zeros(nSB,nSf); 


nSfi    =cumsum([1 r.*p]); 

for b=1:nB
    
    %position of monthly loadings in super block
    MloadPerSB(:,nSfi(b):nSfi(b)+r(b)-1)           =repmat(superB(:,b),1,r(b));
    
end

%idiosyncratic
iMloadPerSB     =zeros(N,nS); 

iMloadPerSB(1:nM,nSf+1:nSf+nSiM)    =kron(eye(nM),[1 zeros(1,s-1)]); 


%build identifiers -- these will pick the right columns of the states
%vector depending on which variables load on which states in each block
MloadPerSB      =logical(MloadPerSB); 
iMloadPerSB     =logical(iMloadPerSB); 


for sb=1:nSB
    

    %find all variables that load on all factors in super block
    selectV     =find(ismember(blocks,superB(sb,:),'rows'));
    
    MinSB       =selectV(selectV<=nM); 
    nMSB        =numel(MinSB); 
    
    
    %vec(C)=[sum(kron(E(FF'),W))]^{-1}[vec(sum(WyE(F')))-vec(sum(WE(F*e')))]
    %vec(C)=[MC1]^{-1}[MC2]
    % *~~~* MONTHLY LOADINGS *~~~*
    if~isempty(MinSB)
        
        MC1     =zeros(nMSB*sum(MloadPerSB(sb,:)),nMSB*sum(MloadPerSB(sb,:))); 
        MC2     =zeros(nMSB,sum(MloadPerSB(sb,:)));
        
        for t=1:T
            
            %track missing values
            Wx =x(t,MinSB)';     W =diag(~isnan(Wx));     Wx(isnan(Wx)) =0;
            
            
            MC1 =MC1 + kron(S_smooth(MloadPerSB(sb,:),t+1)*S_smooth(MloadPerSB(sb,:),t+1)'+...
                 P_smooth{t+1}(MloadPerSB(sb,:),MloadPerSB(sb,:)),W);
            
            %
            MC2 =MC2 + Wx*S_smooth(MloadPerSB(sb,:),t+1)'-...
                 W*(S_smooth(any(iMloadPerSB(MinSB,:),1),t+1)*S_smooth(MloadPerSB(sb,:),t+1)'+...
                 P_smooth{t+1}(any(iMloadPerSB(MinSB,:),1),MloadPerSB(sb,:)));
        end
        
        %update C matrix for monthly variables
        vecC    =MC1\MC2(:); 
        tempC   =reshape(vecC,nMSB,sum(MloadPerSB(sb,:)));
        
        %update C
        C_end(MinSB,MloadPerSB(sb,:))=tempC;

    end
    
end


C_end(1:nM,nSf+1:nSf+nSiM)              =kron(eye(nM),[1 zeros(1,s-1)]);


%update R (small number)
R_end   =diag(1e-04*ones(N,1));


%update P
P_end   =(eye(nS^2) - kron(A_end,A_end))\Q_end(:);
P_end   =reshape(P_end,nS,nS);


%load structure with results
SystemMatrices_end.S =S_smooth(:,1); 
SystemMatrices_end.P =P_end;
SystemMatrices_end.C =C_end; 
SystemMatrices_end.R =R_end;
SystemMatrices_end.A =A_end; 
SystemMatrices_end.Q =Q_end;

% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %




function [converged,decrease] = EMconverged(loglik,previous_loglik,threshold,check_increased,iteration)
%
% We have converged if the slope of the log-likelihood function falls below 'threshold',
% i.e., |f(t) - f(t-1)| / avg < threshold,
% where avg = (|f(t)| + |f(t-1)|)/2 and f(t) is log lik at iteration t.
% 'threshold' defaults to 1e-4.
%
% This stopping criterion is from Numerical Recipes in C p423
%
% If we are doing MAP estimation (using priors), the likelihood can decrase,
% even though the mode of the posterior is increasing.

if nargin<3, threshold=1e-4;     end
if nargin<4, check_increased=1;  end
if nargin<5, iteration=1;        end

converged   =false; 
decrease    =false;

if check_increased
    if loglik-previous_loglik < -1e-3 % allow for a little imprecision

        fprintf(1,'(!) likelihood decreased from %6.4f to %6.4f at iteration %i (!)\n',previous_loglik,loglik,iteration);
        decrease=true;
        
    end
end

delta_loglik =abs(loglik-previous_loglik);
avg_loglik   =(abs(loglik)+abs(previous_loglik)+eps)/2;

if (delta_loglik/avg_loglik)<threshold; converged=true; end

% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %




function SystemMatrices_init = initializeEM(x,modelSpec,optNaN)

r       =modelSpec.r; %number of factors
p       =modelSpec.p; %number of lags in factor VAR
s       =modelSpec.s; %number of lags in idio dynamic

nM      =modelSpec.nM; 


%
blocks  =modelSpec.blocks; %blocks restrictions
nB      =size(blocks,2);



%remove nans for initialization with principal components
%replace with MA(k) after deleting rows with leading or ending nans
optNaN.method   =2;                 
optNaN.k        =3;
[xSpline,~]     =remNaNs_spline(x,optNaN); 

[T,N]           =size(xSpline);



nSf             =sum(r.*p); 
nSiM            =nM*s; 
nS              =nSf+nSiM;     %number of states

C_init          =zeros(N,nS); 
A_init          =zeros(nS,nS); 
Q_init          =zeros(nS,nS); 



% initialize state space w\ principal components & least squares
tempX           =xSpline;           %data without nans

tempF           =NaN(T,nSf); 


nSfi            =cumsum([1 r.*p]);

for b=1:nB
    
    % *~~~*~~~*~~~*~~~*~~~* MEASUREMENT EQUATION *~~~*~~~*~~~*~~~*~~~* 
    rb      =r(b); %# of factors per block
    nlF     =1; %# of effective lags
    
    btemp   =find(blocks(:,b)); 
    MinB    =btemp(btemp<=nM); 
    

    %factors are initialized on the monthly variables
    vX      =cov(tempX(:,MinB)); 
    [V,~]   =eigs(vX,rb);
    F       =tempX(:,MinB)*V; 
    F       =F(:,1:rb); %principal components for factors
    
    
    lagF_Meqn=nan(T-nlF,rb*(nlF+1)); %[F_{t},F_{t-1},...,F_{t-h}]';
    for j=1:nlF+1
        
        lagF_Meqn(:,rb*(j-1)+1:rb*j)=F((nlF+1)-j+1:end-(j-1),:);
    end
    tempC =lagF_Meqn(:,1:rb)\tempX(nlF+1:end,MinB); %loadingsM

    C_init(MinB,nSfi(b):nSfi(b)+rb-1) =tempC';
    %

    
    %orthogonalize for next block extraction
    tempF(:,nSfi(b):nSfi(b)+rb*(nlF+1)-1)=[zeros(size(tempX,1)-size(lagF_Meqn,1),rb*(nlF+1));lagF_Meqn];
    
    %`data' orthogonal to factors in preceding block
    tempX =tempX-tempF(:,nSfi(b):nSfi(b)+rb*(nlF+1)-1)*C_init(:,nSfi(b):nSfi(b)+rb*(nlF+1)-1)';

    
    
    % *~~~*~~~*~~~*~~~*~~~*~~~* STATE EQUATION *~~~*~~~*~~~*~~~*~~~*~~~* 
    lagF_Seqn =nan(T-(p+1),rb*(p+1)); %[F_{t},F_{t-1},...,F_{t-h}]';
    for j=1:p+1

        lagF_Seqn(:,rb*(j-1)+1:rb*j) =F(((p+1)+1)-j:end-j,:);
    end    
    
    
%     lagF_Seqn   =lagF_Meqn(:,1:rb*(p+1)); %[F_{t},F_{t-1},...,F_{t-p}]';
    tempA       =lagF_Seqn(:,rb+1:end)\lagF_Seqn(:,1:rb);
    u           =lagF_Seqn(:,1:rb)-lagF_Seqn(:,rb+1:end)*tempA; 
    tempQ       =cov(u);
    
    %
    blockA                      =zeros(rb*(nlF+1),rb*(nlF+1));
    blockA(1:rb,1:rb*p)         =tempA'; 
    blockA(rb+1:end,1:end-rb)   =eye(rb*nlF);
    
    A_init(nSfi(b):nSfi(b+1)-1,nSfi(b):nSfi(b+1)-1) = blockA;
    
    
    blockQ                      =zeros(rb*(nlF+1),rb*(nlF+1));
    blockQ(1:rb,1:rb)           =tempQ;
    
    Q_init(nSfi(b):nSfi(b+1)-1,nSfi(b):nSfi(b+1)-1) = blockQ;
end


% *~~~*~~~*~~~*~~~*~~~* MEASUREMENT EQUATION *~~~*~~~*~~~*~~~*~~~*
C_init(1:nM,nSf+1:nSf+nSiM)=kron(eye(nM),[1 zeros(1,s-1)]);



R_init=diag(1e-04.*ones(N,1)); %can be improved collecting the states



% *~~~*~~~*~~~*~~~*~~~*~~~* STATE EQUATION *~~~*~~~*~~~*~~~*~~~*~~~* 
E=tempX-tempF*C_init(:,1:nSf)'; %orthogonal to all the factors

%idioM
for i=1:nM
    
    iE =E(:,i); 
    
    lagE_Seqn=nan(T-s,s+1); %[e1_{t},...,e1_{t-s},...,enM_{t},...,enM_{t-s}]';
    for j=1:s+1
        
        lagE_Seqn(:,j)=iE((s+1)-j+1:end-(j-1),:);
    end
    
    tempA   =(lagE_Seqn(:,2:end)\lagE_Seqn(:,1));
    uM      =lagE_Seqn(:,1)-lagE_Seqn(:,2:end)*tempA; 
    tempQ   =cov(uM);
    
    %
    blockA                  =zeros(s,s); 
    blockA(1,:)             =tempA'; 
    blockA(2:end,1:end-1)   =eye(s-1);
    
    A_init(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s) = blockA;
    
    
    blockQ                  =zeros(s,s); 
    blockQ(1,1)             =tempQ; 
    
    Q_init(nSf+1+(i-1)*s:nSf+i*s,nSf+1+(i-1)*s:nSf+i*s) = blockQ;
end



%-------------------------------------------------------------------------%
    
%states & variance
S_init  =zeros(nS,1);
vecP    =(eye(nS^2)-kron(A_init,A_init))\Q_init(:);
P_init  =reshape(vecP,nS,nS); %var(S_init)


%load structure
SystemMatrices_init.S =S_init; 
SystemMatrices_init.P =P_init;

SystemMatrices_init.C =C_init; 
SystemMatrices_init.R =R_init;
SystemMatrices_init.A =A_init; 
SystemMatrices_init.Q =Q_init;

% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %




function [KalmanFilterOutput,SystemMatrices] = kalman_filter(x,SystemMatrices)
% model in state-space form:
% X_{t}=C S_{t} + e_{t}
% S_{t}=A S_{t-1} + u_{t}
% e_{t}~N(0,R); u_{t}~N(0,Q); P_{t}=var(S_{t});

% INPUT: 
% S_init=vector of states S_{0|0}
% P_init=variance of the states P_{0|0}
% C,R,A,Q=system matrices

% OUTPUT:
% S_filter=S_{t|t} t=0,..,T [contains also initial value]
% P_filter=P_{t|t} t=0,..,T [contains also initial value]
% S_forecast=S_{t|t-1} t=1,...,T
% P_forecast=P_{t|t-1} t=1,...,T
% K_{t}=P_{t|t-1}C'[CP_{t|t-1}C'+R]^{-1} t=1,...,T %gain
% J_{t-1}=P_{t-1|t-1}A'[P_{t|t-1}]^{-1} t=1,...,T
% C,R,A,Q=system matrices
% loglklhd=log likelihood
% includes modification for missing values


% miranda 2014 smirandaagrippino@london.edu
%--------------------------------------------------------------------------


S_init  =SystemMatrices.S; 
P_init  =SystemMatrices.P;

C       =SystemMatrices.C; 
R       =SystemMatrices.R;
A       =SystemMatrices.A; 
Q       =SystemMatrices.Q;

T       =size(x,1); 
nS      =length(S_init); %number of states


S_filter        =nan(nS,T+1); 
S_filter(:,1)   =S_init;
S_forecast      =nan(nS,T);

P_filter        =cell(T+1,1); 
P_filter{1}     =P_init;
P_forecast      =cell(T,1); 

K               =cell(T,1); 
J               =cell(T,1);


loglklhd=0;

for t=1:T
    
    %prediction S_{t|t-1} & P_{t|t-1}
    S_forecast(:,t) =A*S_filter(:,t);
    P_forecast{t}   =A*P_filter{t}*A'+Q; 
    P_forecast{t}   =.5*(P_forecast{t}+P_forecast{t}'); %ensures +ve def
    
    
    %remove rows with missing data
    keeprows        =~isnan(x(t,:)'); 
    
    xt              =x(t,keeprows)'; 
    Ct              =C(keeprows,:); 
    Rt              =R(keeprows,keeprows);
        
    B               =Ct*P_forecast{t}; % B_{t|t-1}=E[(X_{t}-X_{t|t-1})(S_{t}-S_{t|t-1})']
    H               =Ct*P_forecast{t}*Ct'+Rt; % H_{t|t-1}=E[(X_{t}-X_{t|t-1})(X_{t}-X_{t|t-1})']
    K{t}            =B'/H; %gain
    J{t}            =P_filter{t}*A'*pinv(P_forecast{t}); %needed for smoother

    %loglikelihood @ each t
    ll              =log(det(inv(H)))-(xt-Ct*S_forecast(:,t))'/H*(xt-Ct*S_forecast(:,t)); 
    loglklhd        =loglklhd+(1/2)*ll;
    
    %update S_{t|t} & P_{t|t}
    S_filter(:,t+1) =S_forecast(:,t)+K{t}*(xt-Ct*S_forecast(:,t));
    P_filter{t+1}   =P_forecast{t}-K{t}*B; 
    P_filter{t+1}   =.5*(P_filter{t+1}+P_filter{t+1}');
end

SystemMatrices.CbigT=Ct;



%load structure
KalmanFilterOutput.S_filter         =S_filter; 
KalmanFilterOutput.P_filter         =P_filter; 
KalmanFilterOutput.S_forecast       =S_forecast; 
KalmanFilterOutput.P_forecast       =P_forecast;

KalmanFilterOutput.SystemMatrices   =SystemMatrices;

KalmanFilterOutput.K                =K; 
KalmanFilterOutput.J                =J; 

KalmanFilterOutput.loglklhd         =loglklhd;

% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %



function KalmanSmootherOutput = kalman_smoother(KalmanFilterOutput)
% model in state-space form:
% X_{t}=C S_{t} + e_{t}
% S_{t}=A S_{t-1} + u_{t}
% e_{t}~N(0,R); u_{t}~N(0,Q); P_{t}=var(S_{t});

% INPUT: (=kalman_filter output)
% S_filter=S_{t|t} t=0,..,T [contains also initial value]
% P_filter=P_{t|t} t=0,..,T [contains also initial value]
% S_forecast=S_{t|t-1} t=1,...,T
% P_forecast=P_{t|t-1} t=1,...,T
% K_{t}=P_{t|t-1}C'[CP_{t|t-1}C'+R]^{-1} t=1,...,T %gain
% J_{t-1}=P_{t-1|t-1}A'[P_{t|t-1}]^{-1} t=1,...,T
% C,R,A,Q=system matrices

% OUTPUT:
% S_smooth=S_{t|T} t=1,...,T
% P_smooth=P_{t|T}=E[(S_{t}-S_{t|T})(S_{t}-S_{t|T})']
% PP_smooth=P_{t,t-1|T}=E[(S_{t}-S_{t|T})(S_{t-1}-S_{t-1|T})']
% C,R,A,Q=system matrices


% miranda 2014 smirandaagrippino@london.edu
%--------------------------------------------------------------------------



%unpack structures
SystemMatrices  =KalmanFilterOutput.SystemMatrices; 
C               =SystemMatrices.CbigT; 
A               =SystemMatrices.A; 

S_filter        =KalmanFilterOutput.S_filter; 
P_filter        =KalmanFilterOutput.P_filter; 
S_forecast      =KalmanFilterOutput.S_forecast; 
P_forecast      =KalmanFilterOutput.P_forecast;

K               =KalmanFilterOutput.K; 
J               =KalmanFilterOutput.J;

%
[nS,T]          =size(S_filter); 

S_smooth        =nan(nS,T); 
P_smooth        =cell(T,1); 
PP_smooth       =cell(T,1);

%t=T
S_smooth(:,T)   =S_filter(:,T); 
P_smooth{T}     =P_filter{T};
PP_smooth{T}    =(eye(nS)-K{end}*C)*A*P_filter{T}; %Shumway(1988)

%t<T
for t=T-1:-1:1
    
    S_smooth(:,t) =S_filter(:,t)+J{t}*(S_smooth(:,t+1)-S_forecast(:,t));
    P_smooth{t}   =P_filter{t}-J{t}*(P_forecast{t}-P_smooth{t+1})*J{t}';
    if t>1 
        
        PP_smooth{t} =P_filter{t}*J{t-1}'+...
                      J{t}*(PP_smooth{t+1}-A*P_filter{t})*J{t-1}'; %Shumway(1988)
    end
end

KalmanSmootherOutput.S_smooth       =S_smooth; 
KalmanSmootherOutput.P_smooth       =P_smooth;
KalmanSmootherOutput.PP_smooth      =PP_smooth;

% -----*-----*-----*-----*-----*-----*-----*-----*-----*-----*-----*----- %