

%DYNAMIC FACTOR MODEL WITH BLOCK RESTRICTIONS

% built on Banbura & Modugno (2010)
% miranda 2014 smirandaagrippino@london.edu

clear
clc

addpath([pwd '/subroutines/'])

countryQ='US'; 


DFMinit.StartEst =[1991 1]; 
load TEST_USdata %already transformed in growth rates


%-------------------------------------------------------------------------
% set parameters

DFMinit.nM          =nM; %number of monthly variables
DFMinit.max_iter    =500; %EM iterations
%

DFMinit.blocks      =blocks;


%--------------------------------------------------------------------------

                
DFMinit.r           =[1 1 1]; %number of factors per block
DFMinit.p           =2; %number of lags in factors VAR
DFMinit.s           =1; %lags in idio dynamic


%MODEL ESTIMATES (in-sample)
% X         =matrix of data at mixed (monthly and quarterly) frequency
% DFMinit   =structure with model parameters
% []        =leave empty to estimate new state space parameters


DFMend = EM_DHDFM_SS_block(X,DFMinit,[]);


Factor = DFMend.States(:,1);
