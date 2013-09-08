function NNIA()
% NNIA.m
% NNIA performs experiment presented in our paper:
% Maoguo Gong, Licheng Jiao,Haifeng Du, Liefeng Bo. Multiobjective Immune Algorithm with Nondominated Neighbor-based Selection.
% Evolutionary Computation Journal, MIT Press, 2007, in press.
% avialable from http://see.xidian.edu.cn/iiip/mggong/publication.htm
%
% Authors: Maoguo Gong and Licheng Jiao
% April 7, 2006
% Copyright (C) 2005-2006 by Maoguo Gong (e-mail: gong@ieee.org)
% revised by fuellee
%
% REFERENCE
% Please refer to the following paper if you use the toolbox NNIA
% Maoguo Gong, Licheng Jiao,Haifeng Du, Liefeng Bo. Multiobjective Immune Algorithm with Nondominated Neighbor-based Selection.
% Evolutionary Computation Journal, MIT Press, 2007, in press.
% NNIA is Copyrighted by the Authors.
%
% RELEASE
% version 1.0, 2006/04/10, tested on Matlab 7.0
% license granted for research use ONLY
% Copyright (C) 2005-2007 by Maoguo Gong (e-mail: gong@ieee.org)
% Please refer to the Readme.txt for a detailed description.
%--------------------------------------------------------------------------
clear all;
%--------------------------------------------------------------------------
fprintf('Nondominated Neighbor Immune Algorithm (NNIA) min\n');
fprintf('Authour: Maoguo Gong and Licheng Jiao\n');
fprintf('Last Modified: Oct. 10, 2007\n');
print_EMOinstruction;
%% display the instruction for running the programming.
%--------------------------------------------------------------------------
% TestNO=input('press the enter key after inputting the serial number of test problem:');
% Trial=input('input the number of independent runs:');
TestNO=18;
Trial=1;
Gmax=500;                              % maximum number of iterations(generations) default:500
n_D=100;                               % (maximum) size of dominant population
n_A=20;                                 % size of active population
CS=100;                                % clonal scale
[bu,bd,testfunction]=getbud(TestNO);   % bu denotes the upper boundary of variable;bd denotes the nether boundary of variable;bu and bd has the same dimensionality.
c=size(bu,2);
pm=1/c;
if bu==bd 
    return;end
%--------------------------------------------------------------------------
NS = zeros(1,Trial);
% paretof=[paretof;MEpa];
runtime = zeros(1,Trial);
Clonetime = zeros(1,Trial);
DAStime = zeros(1,Trial);
PNmtime = zeros(1,Trial);
DONjudtime = zeros(1,Trial);
RNDCDtime = zeros(1,Trial);
paretof=[];

for trial=1:Trial
    timerbegin=clock;
    %--------------------------------------------------------------------------
    % step 1: init population (n_D antibodies)
    POP=bsxfun(@times,rand(n_D,c),(bu-bd)) + ones(n_D,1)*bd;
    %--------------------------------------------------------------------------
    pa=OVcom(POP,TestNO); % pa: the current trade-off pareto front 
    DON = Dominant_Antibodies_index(pa);
    MEpa=pa(DON,:);MEPOP=POP(DON,:);
    
    [ClonePOP,Clonepa,~]=Update_Dominant_Population(MEPOP,MEpa,n_A);

    it=0;Cloneti=0;DASti=0;PNmti=0;DONjudti=0;RNDCDti=0;
    while it<Gmax   
        %--------------------------------------------------------------------------
        % cloneover=[];
        [cloneover,Clonet]=Clonef(ClonePOP,Clonepa,CS);
        [cloneover,DASt]=Recombinationf(cloneover,bu,bd,ClonePOP);
        [cloneover,PNmt]=Mutationf(cloneover,bu,bd,pm);
        %--------------------------------------------------------------------------
        clonepa=OVcom(cloneover,TestNO);
        NPOP=[MEPOP;cloneover];
        Npa=[MEpa;clonepa];
        [NDON,DONjudt]=Dominant_Antibodies_index(Npa);
        Nnodom=(NDON==1);
        NEpa=Npa(Nnodom,:);NEPOP=NPOP(Nnodom,:);
        % Nnumnod=size(NEPOP,1);
        [MEPOP, MEpa, ~]=Update_Dominant_Population(NEPOP,NEpa,n_D);
        numnod=size(MEPOP,1);
        [ClonePOP,Clonepa,RNDCDt]=Update_Dominant_Population(MEPOP,MEpa,n_A);
        %Update Dominant Population
        %--------------------------------------------------------------------------
        it=it+1;
        Cloneti=Cloneti+Clonet;DASti=DASti+DASt;PNmti=PNmti+PNmt;
        DONjudti=DONjudti+DONjudt;RNDCDti=RNDCDti+RNDCDt;

        fprintf('time: %d   generation: %d    number of nodominate:  %d\n',trial,it,numnod);

    end  %the end of iterations
    %--------------------------------------------------------------------------
    %Save the output solutions
    [NS(trial),NF]=size(MEpa); %used in eval! can't be removed
    paretof=[paretof;MEpa];
    runtime(trial)=etime(clock,timerbegin);
    Clonetime(trial)=Cloneti;DAStime(trial)=DASti;PNmtime(trial)=PNmti;
    DONjudtime(trial)=DONjudti;RNDCDtime(trial)=RNDCDti;

    Datime=date;
    Datime(size(Datime,2)-4:size(Datime,2))=[];%%test date e.g: 08-sep
    TestTime=clock;%%test time
    TestTime=[num2str(TestTime(4)),'-',num2str(TestTime(5))];
    Method='NNIA';

    eval(['save ', testfunction Method Datime TestTime ,' Method Gmax n_A CS paretof runtime trial NS NF TestTime Datime testfunction ']) ;
end  %the end of runs
%--------------------------------------------------------------------------
Frontshow(MEpa);% plot the Pareto fronts solved by the last run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DA,time]=Dominant_Antibodies_index(pa) 
% DA: Dominant Antibodies index;
% time: excution time
tic;
[N,C]=size(pa);
DA=true(N,1);
for i=1:N
    % i_antibody = pa(i,:);% i th antibody, to be check wether it's domainant
    temppa=pa;
    temppa(i,:)=[]; % delete i th antibody
    for j=1:C
        temppa=temppa(temppa(:,j)<=pa(i,j), :);
    end; % temppa: antibodies dominanted by or equal to the i th antibody.
    % temppa = temppa(sum(bsxfun(@le, temppa, i_antibody), 2)==C,:); 
    if ~isempty(temppa)
        for k=1:C
            Lessthan=(temppa(:,k)<pa(i,k));%lessthen: antibodies dominanted by the i th antibody.
            if ~isempty(Lessthan)
                DA(i)=false;break;end
        end
        % DA(i) = isempty(find(bsxfun(@lt, temppa, i_antibody),1)); 
    end    
end
time=toc;
%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EPOP,Epa,time]=Update_Dominant_Population(POP,pa,n_D) 
%--------------------------------------------------------------------------
tic;
% remove duplicated antibodies
[POP,index] = unique(POP,'rows');
pa = pa(index,:);
[Ns,~]=size(pa);

if Ns>n_D
    [POP,pa,padis]=CDAf(POP,pa); % pasid : crowding-distance
    [~,ss]=sort(-padis);
    EPOP=POP(ss(1:n_D),:);Epa=pa(ss(1:n_D),:);
else
    EPOP=POP;Epa=pa;
end
time=toc;
%%-------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [NPOP,time]=Clonef(POP,pa,CS) 
%-------------p-------------------------------------------------------------
tic
N = size(POP);
[POP,~,padis]=CDAf(POP,pa); % padis : crowding distance
aa=(padis==inf);
bb=(padis~=inf);
if ~isempty(bb)
    padis(aa)=2*max(padis(bb)); % set inf crowding-distance antibody to 2 times of max remain antibody
    NC=ceil(CS*padis./sum(padis));
else
    NC=ceil(CS/length(aa))+zeros(1,N);
end
NPOP=[];
for i=1:N
    NiPOP=ones(NC(i),1)*POP(i,:);
    NPOP=[NPOP;NiPOP];
end
time=toc;
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [POP,pa,padis]=CDAf(tPOP,tpa) % pasid : crowding-distance
%--------------------------------------------------------------------------
[Ns,C]=size(tpa);
for i=1:C
    [~,l]=sort(tpa(:,i));
    N=tpa(l,:);M=tPOP(l,:);
    tpa=N;tPOP=M;
    tpa(1,C+1)=Inf;tpa(Ns,C+1)=Inf;
    tpai1=tpa(3:Ns,i);tpad1=tpa(1:(Ns-2),i);
    fimin=min(tpa(:,i));fimax=max(tpa(:,i));
    tpa(2:(Ns-1),C+1)=tpa(2:(Ns-1),C+1)+(tpai1-tpad1)/(0.0001+fimax-fimin);
end
pa=tpa(:,1:C);POP=tPOP;padis=tpa(:,C+1);  
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Recombinationf(POP,bu,bd,CPOP)
%--------------------------------------------------------------------------
tic;
eta_c=15;
[N,C]=size(POP);
NPOP=POP;
for i=1:N
    r1=rand;
    if r1<=1.1
        aa=randperm(size(CPOP,1));bb=aa(1);
        for j=1:C
            par1=POP(i,j);par2=CPOP(bb,j);
            yd=bd(j);yu=bu(j);
            r2=rand;
            if r2<=0.5
                if abs(par1-par2)>10^(-14)
                    y1=min(par1,par2);y2=max(par1,par2);
                    if (y1-yd)>(yu-y2)
                        beta=1+2*(yu-y2)/(y2-y1);
                    else
                        beta=1+2*(y1-yd)/(y2-y1);
                    end
                    expp=eta_c+1;beta=1/beta;alpha=2.0-beta^(expp);
                    r3=rand;
                    if r3<=1/alpha
                        alpha=alpha*r3;expp=1/(eta_c+1.0);
                        betaq=alpha^(expp);
                    else
                        alpha=1/(2.0-alpha*r3);expp=1/(eta_c+1);
                        betaq=alpha^(expp);
                    end
                    chld1=0.5*((y1+y2)-betaq*(y2-y1));
                    chld2=0.5*((y1+y2)+betaq*(y2-y1));   
                    if rand<=0.5
                        aa=max(chld1,yd);NPOP(i,j)=min(aa,yu);
                    else
                        aa=max(chld2,yd);NPOP(i,j)=min(aa,yu);
                    end
                end  
            end
        end
    end
end
time=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [NPOP,time]=Mutationf(POP,bu,bd,pm)
%--------------------------------------------------------------------------
tic;
[N,C]=size(POP);
eta_m=20;
NPOP=POP;
for i=1:N
    for j=1:C
        r1=rand;
        if r1<=pm
            y=POP(i,j);
            yd=bd(j);yu=bu(j);
            if y>yd
                if (y-yd)<(yu-y)
                    delta=(y-yd)/(yu-yd);
                else
                    delta=(yu-y)/(yu-yd);
                end
                r2=rand;
                indi=1/(eta_m+1);
                if r2<=0.5
                    xy=1-delta;
                    val=2*r2+(1-2*r2)*(xy^(eta_m+1));
                    deltaq=val^indi-1;
                else
                    xy=1-delta;
                    val=2*(1-r2)+2*(r2-0.5)*(xy^(eta_m+1));
                    deltaq=1-val^indi;
                end
                y=y+deltaq*(yu-yd);
                NPOP(i,j)=min(y,yu);NPOP(i,j)=max(y,yd);
            else
                NPOP(i,j)=rand*(yu-yd)+yd;
            end
        end
    end
end
time=toc;
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Frontshow(Epa)
%--------------------------------------------------------------------------
if size(Epa,2)==2 %%plot pareto fronts in 2-D
    f1=Epa(:,1);f2=Epa(:,2);
    plot(f1,f2,'r*'); grid on;
    xlabel('Function 1');
    ylabel('Function 2');
    hold on   
elseif size(Epa,2)>=3 %%plot pareto fronts in 3-D
    f1=Epa(:,1);f2=Epa(:,2);f3=Epa(:,3);
    plot3(f1,f2,f3,'kd'); grid on;
    xlabel('Function 1');
    ylabel('Function 2');
    zlabel('Function 3');
    hold on   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
