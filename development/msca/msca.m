function [Model]=msca(X,DesignMat,FacB,FacW)

%   Multilevel SCA

%   inputs:
%       X           ((I x Ki) x J)     =   data
%       DesignMat   ((I x Ki) x I)     =   Matrix containing 0's and 1's defining which samples belong to which individual
%       FacB                           =   number of Between-individual factors
%       FacW                           =   number of Within-individual factors


%   outputs:
%       Model structure
%       Model.offset        (1 x J)             = overall offset of the data

%       Model.between.data                      =   Data used for the between-individual level
%       Model.between.scores                    =   Between-individual scores
%       Model.between.loadings                  =   Between-individual loadings
%       Model.between.percexp.total             =   Total percentage of explained variance
%       Model.between.percexp.perpc             =   Percentage of explained variance for each PC
% 
%       Model.within.data                       =   Data used for the within individual level
%       Model.within.scores                     =   Within-individual scores
%       Model.within.loadings                   =   Within-individual loadings
%       Model.within.percexp.total              =   Total amount of explained variance 
%       Model.within.percexp.perpc              =   Percentage of explained variance of each PC
%       Model.within.ind.'nind'.data            =   Data used for the within individual level for each individual
%       Model.within.ind.'nind'.scores          =   Within-individual scores for each individual
%       Model.within.ind.'nind'.percexp.total   =   Total amount of explained variance for each individual
%       Model.within.ind.'nind'.percexp.perpc   =   Percentage of explained variance of each PC for each individual
%   
%       I/O: [Model]=msca(X,DesignMat,FacB,FacW)
%
%   ==========================================================================
%   Copyright 2005 Biosystems Data Analysis Group ; Universiteit van Amsterdam
%   ==========================================================================

I   =   size(DesignMat,2);
J   =   size(X,2);
K   =   size(X,1);
Ki  =   [];
for i=1:I
    Ki(i)   =   length(find(DesignMat(:,i)==1));
end
clear i


%   Part I: Calculate Offset
offset  =   mean(X);

%   Subtract offset from Data
Xoff    =   X-ones(K,1)*offset;

%   Part II: Between individual scores
%   Calculate mean of each individual
for i=1:I
    [smpn,dummy]=find(DesignMat(:,i));
    Xind(i,:)=mean(Xoff(smpn,:));
end

Xindsc=Xind;


W   =   diag(sqrt(Ki));
[scb,ldb,ssqb,resb,reslmb,tsqlmb,tsqb] = pca(W*Xindsc,0,[],FacB);
scb =   inv(W)*scb;

%   Part III: Within individual scores
%   Center data for each individual
Xwit=Xoff-DesignMat*Xind;

Xwitsc=Xwit;


%   Calculate within-individual model
[scw,ldw,ssqw,resw,reslmw,tsqlmw,tsqw] = pca(Xwitsc,0,[],FacW);

Model=[];

Model.dimensions.I              =   I;
Model.dimensions.J              =   J;
Model.dimensions.K              =   K;
Model.dimensions.Ki             =   Ki;
Model.dimensions.FacB           =   FacB;
Model.dimensions.FacW           =   FacW;

Model.offset                    =   offset;

Model.between.data              =   Xindsc;
Model.between.scores            =   scb;
Model.between.loadings          =   ldb;
Model.between.percexp.total     =   (1-(ssq(Xindsc-scb*ldb')/ssq(Xindsc)))*100;
Model.between.percexp.perpc     =   ssqb([1:FacB],3);

Model.within.data               =   Xwitsc;
Model.within.scores             =   scw;
Model.within.loadings           =   ldw;
Model.within.percexp.total      =   (1-(ssq(Xwitsc-scw*ldw')/ssq(Xwitsc)))*100;
Model.within.percexp.perpc      =   ssqw([1:FacW],3);
for nind=1:I
    [indind,dummy]                             =   find(DesignMat(:,nind)==1);
    clear dummy
    Model.within.ind{nind}.data=Xwit(indind,:);
    Model.within.ind{nind}.scores         =   scw(indind,:);
    Model.within.ind{nind}.percexp.total  =   (1-(ssq(Xwitsc(indind,:)-scw(indind,:)*ldw')/ssq(Xwitsc(indind,:))))*100;
    for pc=1:FacW
        Model.within.ind{nind}.percexp.perpc(pc)   =   (1-(ssq(Xwitsc(indind,:)-scw(indind,pc)*ldw(:,pc)')/ssq(Xwitsc(indind,:))))*100;
    end
    Model.within.ind{nind}.percexp.perpc  =   Model.within.ind{nind}.percexp.perpc';
end


function t = ssq(a)
%SSQ    SSQ(A) is the sum of squares of the elements of matrix A.
t = sum(sum(a.^2));