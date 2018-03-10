function [pdfxy,x] = MakeDensity2d(Name,T,n)
%
%  Make artificial normal mixture density and
%  generate associated normal mixture distributed random numbers for
%  the twelve density examples used in Wand and Jones, JASA, 1993.
%  
% 	Input
%		Name   String: 'UncorNormal','CorNormal','Skewed','Kurtotic',
%                      'BimodalI','BimodalII','BimodalIII','BimodalIV',
%                      'SkewedUnimodal','StronglySkewed','Kurtotic',
%                      'TrimodalI','TrimodalII','TrimodalIII','Quadrimodal'
%		n	   Desired sampling length.
%		T	   Output density length.
%
%	Output
%		pdfxy  Theoretical density.
%       x      Random samples.
%
%  Copyright (c) 2018 Fabien Navarro

if nargin > 1,
  t = linspace(-3,3,T);
end
[X1,X2] = meshgrid(t,t);
dom = [X1(:),X2(:)];
  if strcmp(Name,'UncorNormal'),
    mu = [0 0];
    Sigma = [0.25 0; 0 1];
    pdF = mvnpdf(dom,mu,Sigma);
    pdfxy = reshape(pdF,T,T);
    x = mvnrnd(mu,Sigma,n);
  elseif strcmp(Name,'CorNormal'),
    mu = [0 0];
    Sigma = [1 7/10; 7/10 1];
    pdF = mvnpdf(dom,mu,Sigma);
    pdfxy = reshape(pdF,T,T);
    x = mvnrnd(mu,Sigma,n);
  elseif strcmp(Name,'Skewed'),
    mu1 = [0 0];
    mu2 = [1/2 1/2];
    mu3 = [13/12 13/12];
    Sigma(:,:,1) = [1 0;0 1];
    Sigma(:,:,2) = [4/9 0; 0 4/9];
    Sigma(:,:,3) = [25/81 0; 0 25/81];
    mu(:,1)=mu1;mu(:,2)=mu2;mu(:,3)=mu3;
    p = [1 1 3]/5;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2)); 
    pdF3 = mvnpdf(dom,mu3,Sigma(:,:,3));
    pdF  = p(1)*pdF1+p(2)*pdF2+p(3)*pdF3;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);
  elseif strcmp(Name,'Kurtotic'),
    mu1 = [0 0];
    mu2 = [0 0];
    Sigma(:,:,1) = [1 1;1 4];
    Sigma(:,:,2) = [4/9 -1/9; -1/9 1/9];
    mu(:,1)=mu1;mu(:,2)=mu2;
    p = [2 1]/3;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2));        
    pdF  = p(1)*pdF1+p(2)*pdF2;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);
  elseif strcmp(Name,'BimodalI'),
    mu1 = [-1 0];
    mu2 = [1 0];
    Sigma(:,:,1) = [4/9 0;0 4/9];
    Sigma(:,:,2) = [4/9 0;0 4/9];
    mu(:,1)=mu1;mu(:,2)=mu2;
    p = [1 1]/2;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2));  
    pdF  = p(1)*pdF1+p(2)*pdF2;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);
  elseif strcmp(Name,'BimodalII'),
    mu1 = [-3 0]/2;
    mu2 = [3 0]/2;
    Sigma(:,:,1) = [1/16 0;0 1];
    Sigma(:,:,2) = [1/16 0;0 1];
    mu(:,1)=mu1;mu(:,2)=mu2;
    p = [1 1]/2;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2));  
    pdF  = p(1)*pdF1+p(2)*pdF2;
    pdfxy = reshape(pdF,length(t),length(t));
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);  
  elseif strcmp(Name,'BimodalIII'),
    mu1 = [-1 1];
    mu2 = [1 -1];
    Sigma(:,:,1) = [4/9 4/15;4/15 4/9];
    Sigma(:,:,2) = [4/9 4/15;4/15 4/9];
    mu(:,1)=mu1;mu(:,2)=mu2;
    p = [1 1]/2;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2));  
    pdF  = p(1)*pdF1+p(2)*pdF2;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);
  elseif strcmp(Name,'BimodalIV'),
    mu1 = [1 -1];
    mu2 = [-1 1];
    Sigma(:,:,1) = [4/9 14/45;14/45 4/9];
    Sigma(:,:,2) = [4/9 0;0 4/9];
    mu(:,1)=mu1;mu(:,2)=mu2;
    p = [1 1]/2;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2));  
    pdF  = p(1)*pdF1+p(2)*pdF2;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);
  elseif strcmp(Name,'TrimodalI'),
    mu1 = [-6/5 6/5];
    mu2 = [6/5 -6/5];
    mu3 = [0 0];
    Sigma(:,:,1) = [9/25 27/250; 27/250 9/25];
    Sigma(:,:,2) = [9/25 -27/125; -27/125 9/25];
    Sigma(:,:,3) = [1/16 1/80; 1/80 1/16];
    mu(:,1)=mu1;mu(:,2)=mu2;mu(:,3)=mu3;
    p = [9 9 2]/20;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2)); 
    pdF3 = mvnpdf(dom,mu3,Sigma(:,:,3));
    pdF  = p(1)*pdF1+p(2)*pdF2+p(3)*pdF3;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);  
 elseif strcmp(Name,'TrimodalII'),
    mu1 = [-6/5 0];
    mu2 = [6/5 0];
    mu3 = [0 0];
    Sigma(:,:,1) = [9/25 63/250; 63/250 9/25];
    Sigma(:,:,2) = [9/25 63/250; 63/250 9/25];
    Sigma(:,:,3) = [9/25 -63/250; -63/250 9/25];
    mu(:,1)=mu1;mu(:,2)=mu2;mu(:,3)=mu3;
    p = [1 1 1]/3;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2)); 
    pdF3 = mvnpdf(dom,mu3,Sigma(:,:,3));
    pdF  = p(1)*pdF1+p(2)*pdF2+p(3)*pdF3;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);    
 elseif strcmp(Name,'TrimodalIII'),
    mu1 = [-1 0];
    mu2 = [1 2*sqrt(3)/3];
    mu3 = [1 -2*sqrt(3)/3];
    Sigma(:,:,1) = [9/25 63/250; 63/250 49/100];
    Sigma(:,:,2) = [9/25 0; 0 49/100];
    Sigma(:,:,3) = [9/25 0; 0 49/100];
    mu(:,1)=mu1;mu(:,2)=mu2;mu(:,3)=mu3;
    p = [3 3 1]/7;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2)); 
    pdF3 = mvnpdf(dom,mu3,Sigma(:,:,3));
    pdF  = p(1)*pdF1+p(2)*pdF2+p(3)*pdF3;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);       
 elseif strcmp(Name,'Quadrimodal'),
    mu1 = [-1 1];
    mu2 = [-1 -1];
    mu3 = [1 -1];
    mu4 = [1 1];
    Sigma(:,:,1) = [4/9 8/45; 8/45 4/9];
    Sigma(:,:,2) = [4/9 4/15; 4/15 4/9];
    Sigma(:,:,3) = [4/9 -28/90; -28/90 4/9];
    Sigma(:,:,4) = [4/9 -2/9; -2/9 4/9];
    mu(:,1)=mu1;mu(:,2)=mu2;mu(:,3)=mu3;mu(:,4)=mu4;
    p = [1 3 1 3]/8;
    pdF1 = mvnpdf(dom,mu1,Sigma(:,:,1)); 
    pdF2 = mvnpdf(dom,mu2,Sigma(:,:,2)); 
    pdF3 = mvnpdf(dom,mu3,Sigma(:,:,3));
    pdF4 = mvnpdf(dom,mu4,Sigma(:,:,4));    
    pdF  = p(1)*pdF1+p(2)*pdF2+p(3)*pdF3+p(4)*pdF4;
    pdfxy = reshape(pdF,T,T);
    obj = gmdistribution(mu',Sigma,p);
    x = random(obj,n);      
 else 
    disp('Unkwnown Density type.');
    disp('Allowable Names are:')
	disp('UncorNormal'),
    disp('CorNormal'),
    disp('Skewed'),
    disp('Kurtotic'),
    disp('BimodalI'),   
    disp('BimodalII'), 
    disp('BimodalIII'),
    disp('BimodalIV'),
    disp('TrimodalI'),
    disp('TrimodalII'),
    disp('TrimodalIII'),
    disp('Quadrimodal'),
 end
