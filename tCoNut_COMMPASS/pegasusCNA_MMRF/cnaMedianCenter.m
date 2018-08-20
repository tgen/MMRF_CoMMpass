function [cna,amp,del,dlr,log2fc,bafFreq]=cnaMedianCenter(normal,tumor,smWin,fcThresh,assayID,res,readDepth,intronLength,hetFile,hetDepthN,hetDepthT,hetDev,targetsFile)
% [cna,amp,del,noise,log2fc]=cnaMedianCenter(normal,tumor,smWin,fcThresh,assayID,res,readDepth,intronLength,hetFile,hetDepth,hetDev,targetsFile)
%
%   Copy Number Analysis (CNA) for next-generation sequencing (NGS) data.
%   Algorithm is designed to work on exome and long-insert, low-coverage
%   whole genome sequencing (WGS).
%
%   INPUTS:
%       normalDat is a N x 3 matrix where first column is chromosome, second
%       column is position and third column is number of reads at that
%       position.
%
%       tumorDat is a N x 3 matrix where first column is chromosome, second
%       column is position and third column is number of reads at that
%       position.
%
%       smWin is the size of a sliding window (e.g. 6 for exomes and 19 for WGS).
%       fcThresh is the threshold where 'amp' and 'del' are determined for
%       plotting purposes.
%
%   OUTPUT:
%      cna is a N x 3 matrix where the first column is chromosome, second
%      column is position and third column is log2 fold-change. 
%
%      amp is a N x 3 matrix where the first column is chromosome, second
%      column is position and third column is log2 fold-change.  amp
%      contains only those locations where log2 fold-change > fcThresh.
%
%      del is a N x 3 matrix where the first column is chromosome, second
%      column is position and third column is log2 fold-change.  amp
%      contains only those locations where log2 fold-change < fcThresh.
%
%   David Craig
%   Modified by Jessica Aldrich
%   TGen
%   January 16, 2014
%
%   Modified by Jessica Aldrich
%   April 2015
%   Added BAF calculation and output

% Calculate fold change
paira=cat(2,normal(:,1:3),tumor(:,3));


%dat = paira;

%load exomeTargets
switch lower(assayID)
    case 'exome'
        overlap=load(targetsFile);
        paira(~logical(overlap(:,4)),:)=[];
        
    case 'genome'
        
        if nargin==13
            pos=load(targetsFile);
            paira=paira(logical(pos),:);
        
        end
        low=10;
        inc_ploid=0.8;
        %remove aneuploidy regions in normal sample in autosomal chrs
        autochr = paira(:,1)<23;
        autochrs = paira(autochr,:);
        norm_mode=mode(autochrs(autochrs(:,3)>low,3));
        a=autochrs(autochrs(:,3)>norm_mode*inc_ploid & autochrs(:,3) < norm_mode*(1/inc_ploid),:);
        
        %remove aneuploidy regions in normal sample in sex chrs
        sexchr = paira(:,1)>=23;
        sexchrs = paira(sexchr,:);
        norm_mode=mode(sexchrs(sexchrs(:,3)>low,3));
        b=sexchrs(sexchrs(:,3)>norm_mode*inc_ploid & sexchrs(:,3) < norm_mode*(1/inc_ploid),:);
        
        paira = [a;b];
end

%paira=paira(paira(:,3)>0.7*hetDepthN,:);

pairb=zeros(size(paira));
count=1;

%per chromosome
for chr=1:24
    a=1;
    b=1;
    chrom=paira(paira(:,1)==chr,:);
   % t=tic;
    while b<=size(chrom,1) && a<size(chrom,1)

        if a+res<=size(chrom,1)
            b=a+res;

            while b<size(chrom,1) && chrom(b+1,2)-chrom(b,2) < intronLength && sum(chrom(a:b,3)) < readDepth
                b=b+1;
                if b>size(chrom,1)
                    b=size(chrom,1);
                end
             % toc(t)
            end

            npt = [chrom(a,1) round(median(chrom(a:b,2))) round(mean(chrom(a:b,3:4)))];  %%%%MEAN CAN BE CHANGED - SUM, MEAN, MEDIAN
            pairb(count,:) = npt;
            count=count+1;

        else
            b=size(chrom,1);
            npt = [chrom(a,1) round(median(chrom(a:b,2))) round(mean(chrom(a:b,3:4)))];  %%%%MEAN CAN BE CHANGED - SUM, MEAN, MEDIAN
            pairb(count,:) = npt;
            count=count+1;
        end

        a=b;

    end
    %tt=toc(t)

end

pairb(~logical(pairb(:,1)),:)=[];

for j=3:4
   pairb(isnan(pairb(:,j)),j)=0;
   pairb(:,j)=pairb(:,j)./sum(pairb(:,j));
end

% pairb(isnan(pairb(:,3)),3)=0;
% pairb(:,3)=pairb(:,3)./sum(pairb(:,3));
% 
% pairb(isnan(pairb(:,4)),4)=0;
% pairb(:,4)=pairb(:,4)./(sum(pairb(:,4))*.7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



paira=pairb;
paira(:,4)=log2(paira(:,4)./paira(:,3));
% paira(:,4)=roll_med(paira(:,4),smWin);
% 
cna=paira(isfinite(paira(:,4)),[1,2,4]);

%cna = paira(:,[1 2 4]);
%%%%%%%MEDIAN CENTERING%%%%%%%%%

log2fc=[];

if hetFile   %if hetFile is 0 - skip median centering
    switch lower(assayID)
        case 'exome'
            fid = fopen(hetFile);
            vcf = textscan(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f','headerlines',1);
            fclose(fid);
            vcfMat =cell2mat(vcf);
            
            %adjust to expand or contract window of excepted allele frequency
            minAllFreq=0.5-hetDev;
            maxAllFreq=0.5+hetDev;
            
            allFreq = [vcfMat(:,4)./sum(vcfMat(:,4:5),2) vcfMat(:,7)./sum(vcfMat(:,7:8),2)];
            bafFreq = [vcfMat(:,1) round(vcfMat(:,2)/100)*100 vcfMat(:,8)./sum(vcfMat(:,7:8),2)];
            sumDepth = [sum(vcfMat(:,4:5),2) sum(vcfMat(:,7:8),2)];
            ind = sumDepth(:,1)>hetDepthN & sumDepth(:,1)< hetDepthN+(0.5*hetDepthN) & ...
                  sumDepth(:,2)>hetDepthT & sumDepth(:,2)< hetDepthT+(0.5*hetDepthT) & ...
                  allFreq(:,2)>minAllFreq & allFreq(:,2)<maxAllFreq;
%             ind = sumDepth(:,1)>hetDepth & sumDepth(:,2)>hetDepth & allFreq(:,1)>minAllFreq & ...
%                 allFreq(:,1)<maxAllFreq & allFreq(:,2)>minAllFreq & allFreq(:,2)<maxAllFreq;
        case 'genome'
            fid = fopen(hetFile);
            vcf = textscan(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f','headerlines',1);
            fclose(fid);
            vcfMat =cell2mat(vcf);
            
            minAllFreq=0.5-hetDev;
            maxAllFreq=0.5+hetDev;
            
            allFreq = [vcfMat(:,4)./sum(vcfMat(:,4:5),2) vcfMat(:,7)./sum(vcfMat(:,7:8),2)];
            bafFreq = [vcfMat(:,1) round(vcfMat(:,2)/100)*100 vcfMat(:,8)./sum(vcfMat(:,7:8),2)];
            sumDepth = [sum(vcfMat(:,4:5),2) sum(vcfMat(:,7:8),2)];
            %hetDepthN = median(sumDepth(:,1));
            %hetDepthT = median(sumDepth(:,2));
           
	    ind = sumDepth(:,1)>hetDepthN & sumDepth(:,1)< hetDepthN+(0.5*hetDepthN) & ...
                  sumDepth(:,2)>hetDepthT & sumDepth(:,2)< hetDepthT+(0.5*hetDepthT) & ...
                  allFreq(:,2)>minAllFreq & allFreq(:,2)<maxAllFreq; 

	   %ind = sumDepth(:,1)>hetDepthN-(0.3*hetDepthN) & sumDepth(:,1)< hetDepthN+(0.3*hetDepthN) & ...
           %       sumDepth(:,2)>hetDepthT-(0.3*hetDepthT) & sumDepth(:,2)< hetDepthT+(0.3*hetDepthT) & ...
           %       allFreq(:,2)>minAllFreq & allFreq(:,2)<maxAllFreq;
            
%             sumDepth = [sum(vcfMat(:,4:5),2) sum(vcfMat(:,7:8),2)];
%             ind = sumDepth(:,2)>hetDepthT & sumDepth(:,2)<3*hetDepthT & ...
%                 vcfMat(:,9) < hetDev;
                                           
    end

    vcfMatFilt = vcfMat(ind,:);
    vcfMatFilt = [vcfMatFilt(:,1) round(vcfMatFilt(:,2)/10000)*10000 vcfMatFilt(:,3:end)];

    vcfMatFiltVect = vcfMatFilt(:,1)*1000000000+vcfMatFilt(:,2);
    nm = cna(:,1)*1000000000+round(cna(:,2)/10000)*10000;

    [~,hets] = intersect(nm,vcfMatFiltVect);
    log2fc = cna(hets,:); 
   
    if length(hets) > 1
        cna(:,3) = cna(:,3)-mean(log2fc(:,3));
    end
end

% Median smoothing per chromosome
for i = 1:length(unique(cna(:,1)))
    chr = cna(:,1)==i;
    if sum(chr) > 2*smWin
        cna(chr,3) = roll_med(cna(chr,3),smWin);
    end
end


amp = cna(cna(:,3)>fcThresh,:);
del = cna(cna(:,3)<-fcThresh,:);


%%%%%MEASURE OF NOISE%%%%%%%%%
cna = cna(isfinite(cna(:,3)),:);
dlr = dlrs(cna(1:(2*smWin)+1:end,3));

%vcfMatFilt = vcfMat(ind,:);
% vcfMat = [vcfMat(:,1) round(vcfMat(:,2)/100)*100 vcfMat(:,3:end)];
% vcfMatVect = vcfMat(:,1)*1000000000+vcfMat(:,2);
% nm = normal(:,1)*1000000000+normal(:,2);
% [~,hets] = intersect(nm,vcfMatVect);
% norm_hets = normal(hets,3)./sum(normal(:,3));
% tumor_hets = tumor(hets,3)./sum(tumor(:,3));
% log2noise = log2(tumor_hets./norm_hets);
% log2noise = [normal(hets,1:2) log2noise];
% log2noise = log2noise(~isnan(log2noise(:,3)),:);
% log2noise = log2noise(~isinf(log2noise(:,3)),:);
% 
% chrvec = unique(log2noise(:,1));
% for i = 1:length(chrvec)
%     
%     ind = log2noise(:,1)==chrvec(i);
%     log2noise(ind,3) = log2noise(ind,3)-mean(log2noise(ind,3));
%     
% end
%noise=std(log2noise(:,3));
%
% chr1pos1 = dat(dat(:,1)==1&dat(:,2)<25000000,:);
% chr1pos2 = dat(dat(:,1)==1&dat(:,2)>225000000,:);
% log2pos1 = log2(chr1pos1(:,4)./chr1pos1(:,3));
% log2pos2 = log2(chr1pos2(:,4)./chr1pos2(:,3));
% 
% log2pos1 = log2pos1(isfinite(log2pos1));
% log2pos2 = log2pos2(isfinite(log2pos2));
% 
% log2pos1 = log2pos1-mean(log2pos1);
% log2pos2 = log2pos2-mean(log2pos2);
% 
% noise = std([log2pos1;log2pos2]);


end



