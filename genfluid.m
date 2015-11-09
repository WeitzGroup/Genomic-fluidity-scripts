function [gfluid,gfl_var]=genfluid(shared,total)
% function [gfluid,gfl_var]=genfluid(shared,total);
% computes fuidity estimate and its variance estimate
% input: shared: NxN matrix, with N number of sampled genomes
%         - element (i,j) contains number of shared genes (gene families)
%                                  between genome i and j
%        total: NxN matrix, with N number of sampled genomes
%         - element (i,j) contains number of total genes (gene families)
%                                  between genome i and j
% output: gfluid: estimate for fluidity of genomes group
%         gfl_var: estimate for variance of fluidity estimate
%
% Contact: Bart Haegeman (bart.haegeman@inria.fr)
% Contact: Joshua Weitz (jsweitz@gatech.edu)
%
% Fluidity estimator for Kislyuk et al. [1]
% [1] Genomic fluidity: an integrative view of gene diversity within microbial populations. 
% Andrey O Kislyuk, Bart Haegeman, Nicholas H Bergman and Joshua S Weitz. 
% BMC Genomics 2011, 12:32. doi:10.1186/1471-2164-12-32 
% http://www.biomedcentral.com/1471-2164/12/32

% compute number of genomes N
ss=size(shared);
ss2=size(total);
if ss(1)~=ss(2) || ss2(1)~=ss2(2) || ss(1)~=ss2(1),
    gfluid=NaN;
    gfl_var=NaN;
    disp('matrices shared and total are not compatible')
    return
end

% compute similarity matrix simi
N=ss(1);
simi=zeros(N,N);
for cnt1=1:N-1,
    for cnt2=cnt1+1:N,
        simi(cnt1,cnt2)=2*shared(cnt1,cnt2) ...
            /(total(cnt1,cnt2)+shared(cnt1,cnt2));
    end
end
gfluid=1-sum(sum(simi))*2/N/(N-1);

% compute one-leave-out statistics olos
olos=zeros(1,N);
simia=zeros(N-1,N-1);
for cnt=1:N,
    simia=simi([1:cnt-1 cnt+1:N],[1:cnt-1 cnt+1:N]);
    olos(cnt)=sum(sum(simia))*2/(N-1)/(N-2);
end
gfl_var=(N-1)*var(olos,1);
