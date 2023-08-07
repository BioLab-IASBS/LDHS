% % This code is a modified version of HS_2019_multiCRITERIA5.m from Tou, et. al. 2020
   function [Candidate,canSize,NC,totaltime,flag] =NEW_HS(data,dim_epi,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX)








HMCR=0.98;
PAR=0.35;
n=size(data,2);
State=data(:,n);
Mc = 4; % ??????
Candidate1=ones(CandidateSize,dim_epi+ Mc); %% ?????
canSize1=0;
Candidate2=ones(CandidateSize,dim_epi+ Mc); %% ?????
canSize2=0;
Candidate3=ones(CandidateSize,dim_epi+ Mc); %% ?????
canSize3=0;

Candidate4=ones(CandidateSize,dim_epi+ Mc); %% ?????
canSize4=0;

%% ---------------------------------------------------------------

SNPs=n-1;  %% ?SNP??

EliteSize=fix(HMS/5);  %% HMS????
Elite1=[]; %% ???????
Efit1=[];

Elite2=[]; %% ???????
Efit2=[];

Elite3 = []; %% ???????
Efit3 = [];


Center=[];%zeros(1,dim_epi+1); %% ???????
% maxIterForLocalSearch=min(fix(max_iter/5),3000);  %% ????????????????????????


%% ???
flag = -1;
X=zeros(HMS,dim_epi);
snp=[];

k=dim_epi;
[hl2,uperbound]=size(data);
lowerbound=-1*uperbound;
vardef(1,1:k)=uperbound;
vardef(2,1:k)=lowerbound;


[dummy,NoVar] = size(vardef);

diff = zeros(1,NoVar);
N2=4*HMS;
X = zeros(N2,NoVar);


%Find difference for scaling
for i=1:NoVar
    diff(1,i) = (vardef(2,i) - vardef(1,i));
end





name_snip =["rs6269","rs4680","rs10046","rs3020314","rs2234693","rs1543404","rs3798577","rs2747652","rs207747","rs2175898","rs19340799","rs1709182","rs9478249","rs1514348","rs532010","rs566351","rs660149","rs11571171","rs500760","rs858518","272428","rs858524","rs85852"];
gene_snip=["COMT","COMT","CYP19A1","ESR1","ESR1","ESR1","ESR1","ESR1","ESR1","ESR1","ESR1","ESR1","ESR1","PGR","ESR1","PGR","PGR","PGR","PGR","SHBG","TBC1D9B","FXR2","STS"];
location=["upstream","upstream","intron","intron","upstream","downstream","downstream","downstream","noncoding","upstream","upstream","upstream","intron","intron","upstream","intron","intron","downstream","downstream","upstream","intron","intron","intron"];
kk=2;





my_snp=[];
new_location=[];
new_gene_snip=[];
my_data=[];
temp(1)=gene_snip(1);

for i=1 : size(gene_snip,2)

    a=strcmp(gene_snip(i),temp);
    if a==0
        temp(kk)=gene_snip(i);
        kk=kk+1;
    end
end

tnn=0;
XX=[1,2];
%
size(temp , 2);
% temp(2);
for i9=1:size(temp , 2)
      temp(i9);

    number=0;
    total_number=0;
    gg_int=[];
    ex_int=[];
    number_intron=0;
    number_exon=0;
    for j=1:size(gene_snip , 2)
        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"intron")==1 )

%             size(gg_int)
%             size(gg_int , 2)+1

            gg_int(size(gg_int , 2)+1)=j;
%             gg_int;
            number_intron=number_intron+1;

        end

        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"exon")==1 )


            ex_int(size(ex_int , 2)+1)=j;

            number_exon=number_exon+1;
        end




        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"noncoding")==1 )


            gg_int(size(gg_int , 2)+1)=j;
%             gg_int;
            number_intron=number_intron+1;
        end



        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"coding")==1 )

            ex_int(size(ex_int , 2)+1)=j;
            number_exon=number_exon+1;
        end

        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"upstream")==1 )
%
            ex_int(size(ex_int , 2)+1)=j;
            number_exon=number_exon+1;
        end


        if(strcmp(gene_snip(j),temp(i9))==1 && strcmp(location(j),"downstream")==1 )

            ex_int(size(ex_int , 2)+1)=j;
            number_exon=number_exon+1;
        end
        if(strcmp(gene_snip(j),temp(i9))==1)

            total_number=total_number+1;
        end



    end

    number_exon;
    number_intron;

%     gg_int
%     ex_int
     total_number;

    number=ceil((total_number *92)/23);

    number;

    total_number;
    temp_ex =((number_exon *100 )/ total_number );


    temp_int  =((number_intron *100 )/ total_number );
    temp_int;
    temp_ex;
    if(temp_ex >0 && temp_int>0)
        d=444+1;
        d;
        temp_int =((number_intron *100 )/ total_number )/2 ;
        temp_ex  =((number_exon *100 )/ total_number ) + temp_int;
    end
    temp_ex;
    temp_int;
    g_ex= (temp_ex/100) * number;
    g_int=(temp_int/100)  * number;
    g_ex;
    g_int;

%     g_int;
%     gg_int;
    if(g_int>0)
        for dd1=1:g_int
            a1=ceil(rand*size(gg_int , 2));
            if a1==0
                a1=a1+1;
            end
            b=size(gg_int , 2);
    %         b;
    %         a1;
    %         gg_int(1,a1);
            my_snp(1)=gg_int(1,a1);
            for dd2=2:dim_epi
                a2=ceil(rand*size(gg_int , 2));
                my_snp(2)=gg_int(1,a2);
                while ismember(my_snp(dd2),my_snp(1:dd2-1))

                    my_snp(dd2)=ceil(rand*22);
                end
            end

            temp1=my_snp;
            my_snp=sort(my_snp);
            while ismember(my_snp,XX,"rows")
                jj=ceil(rand*dim_epi);
                my_snp(jj)=ceil(rand*22);
                temp1=my_snp;
                my_snp=sort(my_snp);
            end
            l=size(XX)+1;
            for i=1:dim_epi
                XX(l(1) ,i)=my_snp(i);
            end

            if my_snp == CX
                flag = 1111;
            end

            my_snp=[];
        end
    end


%     XX;


    g_ex;
    ex_int;
    if (g_ex>0)


        for dd1=1:g_ex
            a1=ceil(rand*size(ex_int , 2));
    %         b=size(ex_int , 2);
    %         b;
    %         a1;
            ex_int;

            my_snp(1)=ex_int(1,a1);
            for dd2=2:dim_epi
                a2=ceil(rand*size(ex_int , 2));
                my_snp(2)=ex_int(1,a2);
                while ismember(my_snp(dd2),my_snp(1:dd2-1))
                    my_snp(dd2)=ceil(rand*22);



                end
            end

            temp1=my_snp;
            my_snp=sort(my_snp);
            while ismember(my_snp,XX,"rows")
                jj=ceil(rand*dim_epi);
                my_snp(jj)=ceil(rand*22);
                temp1=my_snp;
                my_snp=sort(my_snp);
            end
            l=size(XX)+1;
            for i=1:dim_epi
                XX(l(1) ,i)=my_snp(i);
                XX;
            end

            if my_snp == CX
                flag = 1111;
            end

    %         HM(gg+1,:)=X(gg+1,:);  %% HM????????
            my_snp=[];
        end
    end


    X=XX(2,:);
    HM=XX(2,:);
    [Fit(1,1),Fit(1,2),Fit(1,3),Fit(1,4)] = multi_criteriaEvaluationFuns3(data(:,X(1,:)),State);

    for i=3:size(XX)

        X(i-1,:)=XX(i,:);

        HM(i-1,:)=XX(i,:);

        [Fit(i-1,1),Fit(i-1,2),Fit(i-1,3),Fit(i-1,4)] = multi_criteriaEvaluationFuns3(data(:,X(i-1,:)),State);

    end
    X;
    size(X);
end
%
% [t1,t2]=size(X);
%
% if(t1<400)
%     for k=1:400-t1
%         snp(1)=ceil(rand*22);
%         for j=2:dim_epi
%           snp(j)=ceil(rand*22);
%
%           while ismember(snp(j),snp(1:j-1))
%              snp(j)=ceil(rand*22);
%           end
%         end
%         temp_d=snp;
%         snp=sort(snp);
%         while ismember(snp,X,'rows')
%             j=ceil(rand*dim_epi);
%             snp(j)=ceil(rand*22);
%             temp_d=snp;
%             snp=sort(snp);
%         end
%
%         [v1,v2]=size(X);
%
%         X(v1+1 , :)=snp;
%         HM(v1+1 , :)=snp;
%
%
%
%
%         if snp == CX
%             flag = 1111;
%         end
%         [Fit(t1+1,1),Fit(t1+1,2),Fit(t1+1,3),Fit(t1+1,4)] = multi_criteriaEvaluationFuns3(data(:,X(v1+1 , :)),State);
%
%
%         snp=[];
%     end
%
% end
% X
% size(X)
%









T1 = max(Fit(:,1)) - min(Fit(:,1));
T2 = max(Fit(:,2)) - min(Fit(:,2));
T3 = max(Fit(:,3)) - min(Fit(:,3));
T4 = max(Fit(:,4)) - min(Fit(:,4));




X2=X(HMS+1:2*HMS,:);

HM2=HM(HMS+1:2*HMS,:);
Fit2=Fit(HMS+1:2*HMS,:);

X3 = X(2*HMS+1:3*HMS,:);
HM3 = HM(2*HMS+1:3*HMS,:);
Fit3 = Fit(2*HMS+1:3*HMS,:);

X4 = X(3*HMS+1:4*HMS,:);
HM4 = HM(3*HMS+1:4*HMS,:);
Fit4 = Fit(3*HMS+1:4*HMS,:);

X = X(1:HMS,:);
Fit = Fit(1:HMS,:);
HM = HM(1:HMS,:);


NC=HMS;
 LT=0;
%%-------------------------------------------------------------------------
tic;
kk=1;


%%-------------------------------------------------------------------------
tic;
while NC <= max_iter
    R4 = ceil(rand*4);
    if R4 == 1
        [SF1, sind1] = sort(Fit(:,1));
          Xbest1 = HM(sind1(1),:);
    elseif R4 == 2
        [SF2, sind2] = sort(Fit2(:,2));
        Xbest2 = HM2(sind2(1),:);
    elseif R4 == 3
        [SF3, sind3] = sort(Fit3(:,3));
         Xbest3 = HM3(sind3(1),:);
    else
        [SF4, sind4] = sort(Fit4(:,4));
        Xbest4 = HM4(sind4(1),:);
    end

     i=1;
     while i<=dim_epi
         if rand<HMCR
               %sL = length(sind1(:,1));
                 %% ???? ?? ?????????
                 a = 1 + abs( ceil(normrnd(0,HMS/3,1)));
                 while a > HMS
                     a = 1 + abs(ceil(normrnd(0,HMS/3,1)));
                 end
             %%

              if R4 == 1
                  Xnew(i)=HM(sind1(a) ,i);
              elseif R4 == 2
                  Xnew(i)=HM2(sind2(a) ,i);
              elseif R4 == 3
                  Xnew(i)=HM3(sind3(a) ,i);
              else
                  Xnew(i)=HM4(sind4(a) ,i);
              end
              if rand < PAR
                      rh = ceil(rand(1,3)*HMS);
                      while rh(2) == rh(3)
                          rh(3) = ceil(rand*HMS);
                      end
                      if R4 == 1
                          Xnew(i)=HM(rh(1),i)+ 2*(rand-0.5)*(HM(rh(2),i) - HM(rh(3),i));
                      elseif R4 == 2
                          Xnew(i)=HM2(rh(1),i)+ 2*(rand-0.5)*(HM2(rh(2),i) - HM2(rh(3),i));
                      elseif R4 == 3
                          Xnew(i)=HM3(rh(1),i)+ 2*(rand-0.5)*(HM3(rh(2),i) - HM3(rh(3),i));
                      else
                          Xnew(i)=HM4(rh(1),i)+ 2*(rand-0.5)*(HM4(rh(2),i) - HM4(rh(3),i));
                      end

                        Xnew(i) = round(Xnew(i));
                        Xnew(i)=max(min(Xnew(i),SNPs),1);

                end
          else
                Xnew(i)=ceil(rand*SNPs);
          end
             %% ??????
             cc = 0;
          while i>1 && ismember(Xnew(i),Xnew(1:i-1))

              if R4 == 1
                  rr = ceil(rand*HMS);
                Xnew(i)=HM(rr,i); %ceil(rand*SNPs);
              elseif R4==2
                  rr = ceil(rand*HMS);
                 Xnew(i)=HM2(rr,i);
              elseif R4 == 3
                  rr = ceil(rand*HMS);
                 Xnew(i)=HM3(rr,i);
              else
                  rr = ceil(rand*HMS);
                  Xnew(i)=HM4(rr,i);
              end
              cc = cc + 1;
              if cc > 2
                  Xnew(i) =  ceil(rand*SNPs);
              end

          end
          i=i+1;


     end

      %% ??????
              Xtemp=Xnew;
              Xnew=sort(Xnew);
              c2 = 0;
              while ( ismember(Xnew,X,"rows") || ismember(Xnew,X2,"rows") || ismember(Xnew,X3,"rows")|| ismember(Xnew,X4,"rows")||~isempty(Center) && getMinDistance2(Xnew,Center(:,1:dim_epi),dim_epi)< dim_epi - 1)
                 J=ceil(rand*dim_epi);
                  r=ceil(rand*SNPs);
                  while ismember(r,Xnew)
                      r=ceil(rand*SNPs);
                  end
                  Xnew(J)=r;
                  Xtemp=Xnew;
                 Xnew=sort(Xnew);
                 c2 = c2 + 1;
%                  if c2 > 2
%                      break;
%                  end
              end
       %%

%      if length(unique(Xnew))<dim_epi
%          fprintf("%d, rep,rep!!!!!!!!!!!\n",Xnew);
%      end



  [score,score2,score3,score4] = multi_criteriaEvaluationFuns3(data(:,Xnew),State);

   %%
   Flag = 0;

   [fworst,idworst] = max(Fit(:,1));
   [fworst2,idworst2] = max(Fit2(:,2));
   [fworst3,idworst3] = max(Fit3(:,3));
   [fworst4,idworst4] = max(Fit4(:,4));
        if score<fworst || score2<fworst2 || score3<fworst3 || score4<fworst4
            if  score<=fworst || rand < 0.2*exp(-(score - fworst)/T1)
                Fit(idworst,:)=[score,score2,score3,score4];
                X(idworst,:)=Xnew;
                HM(idworst,:)=Xtemp;
                Flag = Flag + 1000;

            end
            if score2<=fworst2 || rand <0.2*exp(-(score2 - fworst2)/T2)
               Fit2(idworst2,:)=[score,score2,score3,score4];
               X2(idworst2,:)=Xnew;
               HM2(idworst2,:)=Xtemp;
                Flag = Flag + 100;

            end

            if score3 <= fworst3 || rand < 0.2*exp(-(score3 - fworst3)/T3)
               Fit3(idworst3,:)=[score,score2, score3,score4];
               X3(idworst3,:)=Xnew;
               HM3(idworst3,:)=Xtemp;
                Flag = Flag + 10;
            end
             if score4 <= fworst4 || rand < 0.2*exp(-(score4 - fworst4)/T4)
               Fit4(idworst4,:)=[score,score2, score3,score4];
               X4(idworst4,:)=Xnew;
               HM4(idworst4,:)=Xtemp;
                Flag = Flag + 1;
            end
        end
       NC=NC+1;
 %%  The program is terminted if the Xnew is the solution.
% flag = -1;
%    if Xnew == CX
%           canSize1 = canSize1+1;
%           Candidate1(canSize1,:) = [CX,score,score2,score3,score4];
%           if Flag > 0
%               flag = Flag;
%           else
%               flag = 1111;
%           end
%
%             break;
%    end
  %% ??????
  if flag > 0 || mod(NC,maxIterForLocalSearch)==0
                  for i=1:HMS
                      if isempty(Elite1)  %% ???2?
                          Elite1 = X(i,:);
                          Efit1 = Fit(i,:);

                          Elite2 = X2(i,:);
                          Efit2 = Fit2(i,:);

                          Elite3 = X3(i,:);
                          Efit3 = Fit3(i,:);

                          Elite4 = X4(i,:);
                          Efit4 = Fit4(i,:);


                      else
                              if length(Elite1(:,1))<EliteSize
                                  Elite1=[Elite1;X(i,:)];
                                  Efit1=[Efit1;Fit(i,:)];

                                  Elite2=[Elite2;X2(i,:)];
                                  Efit2=[Efit2;Fit2(i,:)];

                                  Elite3=[Elite3;X3(i,:)];
                                  Efit3=[Efit3;Fit3(i,:)];

                                  Elite4=[Elite4;X4(i,:)];
                                  Efit4=[Efit4;Fit4(i,:)];

                              else
                                %%  ?????????
                                      [~,eidworst]=max(Efit1(:,1));
                                      [~,eidworst2]=max(Efit2(:,2));
                                      [~,eidworst3]=max(Efit3(:,3));
                                      [~,eidworst4]=max(Efit4(:,4));
%                                       [size(Efit1),size(Fit)]
                                          if Efit1(eidworst,1)> Fit(i,1) %%&& Efit(eidworst,2)> Fit(i,2)
                                              Elite1(eidworst,:)=X(i,:);
%                                               size(Elite1)
%                                               eidworst
                                              Efit1(eidworst,:)=Fit(i,:);
                                          end
                                          if Efit2(eidworst2,2)> Fit2(i,2)
                                              Elite2(eidworst2,:)=X2(i,:);
                                              Efit2(eidworst2,:)=Fit2(i,:);
                                          end

                                          if Efit3(eidworst3,3)> Fit3(i,3)
                                              Elite3(eidworst3,:)=X3(i,:);
                                              Efit3(eidworst3,:)=Fit3(i,:);
                                          end

                                          if Efit4(eidworst4,4)> Fit4(i,4)
                                              Elite4(eidworst4,:)=X4(i,:);
                                              Efit4(eidworst4,:)=Fit4(i,:);
                                          end


                              end
                      end

                  end



      %% ??????
      [~,idebest1]=min(Efit1(:,1));
      [~,idebest2]=min(Efit2(:,2));
      [~,idebest3]=min(Efit3(:,3));
      [~,idebest4]=min(Efit4(:,4));

      E1=Elite1(idebest1,:);
      E2=Elite2(idebest2,:);
      E3=Elite3(idebest3,:);
      E4=Elite4(idebest4,:);

          CCC1=Elite1(idebest1,:);
          CCC2=Elite2(idebest2,:);
          CCC3=Elite3(idebest3,:);
          CCC4=Elite4(idebest4,:);

          minDist1=getMinDistance2(CCC1,Elite1,dim_epi) ;  %% ??????,?????
           Center=[Center;[CCC1,minDist1]];
%            disp(Center)

   %%
        if CCC1(1,:)~=CCC2(1,:)
           minDist2=getMinDistance2(CCC2,Elite2,dim_epi) ;
           Center=[Center;[CCC2,minDist2]];
        end


        if CCC3(1,:) ~= CCC1(1,:) & CCC3(1,:) ~= CCC2(1,:)
             minDist3=getMinDistance2(CCC3(1,:),Elite3,dim_epi) ;
           Center=[Center;[CCC3(1,:),minDist3]];
        end

        if CCC4(1,:) ~= CCC1(1,:) & CCC4(1,:) ~= CCC2(1,:) & CCC4(1,:) ~= CCC3(1,:)
             minDist4=getMinDistance2(CCC4(1,:),Elite4,dim_epi) ;
           Center=[Center;[CCC4(1,:),minDist4]];
        end


          Alen = min([EliteSize,length(Elite1(:,1)),length(Elite2(:,1)),length(Elite3(:,1)),length(Elite4(:,1))]);
           for i=1:Alen
                      if canSize1==0
                          canSize1=canSize1+1;
                          Candidate1(canSize1,:)=[Elite1(i,:),Efit1(i,:)]
                      else
                          if ~ismember(Elite1(i,:),Candidate1(1:canSize1,1:dim_epi),"rows")
                              if canSize1<CandidateSize
                                  canSize1=canSize1+1;
                                  Candidate1(canSize1,:)=[Elite1(i,:),Efit1(i,:)]
                              else
                                  [Fitworst,Cind]=max(Candidate1(:,dim_epi+1));  %% dim_epi+1 ???1? dim_epi+2??2
                                  if Fitworst>Efit1(i,1)
                                      Candidate1(Cind,:)=[Elite1(i,:),Efit1(i,:)]
                                  end
                              end
                          end
                      end

                    if canSize2==0
                          canSize2=canSize2+1;
                          Candidate2(canSize2,:)=[Elite2(i,:),Efit2(i,:)]
                    else
                          if ~ismember(Elite2(i,:),Candidate2(1:canSize2,1:dim_epi),"rows")
                              if canSize2<CandidateSize
                                  canSize2=canSize2+1;
                                  Candidate2(canSize2,:)=[Elite2(i,:),Efit2(i,:)]
                              else
                                  [Fitworst,Cind]=max(Candidate2(:,dim_epi+2));
                                  if Fitworst>Efit2(i,2)
                                      Candidate2(Cind,:)=[Elite2(i,:),Efit2(i,:)]
                                  end
                              end
                          end
                    end

                     if canSize3==0
                          canSize3=canSize3+1;
                          Candidate3(canSize3,:)=[Elite3(i,:),Efit3(i,:)]
                     else

                          if ~ismember(Elite3(i,:),Candidate3(1:canSize3,1:dim_epi),"rows")
                              if canSize3<CandidateSize
                                  canSize3=canSize3+1;
                                  Candidate3(canSize3,:)=[Elite3(i,:),Efit3(i,:)]
                              else
                                  [Fitworst,Cind]=max(Candidate3(:,dim_epi+3));
                                  if Fitworst>Efit3(i,3)
                                      Candidate3(Cind,:)=[Elite3(i,:),Efit3(i,:)]
                                  end
                              end
                          end
                     end


                     if canSize4==0
                          canSize4=canSize4+1;
                          Candidate4(canSize4,:)=[Elite4(i,:),Efit4(i,:)]
                     else

                          if ~ismember(Elite4(i,:),Candidate4(1:canSize4,1:dim_epi),"rows")
                              if canSize4<CandidateSize
                                  canSize4=canSize4+1;
                                  Candidate4(canSize4,:)=[Elite4(i,:),Efit4(i,:)]
                              else
                                  [Fitworst,Cind]=max(Candidate4(:,dim_epi+4));
                                  if Fitworst>Efit4(i,4)
                                      Candidate4(Cind,:)=[Elite4(i,:),Efit4(i,:)]
                                  end
                              end
                          end
                      end
           end
%            Elite1
% Efit1
%
% Elite2
% Efit2
%
% Elite3
% Efit3
%
% Elite4
% Efit4

%     Candidate1(:,1:dim_epi)
    % Center
           Elite1=[];
           Efit1=[];

           Elite2=[];
           Efit2=[];

           Elite3=[];
           Efit3=[];

           Elite4=[];
           Efit4=[];
           %%???????
          %% ?????
       fprintf("r!  ");
%           fprintf("initHM");
          [X,HM,Fit]=InitHM102(data,HMS,dim_epi,Center);

         X2=X;
         HM2=HM;
         Fit2=Fit;

          X3=X;
         HM3=HM;
         Fit3=Fit;
         X4=X;
         HM4=HM;
         Fit4=Fit;

%          if ismember(CX,X,"rows")
%               canSize1 = canSize1+1;
%               Candidate1(canSize1,:) = [CX,score,score2,score3,score4];
%               if Flag > 0
%                   flag = 1111;%Flag;
%                    break;
%               elseif flag<0
%                   flag = 0;
%               end
%          end
%          fprintf("end Init");
  end

%   if flag > 0
%       break;
%   end
end
Candidate1
Candidate2
Candidate3
Candidate4
Candidate=[Candidate1(1:canSize1,:);Candidate2(1:canSize2,:);Candidate3(1:canSize3,:);Candidate4(1:canSize4,:)];

Candidate = unique(Candidate,'rows')
canSize = length(Candidate(:,1));
totaltime=toc;
%
