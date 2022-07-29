%%
clear
close all


DATADIR=" dir to the dataset folder";
load([DATADIR 'DATA_TUS_All_02.mat'])
FUS_DS=readtable([DATADIR 'TUS_DataCollection_BothMonkey.xlsx']); %reading excel sheet

cd(DATADIR)

DimFLabel={'Shape','Pattern','Color','Arm'};
FUSCnd_Label={'S-ACC','H-ACC','S-aSTR','H-aSTR'};
DVt=[1,2,3,5]; %dimensions

%%
Trialdata = FLUData;
EventData = EV;
if 0
Event.OBjects={Trialdata.Object1,Trialdata.Object2,Trialdata.Object3};
Event.OBjectsLoc={Trialdata.Object1loc,Trialdata.Object2loc,Trialdata.Object3loc};
for i=1:length(EventData.Trial)
Event.Time_fix(i,:)=[EventData.Trial{i, 1}.Time_Fix_Obj1,EventData.Trial{i, 1}.Time_Fix_Obj2,EventData.Trial{i, 1}.Time_Fix_Obj3,EventData.Trial{i, 1}.Time_Fix_TokenBar,EventData.Trial{i, 1}.Time_Fix_ChosenObj];
Event.Dilation(i,:)=nanmean(EventData.Trial{i, 1}.Event_c.Dilation(logical(EventData.Trial{i, 1}.Event_c.ChoiceEvent)));  
Event.ChosenObjGaze(i,:)=Event.Time_fix(i,Trialdata.ChosenObject(i));
Event.ObjFeat(i,:)=Event.OBjects{1,Trialdata.ChosenObject(i)}(i,:);
Event.ObjLoc(i,:)=Event.OBjectsLoc{1,Trialdata.ChosenObject(i)}(i,:);
Event.Num_fix(i,:)=[EventData.Trial{i, 1}.Num_Fix_Obj1,EventData.Trial{i, 1}.Num_Fix_Obj2,EventData.Trial{i, 1}.Num_Fix_Obj3,EventData.Trial{i, 1}.Num_Fix_TokenBar];
Event.NumFixChosenObjGaze(i,:)=Event.Num_fix(i,Trialdata.ChosenObject(i));
Event.After_FB_Time_fix(i,:)=[EventData.Trial{i, 1}.After_FB_Time_Fix_Obj1,EventData.Trial{i, 1}.After_FB_Time_Fix_Obj2,EventData.Trial{i, 1}.After_FB_Time_Fix_Obj3,EventData.Trial{i, 1}.After_FB_Time_Fix_TokenBar];
Event.After_FB_ChosenObjGaze(i,:)=Event.After_FB_Time_fix(i,Trialdata.ChosenObject(i));
end

Event.Explore=nansum(Event.Time_fix(:,1:3),2);
Event.Exploite=nansum(Event.Time_fix(:,5),2);
Event.TokenGaze=nansum(Event.Time_fix(:,4),2);
end
%%
for j=1:length(Trialdata.BlockLabel),
    if ~isempty(findstr(Trialdata.BlockLabel{j},'0_irrel')), iCndDim(j)=1; iCndDimLabel{j} = '1Dim'; end
    if ~isempty(findstr(Trialdata.BlockLabel{j},'1_irrel')), iCndDim(j)=2; iCndDimLabel{j} = '2Dim'; end
    if ~isempty(findstr(Trialdata.BlockLabel{j},'2_irrel')), iCndDim(j)=3; iCndDimLabel{j} = '3Dim'; end
    
    if ~isempty(findstr(Trialdata.BlockLabel{j},'2.100_Inc.-1.100')), iCndTok(j)=1; iCndTokLabel{j} = '1G0L'; end
    if ~isempty(findstr(Trialdata.BlockLabel{j},'3.100_Inc.-0.100')), iCndTok(j)=2; iCndTokLabel{j} = '2G3L'; end
    
    if ~isempty(findstr(Trialdata.BlockLabel{j},'Intra')), iCndSwitch(j)=1; iCndSwitchLabel{j} = 'Intra'; end
    if ~isempty(findstr(Trialdata.BlockLabel{j},'Extra')), iCndSwitch(j)=2; iCndSwitchLabel{j} = 'Extra'; end
    
end

SessN=nan(size(Trialdata.DatasetName));
FUS_DS=readtable([DIR.Sheet 'TUS_DataCollection_all_06.xlsx']); %reading excel sheet
FUS=[];
Sess_Date=FUS_DS.Date;
FUS_Sess=FUS_DS.Session;
for i=1:length(Sess_Date)    
Vsdate=contains(Trialdata.DatasetName,Sess_Date{i,1});
SessN(Vsdate)=FUS_Sess(i);
ACC_Sess(i)=nanmean(Trialdata.Accuracy(Vsdate));
end

TIB=Trialdata.TrialInBlock;
NumT=length(MnkID);
OC=Trialdata.Outcome;
Blkbrd=find(diff(TIB)<0);
Blkbrd=[0; Blkbrd; NumT];
FullBarr=zeros(length(OC),1);
FeatDim=Trialdata.dimensionVector(Trialdata.TargetFeature);
full_barr=find(strcmp(Trialdata.AllTokensCompleted,'True'));
FullBarr(full_barr)=1;
Rt=Trialdata.ReactionTime;
Rt(Rt>1 | Rt<.075)=nan;
Tocc=nan(1,NumT);
NumBlk=length(Blkbrd)-1;
NumSess=length(unique(SessN));
ITI=diff(Trialdata.ITI);
ITI=(Trialdata.ITI);
% ITI=[0;ITI];
ITI(ITI<0)=0;


FUS.Inc=zeros(NumT,1);
FUS.TIE=nan(NumT,1);
FUS.Cnd_Label=cell(NumT,1); %% High Low Sham
FUS.Cnd=nan(NumT,1); %High=1 Low=2 Sham=3
FUS.PP=nan(NumT,1); %% pre-stim=1 post-stim=2
FUS.Time=Trialdata.ITI; %% time from stim
FUS.Area=cell(NumT,1); %% targeted area ACC, aSTR
U_Sess=unique(SessN);
V_Sess=FUS_DS.Session(~isnan(FUS_DS.Session));

for i=1:length(V_Sess)
    
        fus_sess=find(SessN==V_Sess(i));
        if ~isempty(fus_sess)
            FUS.TIE(fus_sess,1)=1:length(fus_sess);
            FUS.Cnd_Label(fus_sess,1)=FUS_DS.Condition(i);
            FUS.Area(fus_sess,1)=FUS_DS.Area(i);
            trialp=find(SessN==V_Sess(i) & TIB==FUS_DS.Trial_resume(i) & Trialdata.BlockNum==FUS_DS.Block_resume(i));
            Extrial=find(SessN==V_Sess(i)  & (Trialdata.BlockNum>FUS_DS.End_block(i) | FUS_DS.Inclusion(i)==0));
            inc=find(SessN==V_Sess(i) & Trialdata.BlockNum<=FUS_DS.End_block(i));
            FUS.Inc(inc,1)=1;
            TIEp=FUS.TIE(trialp);
%             FUS.Inc(fus_sess,1)=1:length(fus_sess);
            FUS.Inc(Extrial,1)=0;
            pre_stim=find(FUS.TIE<TIEp  & SessN==V_Sess(i) );
            post_stim=find(FUS.TIE>=TIEp & SessN==V_Sess(i) );
            FUS.PP(pre_stim,1)=1;
            FUS.PP(post_stim,1)=2;
            FUS.Time(fus_sess,1)=FUS.Time(fus_sess,1)-Trialdata.ITI(TIEp);
        end
end
FUS.AreaCnd=nan(size(FUS.Cnd));
FUS.Cnd(find(strcmp(FUS.Cnd_Label,'High')),1)=1;
FUS.Cnd(find(strcmp(FUS.Cnd_Label,'Sham')),1)=0;

FUS.AreaCnd(find(strcmp(FUS.Area,'ACC')),1)=1;
FUS.AreaCnd(find(strcmp(FUS.Area,'aSTR')),1)=2;

FUS.Cnd4=nan(size(FUS.Cnd));
FUS.Cnd4(FUS.AreaCnd==1 & FUS.Cnd==0)=1;
FUS.Cnd4(FUS.AreaCnd==1 & FUS.Cnd==1)=2;
FUS.Cnd4(FUS.AreaCnd==2 & FUS.Cnd==0)=3;
FUS.Cnd4(FUS.AreaCnd==2 & FUS.Cnd==1)=4;




%% Getting Block data
col={[.8,.1,.1], [.1,.1,.8],[.1,.8,.1]};
Col4Cnd=[.8 .5 .5; 1 .1 .1 ;.5 .8 .5; .1 1 .1];
clear BlockDATA    
for i=1:NumBlk
    T_occ(Blkbrd(i)+1:Blkbrd(i+1))=mod(cumsum([0; OC(Blkbrd(i)+1:Blkbrd(i+1)-1)]),5);
 
    Acc_OC=Trialdata.Accuracy(Blkbrd(i)+1:Blkbrd(i+1));
    CONVAC=conv(Acc_OC,ones(1,12))/12;
    FF_Acc=CONVAC(12:end-11);

    learning_point=find(FF_Acc>.80 ,1,'first');
    if learning_point<1
        learning_point=nan;
        BlockDATA.Asymp(i,1)=nan;
    end
    if ~isempty(learning_point) && learning_point<60 && ~isnan(learning_point)
        
        BlockDATA.LP(i,1)=learning_point;
        BlockDATA.Asymp(i,1)=mean(Trialdata.Accuracy(learning_point+Blkbrd(i):Blkbrd(i+1)));
        BlockDATA.EplrPr(i,1)=(nanmean(nansum(Event.Time_fix(Blkbrd(i)+1:learning_point+Blkbrd(i),1:3),2)));
        BlockDATA.EplrPs(i,1)=(nanmean(nansum(Event.Time_fix(learning_point+Blkbrd(i):Blkbrd(i+1),1:3),2)));
    elseif isempty(learning_point)
        blkL=Blkbrd(i+1)-Blkbrd(i);
        mDl =fitlm(FF_Acc(end-12:end)',[blkL-12:blkL]);
        lp=round(predict(mDl,.8));
        BlockDATA.EplrPr(i,1)=(nanmean(nansum(Event.Time_fix(Blkbrd(i)+1:Blkbrd(i+1),1:3),2)));
        BlockDATA.EplrPs(i,1)=nan;
        BlockDATA.LP(i,1)=lp;
        BlockDATA.Asymp(i,1)=nan;%nanmean(Trialdata.Accuracy(Blkbrd(i+1)-12:Blkbrd(i+1)));
    else
        BlockDATA.LP(i,1)=nan;
        BlockDATA.Asymp(i,1)=nan; 
    end
    BlockDATA.Length(i,1)=Blkbrd(i+1)-Blkbrd(i)+1;
    BlockDATA.SessionNum(i,1)=mean(SessN(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.BlockNum(i,1)=mean(Trialdata.BlockNum(Blkbrd(i)+1:Blkbrd(i+1)));
    BlkInSess=find(V_Sess==BlockDATA.SessionNum(i,1));
    BlockDATA.Inc(i,1)=(BlockDATA.BlockNum(i,1)<= FUS_DS.End_block(BlkInSess) && FUS_DS.Inclusion(BlkInSess));
    BlockDATA.RT(i,1)=nanmean(Rt(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.DimCond(i,1)=mean(iCndDim(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.DimN(i,1)=Trialdata.dimensionVector((mean(Trialdata.TargetFeature(Blkbrd(i)+1:Blkbrd(i+1),:),1)));
    BlockDATA.RR(i,1)=nansum(FullBarr(Blkbrd(i)+1:Blkbrd(i+1)))/(Blkbrd(i+1)-Blkbrd(i));
    BlockDATA.TokCond(i,1)=mean(iCndTok(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.SwitchCond(i,1)=mean(iCndSwitch(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.Accuracy(i,1)=mean(Trialdata.Accuracy(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.ITI(i,1)=mean(ITI(Blkbrd(i)+1:Blkbrd(i+1)));
    BlockDATA.Mnk(i,1)=mean(MnkID(Blkbrd(i)+1:Blkbrd(i+1)));
   
    BlockDATA.FUStime(i,1)=(nanmean(FUS.Time(Blkbrd(i)+1:Blkbrd(i+1))))/60;
    BlockDATA.FUSCnd(i,1)=(nanmean(FUS.Cnd(Blkbrd(i)+1:Blkbrd(i+1))));
    BlockDATA.FUSCnd4(i,1)=(nanmean(FUS.Cnd4(Blkbrd(i)+1:Blkbrd(i+1))));
    BlockDATA.FUSAreaCnd(i,1)=(nanmean(FUS.AreaCnd(Blkbrd(i)+1:Blkbrd(i+1))));
    
    BlockDATA.Eplr(i,1)=(nanmean(nansum(Event.Time_fix(Blkbrd(i)+1:Blkbrd(i+1),1:3),2)));
    BlockDATA.Eplt(i,1)=(nanmean(nansum(Event.Time_fix(Blkbrd(i)+1:Blkbrd(i+1),5),2)));
    BlockDATA.ObjG(i,1)=(nanmean((Event.ChosenObjGaze(Blkbrd(i)+1:Blkbrd(i+1),1))));
    BlockDATA.TokG(i,1)=(nanmean((Event.TokenGaze(Blkbrd(i)+1:Blkbrd(i+1),1))));
    BlockDATA.AFB_Eplr(i,1)=(nanmean(nansum(Event.After_FB_Time_fix(Blkbrd(i)+1:Blkbrd(i+1),1:3),2)));
    BlockDATA.AFB_ObjG(i,1)=(nanmean((Event.After_FB_ChosenObjGaze(Blkbrd(i)+1:Blkbrd(i+1),1))));
    BlockDATA.AFB_TokG(i,1)=(nanmean((Event.After_FB_Time_fix(Blkbrd(i)+1:Blkbrd(i+1),4))));
       
end

BlockDATA.StimCnd=nan(NumBlk,1); % pre-stim=1 , post-stim=2
BlockDATA.StimTCnd=nan(NumBlk,1); % pre-stim=1 , post-stim=2

for i=1:length(U_Sess)
    sess_id=find(V_Sess==U_Sess(i));
    blkp=FUS_DS.Block_pause(sess_id);
    blkr=FUS_DS.Block_resume(sess_id);
    trialr=FUS_DS.Trial_resume(sess_id);
    pre_stim_blk=find(BlockDATA.SessionNum==U_Sess(i) & BlockDATA.BlockNum<=blkp);
    post_stim_blk=find(BlockDATA.SessionNum==U_Sess(i) & BlockDATA.BlockNum>=blkr);
    BlockDATA.StimTCnd(pre_stim_blk,1)=1;
    BlockDATA.StimTCnd(post_stim_blk,1)=2;
end


   


%% Getting the session data
pind=[3,2,1];
YLabel={'Learning Point','Exploration','Asymptote'}
Title={'Both','Igor','Wotan'};
% close all
figure
for Mnki=1:3
    
Sess=[];
Sess.SeshNum=V_Sess;
for i=1:length(V_Sess)
    for j=1:4
    idn=find(BlockDATA.SessionNum==V_Sess(i));
    idnBs=find(BlockDATA.SessionNum==V_Sess(i) & BlockDATA.StimTCnd==1 & BlockDATA.FUSCnd4==j & BlockDATA.TokCond~=2 & BlockDATA.BlockNum<36 & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc);
    idnTUS=find(BlockDATA.SessionNum==V_Sess(i) & BlockDATA.StimTCnd==2 & BlockDATA.FUSCnd4==j & BlockDATA.TokCond~=2 & BlockDATA.BlockNum<36 & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc);
    Sess.LP_Bs(i,j)=nanmean(BlockDATA.LP(idnBs));
    Sess.LP_TUS(i,j)=nanmean(BlockDATA.LP(idnTUS));
    Sess.Eplr_Bs(i,j)=nanmean(BlockDATA.Eplr(idnBs));
    Sess.Eplr_TUS(i,j)=nanmean(BlockDATA.Eplr(idnTUS));
    Sess.Asymp_Bs(i,j)=nanmean(BlockDATA.Asymp(idnBs));
    Sess.Asymp_TUS(i,j)=nanmean(BlockDATA.Asymp(idnTUS));
    end
end
Sess.LP_Bs(isempty(Sess.LP_Bs))=nan;
Sess.LP_TUS(isempty(Sess.LP_Bs))=nan;
Sess.Eplr_Bs(isempty(Sess.LP_Bs))=nan;
Sess.Eplr_TUS(isempty(Sess.LP_Bs))=nan;
Sess.Asymp_Bs(isempty(Sess.LP_Bs))=nan;
Sess.Asymp_TUS(isempty(Sess.LP_Bs))=nan;
VecSess{1,1}=Sess.LP_Bs;
VecSess{2,1}=Sess.Eplr_Bs;
VecSess{3,1}=Sess.Asymp_Bs;
VecSess{1,2}=Sess.LP_TUS;
VecSess{2,2}=Sess.Eplr_TUS;
VecSess{3,2}=Sess.Asymp_TUS;

end


%% plot block resolved figs
YLabel={'Learning Point','Exploration','Asymptote','Reaction time'}
Title={'Both','Igor','Wotan'};
close all
clear Blk
figure
pind=[3,2,1];
arr1=[1,6:6:24];
arr2=[6,12:6:30];
arr1=[0,.25,.5,.75];
arr2=[.25,1,1.25,1.75];
l0 = length(arr1);
for Mnki=1:3
    
    for j=1:4
        idn=find(BlockDATA.BlockNum==1 & BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc);
        Blk.Vec{1,1}{Mnki,j}=nan(length(idn),l0);
        Blk.Vec{1,2}{Mnki,j}=nan(length(idn),l0);
        Blk.Vec{1,3}{Mnki,j}=nan(length(idn),l0);
        Blk.Vec{1,4}{Mnki,j}=nan(length(idn),l0);
        for i=1:l0
            idn=find((BlockDATA.ITI/3600)>=arr1(i) & (BlockDATA.ITI/3600)<arr2(i) & BlockDATA.TokCond~=2 & BlockDATA.FUSCnd4==j  & BlockDATA.BlockNum<36 & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc)
            l1=length(idn);
            
            Blk.Vec{1,1}{Mnki,j}(1:l1,i)=BlockDATA.LP(idn);
            Blk.Vec{1,2}{Mnki,j}(1:l1,i)=BlockDATA.Eplr(idn);
            Blk.Vec{1,3}{Mnki,j}(1:l1,i)=BlockDATA.Asymp(idn);
            Blk.Vec{1,4}{Mnki,j}(1:l1,i)=BlockDATA.RT(idn);
            
            Blk.Mn{1,1}(j,i)=nanmean(BlockDATA.LP(idn));
            Blk.Mn{1,2}(j,i)=nanmean(BlockDATA.Eplr(idn));
            Blk.Mn{1,3}(j,i)=nanmean(BlockDATA.Asymp(idn));
            Blk.Mn{1,4}(j,i)=nanmean(BlockDATA.RT(idn));
            
            Blk.Se{1,1}(j,i)=nanstd(BlockDATA.LP(idn))./sqrt(l1);
            Blk.Se{1,2}(j,i)=nanstd(BlockDATA.Eplr(idn))./sqrt(l1);
            Blk.Se{1,3}(j,i)=nanstd(BlockDATA.Asymp(idn))./sqrt(l1);
            Blk.Se{1,4}(j,i)=nanstd(BlockDATA.RT(idn))./sqrt(l1);
            
        end
    end
  
        for v=1:4
        P_ind = sub2ind([ 3 4],Mnki,v);
        subplot(4,3,P_ind)
        hold on
        for k=1:4
            x=1:l0;
        Mn=Blk.Mn{1,v}(k,:);    
        SE=Blk.Se{1,v}(k,:);

        shadedErrorBar(x,Mn,SE,{'color',squeeze(Col_A(3,k,:))},1)
        end
        if Mnki==1
        ylabel(YLabel{v})
        end
        if v==1
        title(['Monkey:  ', Title{Mnki}])
        end
        set(gca,'tickdir','out')
    end
    
end 


%% Trial data 


Nt=5;
dTIB=find([0;diff(TIB)]<0);
NTIB=TIB;
dTIB=reshape(dTIB-repmat(1:Nt,length(dTIB),1),Nt*length(dTIB),1);
NTIB(dTIB)=nan;
Vec{1}=Trialdata.Accuracy;
Vec{2}=Event.Explore;
Vec{3}=Event.TokenGaze;
Vec{4}=Event.ChosenObjGaze;

TUS_Col = [.8 .5 .5; 1 .1 .1; .5 .8 .5; .1 1 .1];
YLabel={'Accuracy','Exploration','reaction time','NumFix'};
ocArr=[0, 0, 1, 1];
PrevTrl_Outcome=[0;Trialdata.Accuracy(1:end-1)];
PrevTrl_Outcome2=[0;0;Trialdata.Accuracy(1:end-2)];
close all
figure 
for p=1:4
for i=1:4

        for k=1:4


    Y1=find( ~isnan(NTIB) & PrevTrl_Outcome==ocArr(p)  & FUS.PP==2 & FUS.Cnd4==k & FUS.Inc & TokCnd==TokArr(p) & MnkID~=3 & iCndDim'>0);
    Y0=find( ~isnan(NTIB) & PrevTrl_Outcome==ocArr(p)  & FUS.PP==1 & FUS.Cnd4==k & FUS.Inc & TokCnd==TokArr(p) & MnkID~=3 & iCndDim'>0);          

            
            for nt=1:Nt
                if nt==1
                    Xvec  = Vec{i}(Y1);
                    Xvec0 = Vec{i}(Y0);
                else
                    Xvec  = [Xvec,Vec{i}(Y1+nt-1)];
                    Xvec0 = [Xvec0,Vec{i}(Y0+nt-1)];
                end
            end
            sd = nanstd([Xvec])/sqrt(length(Y1));
            
            Md0 = nanmedian(Xvec0);
            Mn0 = nanmean(Xvec0);
            sd0 = nanstd(Xvec0);%/sqrt(length(Y0));
            Mn2=(Xvec-Mn0)./sd0;
            Mn=nanmean(Mn2);
            sd=nanstd(Mn2)./sqrt(length(Y1));
            P_ind = sub2ind([ NCond NPar],p,i);
            subplot(NPar,NCond, P_ind)
            hold on

            
%             if j<3
%             scatter(x, Mn,12,'marker', Marker{j}, 'markerfacecolor',lCol,'markeredgecolor',lCol)
%             scatter(x0, Mn0,12,'marker', Marker{j}, 'markerfacecolor',lCol,'markeredgecolor',lCol)
%             
%             errorbar(Mn0, sd0, 'color',lCol,'linewidth',LW(j),'linestyle',':')
            errorbar(Mn, sd, 'color',TUS_Col(k,:),'linewidth',2,'linestyle','-')
            
%             end
        end
        
    set(gca,'tickdir','out','xlim',[0 Nt+1],'xtick',[1:Nt],'ylim',[-.3 .3])
    
    if p==1
    ylabel(YLabel{i})
    end
    if i==1
        title(['Prev-Trial Outcome:  ', Title{p}])
    end
   end 
end


%%
Nt=4;
dTIB=find([0;diff(TIB)]<0);
NTIB=TIB;
dTIB=reshape(dTIB+repmat(1:Nt,length(dTIB),1),Nt*length(dTIB),1);
NTIB(dTIB)=nan;
Vec{1}=Trialdata.Accuracy;
Vec{2}=Event.Explore;

Vec{3}=Event.TokenGaze;
% Vec{4}(Vec{4}>.5)=nan;
% Vec{3}(Vec{3}>0.25)=nan;
% Vec{5}=BlockDATA.EplrPr;
Vec{4}=Event.ChosenObjGaze;
YLabel={'Accuracy','Exploration','reaction time','NumFix'};
ocArr=[-1, 0, 2, 3];
GrossOC=[];
for i=1:Nt
    if i==1;
        GrossOC=[0;OC(1:end-1)];
    else
       GrossOC=[GrossOC,[zeros(i,1);OC(1:end-i)]];
    end
end

GrossOC=nansum(GrossOC,2);
GrArr=unique(GrossOC);
GrArr(GrArr<-3 | GrArr>12)=[];
close all
figure
Mn=[];
sd=[];
for iMnk=1:3
for p=1
    subplot(3,1,iMnk)
    hold on
    for k=1:4
        for g=1:length(GrArr)-5
         Y1=find( ~isnan(NTIB) & GrossOC>=GrArr(g) & GrossOC<GrArr(g+4) & iCndTok'==1 & iCndDim'>0 & FUS.PP==2 & FUS.Cnd4==k & FUS.Inc & MnkID~=iMnk );
         Y0=find( ~isnan(NTIB) & GrossOC>=GrArr(g) & GrossOC<GrArr(g+4) & iCndTok'==1 &  iCndDim'>0 & FUS.PP==1  & FUS.Cnd4==k  & FUS.Inc & MnkID~=iMnk);
         
         Xvec  = Vec{p}(Y1);
         Xvec0 = Vec{p}(Y0);
         
         Mn0 = nanmean(Xvec0);
         sd0 = nanstd(Xvec0); 
         
         Mn2=(Xvec-Mn0)./sd0;
         Mn(g)=nanmean(Mn2);
         sd(g)=nanstd(Mn2)./sqrt(length(Y1));
        end
        lCol=squeeze(Col_A(3,k,:));
        errorbar(Mn, sd, 'color',lCol,'linewidth',2,'linestyle','-')
    end
    set(gca,'tickdir','out','xlim',[0 length(GrArr)-4],'xtick',[1: length(GrArr)-5],'xticklabel',GrArr(1:end-4),'ylim',[-.3 .3])
    ylabel(YLabel{p})
    xlabel(['Gross Token Income =< x'])
end
end           
     
    

 