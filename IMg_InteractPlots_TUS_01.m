
clear
close all
DATADIR   = "direction to the folder containing the data file"
load([DATADIR 'DATA_TUS_All_02.mat'])
FUS_DS=readtable([DATADIR 'TUS_DataCollection_BothMonkey.xlsx']); %reading excel sheet

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
        trialp=find(SessN==V_Sess(i) & TIB==1 & Trialdata.BlockNum==FUS_DS.Block_resume(i));
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
    
    learning_point=find(FF_Acc>.8 ,1,'first');
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
        BlockDATA.Asymp(i,1)=nanmean(Trialdata.Accuracy(Blkbrd(i+1)-12:Blkbrd(i+1)));
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






Npar=1;
ExplrLabel=sprintf(['Explorative \n Sampling (msec.) ']);
AsympLabel=sprintf(['Post-Learning \n Accuracy']);
ChoiceLabel=sprintf(['Choice \n Sampling(msec.)']);
RtLabel=sprintf(['Reaction \n time(msec.)']);

YLabel={'Trial-to-criterion',ExplrLabel,AsympLabel,RtLabel,ChoiceLabel}
Title={'Wotan','Igor','Both'};
DimLabel  = {'1D','2D','3D'};
TUSLabel  = FUSCnd_Label;
TokLabel = {'GL-framed','G-framed'};
close all
clear Blk
figure

ColDim=[.6 .6 .6; .4 .4 .4; .1 .1 .1];
ColTok=[.8 .2 .2; .2 .2 .8];
ColTUS = [.8 .5 .5; 1 .1 .1; .5 .8 .5; .1 1 .1];
pind=1:3;

%perm test
Blk.Par{1,1}=BlockDATA.LP;
Blk.Par{1,2}=BlockDATA.Eplr;
Blk.Par{1,3}=BlockDATA.Asymp;

clear RES;
clear Rand;

tic,
for ip=1:1000
    
    for Mnki=1:3
        
        for j=1:4
            
            idnG1 = find( BlockDATA.FUSCnd4~=j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2);
            idn1  =  find( BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2 );
            Lg = length(idnG1);
            L = length(idn1);
            Perm = randperm(Lg,L);
            
            for ipar=1:3
                RES.TUS{1,ipar}{1,Mnki}{j,1}=nanmean(Blk.Par{1,ipar}(idn1));
                Rand.TUS{1,ipar}{1,Mnki}{j,1}(ip,1)=nanmean(Blk.Par{1,ipar}(idnG1(Perm)));
                Rand.hTUS{1,ipar}{1,Mnki}{j,1}(ip,1)=prctile(Blk.Par{1,ipar}(idnG1(Perm)),95);
            end
            
            for iT=1:2
                
                idnG1 = find(   BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2);
                idn1  =  find( BlockDATA.TokCond==iT & BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2 );
                Lg = length(idnG1);
                L = length(idn1);
                Perm = randperm(Lg,L);
                
                for ipar=1:3
                    RES.TUSGL{1,ipar}{1,Mnki}{j,1}{iT,1}=nanmean(Blk.Par{1,ipar}(idn1));
                    Rand.TUSGL{1,ipar}{1,Mnki}{j,1}{iT,1}(ip,1)=nanmean(Blk.Par{1,ipar}(idnG1(Perm)));
                    Rand.hTUSGL{1,ipar}{1,Mnki}{j,1}{iT,1}(ip,1)=prctile(Blk.Par{1,ipar}(idnG1(Perm)),95);
                end
                
                for id=1:3
                    
                    idnG1 = find( BlockDATA.DimCond==id  & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2);
                    idn1  =  find( BlockDATA.DimCond==id & BlockDATA.TokCond==iT & BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2 );
                    Lg = length(idnG1);
                    L = length(idn1);
                    Perm = randperm(Lg,L);
                    
                    for ipar=1:3
                        RES.TUSGLD{1,ipar}{1,Mnki}{j,1}{iT,1}{id,1} = nanmean(Blk.Par{1,ipar}(idn1));
                        Rand.TUSGLD{1,ipar}{1,Mnki}{j,1}{iT,1}{id,1}(ip,1) = nanmean(Blk.Par{1,ipar}(idnG1(Perm)));
                        Rand.hTUSGLD{1,ipar}{1,Mnki}{j,1}{iT,1}{id,1}(ip,1) = prctile(Blk.Par{1,ipar}(idnG1(Perm)),95);
                    end
                    
                end
            end
        end
    end
end
toc
%%
clear Blk
for Mnki=1:3
    for iT=1:2
        for j=1:4
            for id=1:3
                
                idn=find( BlockDATA.DimCond==id & BlockDATA.TokCond==iT & BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==2 & BlockDATA.BlockNum<37);
                idn0=find( BlockDATA.DimCond==id & BlockDATA.TokCond==iT & BlockDATA.FUSCnd4==j & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc & BlockDATA.StimTCnd==1 & BlockDATA.BlockNum<37);
                
                Blk.Vec{1,1}{1,Mnki}{iT,j}{1,id}=BlockDATA.LP(idn);
                Blk.Vec{1,2}{1,Mnki}{iT,j}{1,id}=BlockDATA.Eplr(idn);
                Blk.Vec{1,3}{1,Mnki}{iT,j}{1,id}=BlockDATA.Asymp(idn);
                Blk.Vec{1,4}{1,Mnki}{iT,j}{1,id}=BlockDATA.RT(idn);
                Blk.Vec{1,5}{1,Mnki}{iT,j}{1,id}=BlockDATA.ObjG(idn);
                
                Blk.Vec0{1,1}{1,Mnki}{iT,j}{1,id}=BlockDATA.LP(idn0);
                Blk.Vec0{1,2}{1,Mnki}{iT,j}{1,id}=BlockDATA.Eplr(idn0);
                Blk.Vec0{1,3}{1,Mnki}{iT,j}{1,id}=BlockDATA.Asymp(idn0);
                Blk.Vec0{1,4}{1,Mnki}{iT,j}{1,id}=BlockDATA.RT(idn0);
                Blk.Vec0{1,5}{1,Mnki}{iT,j}{1,id}=BlockDATA.ObjG(idn0);
                
            end
        end
    end
end
xtik=[-5,0,5];
xtiks=[xtik+10,xtik+30,xtik+50,xtik+70];
XtickLabel={'_G_L','S-ACC','_G','_G_L','H-ACC','_G','_G_L','S-aSTR','_G','_G_L','H-aSTR','_G'};
Npar=1;


for v=1:Npar
    for Mnki=1:3
        
        k2=0;
        for j=1:4
            for iT=1:2
                k2=k2+1;
                RankP=sort(Rand.TUSGL{1,v}{1,Mnki}{j,1}{iT,1});
                M1=RES.TUSGL{1,v}{1,Mnki}{j,1}{iT,1};
                PvalsTok{v,Mnki}(1,k2)=(1-find([0;RankP;M1+1]>=M1,1,'first')/1000);
                PvalsTokL{v,Mnki}(1,k2)=(find([0;RankP;M1+1]>=M1,1,'first')/1000);
                for id=1:3
                    M1=RES.TUSGLD{1,v}{1,Mnki}{j,1}{iT,1}{id,1};
                    RankP=sort(Rand.TUSGLD{1,v}{1,Mnki}{j,1}{iT,1}{id,1});
                    Pvals{v,Mnki}(id,k2)=(1-find([0;RankP;M1+1]>=M1,1,'first')/1000);
                    PvalsL{v,Mnki}(id,k2)=(find([0;RankP;M1+1]>=M1,1,'first')/1000);
                end
                
                
            end
        end
    end
end
    
for v=1:Npar
    for Mnki=1:3
        PvalsTok{v,Mnki}(PvalsTok{v,Mnki}<0)=0;
        Pvals{v,Mnki}(Pvals{v,Mnki}<0)=0;
        PvalsTokL{v,Mnki}(PvalsTokL{v,Mnki}>1)=1;
        PvalsL{v,Mnki}(PvalsL{v,Mnki}>1)=1;
        h1=fdr_bh(PvalsTok{v,Mnki},.05,'pdep');
        h2=fdr_bh(Pvals{v,Mnki},.05,'pdep');
        h1l=fdr_bh(PvalsTokL{v,Mnki},.05,'pdep');
        h2l=fdr_bh(PvalsL{v,Mnki},.05,'pdep');
        P_ind = sub2ind([ 3 Npar],Mnki,v);
        subplot(Npar,3, P_ind)
        hold on
        k=0
        YG=[];
        YG=cell(1,3);
        for id=1:3
            for j=1:4
                for iT=1:2
                    YG{1,id}=[YG{1,id},Blk.Vec{1,v}{1,Mnki}{iT,j}{1,id}',Blk.Vec0{1,v}{1,Mnki}{iT,j}{1,id}'];
                end
            end
        end
        for j=1:4
            for iT=1:2
                k=k+1;
                
                for id=1:3
                    
                    Y1=Blk.Vec{1,v}{1,Mnki}{iT,j}{1,id};
                    Y0=Blk.Vec0{1,v}{1,Mnki}{iT,j}{1,id};
                    Y1=(Y1-nanmean(YG{1,id}))./nanstd(YG{1,id});
                    
                    Mn=nanmean(Y1);
                    SE=nanstd(Y1)./sqrt(sum(~isnan(Y1)));
                    imSp=linspace(-SE,SE,10);
                    %         plot([id id],[Mn-SE Mn+SE],'color',ColTUS(j,:),'linewidth',2)
                    %         scatter(id,Mn,18,'markerfacecolor',ColTUS(j,:),'markeredgecolor',ColTUS(j,:))
                    Yim{k,1}(:,id)=Mn+imSp;
                end
                
                
            end
        end
        
        im=cell2mat(Yim);
        imagesc(1:80,1:3,im','AlphaData',0.8)
        colormap jet
        c=colorbar;
        for ix=20.5:20:80.5
            plot([ix ix],[0.5 3.5],'color','k','linewidth',3)
            plot([ix-10 ix-10],[0.5 3.5],'color','k','linewidth',.5,'linestyle','--')
        end
        for id=1:3
            plot([.5 80.5],[id+.5 id+.5],'color','k','linewidth',1)
            
        end
        set(gca,'xlim',[.5 80.5] ,'ylim',[0.5 3.5],'tickdir','out','xtick',xtiks,'xticklabel',XtickLabel,'ytick',1:3,'yticklabel',DimLabel)
        if Mnki==1
            c.Label.String = ['Diff. Norm.' YLabel{v}];
            
           
        end
        if v==1
            title(['Monkey:  ', Title{Mnki}])
        end
        k1=0;
        for j=1:4
            
            for iT=1:2
                k1=k1+1;

                if (h1(1,k1) || h1l(1,k1))
                    plot([10*k1-9 10*k1],[.5 .5] ,'linewidth',5,'color','w')
                    plot([10*k1-9 10*k1] ,[3.5 3.5],'linewidth',5,'color','w')
                    plot([10*k1-9 10*k1-9] ,[.5 3.5],'linewidth',5,'color','w')
                    plot([10*k1 10*k1] ,[.5 3.5],'linewidth',5,'color','w')
                    
                end
                    
                for id=1:3

                    if (h2(id,k1) || h2l(id,k1))
                        scatter (10*k1-5,id,65,'*','w')
                    end
                end
            end
        end
        
    end
end


set(gcf, 'Position',[    110         806        1152         195])

%% Session plot for load


Title={'Wotan','Igor','Both'};
close all
figure
for Mnki=1:3
    
    Sess=[];
    Sess.SeshNum=V_Sess;
    
    for i=1:length(V_Sess)
        for j=1:4
            for id=1:3
                idn=find(BlockDATA.SessionNum==V_Sess(i));
                idnBs=find(BlockDATA.SessionNum==V_Sess(i) & BlockDATA.StimTCnd==1 & BlockDATA.TokCond~=2 & BlockDATA.FUSCnd4==j & BlockDATA.DimCond==id & BlockDATA.BlockNum<37 & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc);
                idnTUS=find(BlockDATA.SessionNum==V_Sess(i) & BlockDATA.StimTCnd==2 & BlockDATA.TokCond~=2 & BlockDATA.FUSCnd4==j & BlockDATA.DimCond==id & BlockDATA.BlockNum<37 & BlockDATA.Mnk~=pind(Mnki) & BlockDATA.Inc);
                Sess.LP_Bs(id,i,j)=nanmean(BlockDATA.LP(idnBs));
                Sess.LP_TUS(id,i,j)=nanmean(BlockDATA.LP(idnTUS));
                Sess.Eplr_Bs(id,i,j)=nanmean(BlockDATA.Eplr(idnBs));
                Sess.Eplr_TUS(id,i,j)=nanmean(BlockDATA.Eplr(idnTUS));
                Sess.Asymp_Bs(id,i,j)=nanmean(BlockDATA.Asymp(idnBs));
                Sess.Asymp_TUS(id,i,j)=nanmean(BlockDATA.Asymp(idnTUS));
            end
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
    
    
    for i=1
        P_ind = sub2ind([ 3 1],Mnki,i);
        subplot(1,3,P_ind)
        hold on
        
        
        for j=1:4
            for id=1:3
                Y1=squeeze(VecSess{i,2}(id,:,j));
                Y0=squeeze(VecSess{i,1}(id,:,j));
                
                Mn=nanmean(Y1);
                SE=nanstd(Y1)./sqrt(sum(~isnan(Y1)));
                plot([id id],[Mn-SE Mn+SE],'color',ColTUS(j,:),'linewidth',2)
                scatter(id,Mn,18,'markerfacecolor',ColTUS(j,:),'markeredgecolor',ColTUS(j,:))
                Yline(id)=Mn;
                Xline(id)=id;
            end
            plot(Xline,Yline,'color',ColTUS(j,:),'linewidth',1)
        end
        
        set(gca,'xlim',[0 4],'tickdir','out','xtick',1:3,'xticklabel',DimLabel)
        if Mnki==1
            ylabel(YLabel{i})
        end
        if i==1
            title(['Monkey:  ', Title{Mnki}])
        end
    end
    % sum(~isnan(VecSess{i,j}(:,2)))
end
set(gcf, 'Position',[   584   844   678   157])