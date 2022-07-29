%%
clear
close all

DATADIR=" dir to the dataset folder";
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

VECFit=nan(NumT,1);



%% Getting Block data

lb=[-1,-1,1,.7,0];
hb=[1,1,50,1,3];
a=[0,0,15,.8,1];

col={[.8,.1,.1], [.1,.1,.8],[.1,.8,.1]};
Col4Cnd=[.8 .5 .5; 1 .1 .1 ;.5 .8 .5; .1 1 .1];
clear BlockDATA    
for i=1:NumBlk
    T_occ(Blkbrd(i)+1:Blkbrd(i+1))=mod(cumsum([0; OC(Blkbrd(i)+1:Blkbrd(i+1)-1)]),5);
 
    Acc_OC=Trialdata.Accuracy(Blkbrd(i)+1:Blkbrd(i)+35);
    CONVAC=conv(Acc_OC,ones(1,12))/12;
    FF_Acc=CONVAC(12:end-11);

        %
    Y=Trialdata.Accuracy(Blkbrd(i)+1:Blkbrd(i+1));
    [logitCoef,dev] = glmfit(1:length(Y),Y,'binomial','logit');
    logitFit = glmval(logitCoef,1:length(Y),'logit');
    LCss(i,1:length(Y))=logitFit;
%     plot(1:length(Y),Y,'bs', 1:length(Y),logitFit,'r-');
%     xlabel('Weight'); ylabel('Proportion');
    BlockDATA.CFC(i,1:2)=(logitCoef');
    VECFit(Blkbrd(i)+1:Blkbrd(i+1))=logitFit;
     [cf,G]=L5P([1:length(logitFit)]',logitFit,a,lb,hb);
     
     BlockDATA.EL(i,1)=cf.C;
     BlockDATA.RRate(i,1)=cf.A;
     BlockDATA.RAsymp(i,1)=cf.B;


    learning_point=find(logitFit>.80 ,1,'first');
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
        lp=round(predict(mDl,.85));
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



%% Plot session data
Col_A(1,:,:) = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];
Col_A(2,:,:) = [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];
Col_A(3,:,:) = [.8 .5 .5; 1 .1 .1; .5 .8 .5; .1 1 .1];
ylims=[5 20; .88  .94; .15 .25; .4 .55]; 
Marker={'^','s','o'};
LW=[2,2,3];
Vec{1}=BlockDATA.LP;
Vec{2}=BlockDATA.Asymp;
Vec{3}=BlockDATA.Eplr;
Vec{4}=BlockDATA.RT;


%% Plot learning curve and exploration rate
ColDim=[.6 .6 .6; .4 .4 .4; .1 .1 .1];
ColTok=[.8 .2 .2; .2 .2 .8];
ColTUS = [.8 .5 .5; 1 .1 .1; .5 .8 .5; .1 1 .1];
Ylim=[.2 1; 0 .3; 0 .15; 0 1]; 
Title={'Wotan:  ','Igor:  ','Both:  '};
YLabel={'GLMfit','Exploration (msec.)','ChosenObje_Pre_sampling (msec.)', 'Proportion Correct'} ;
Vec{1}=VECFit;
Vec{2}=Event.Explore;
Vec{3}=Event.ChosenObjGaze;
Vec{4}=Trialdata.Accuracy;

close all
for iv=1:4
Nt=60;
X=1:30;  %Trial Numbers
y=Vec{iv};
smw=([0:1:59]/40).^2;
flw=round(exp(smw));
flw=fliplr(max(flw)-flw)+1;
%learning curve for different loads
%learning curve for different Token Conditions


figure
pind=1:3;

for i=1:3
    subplot(5,3,i)
    hold on
    for j=1:3
        for it=X
            
            inds=find(iCndDim'==j & MnkID~=pind(i) & TIB==it);
            l2 = length(inds);
            if it==1
                l1=length(inds);
                Yout=nan(l1,Nt);
                Ys=nan(l1,Nt);
            end
            Yout(1:l2,it)=y(inds);
            
            
        end
        for is=X
            Ys(:,is)=nanmean(Yout(:,is:is+flw(is)),2);
        end
        Mn=nanmean(Ys);
        SE=nanstd(Ys)./sqrt(l1);
        shadedErrorBar(X,Mn(X),SE(X),{'color',ColDim(j,:)},1)
    end
    set(gca,'tickdir','out','xlim',[X(1) X(end)],'ylim',Ylim(iv,:))
    plot(get(gca,'xlim'),[.8 .8],'linestyle','--','linewidth',1,'color','k')
    title ([Title{i} , 'Att. load'])
        if i==1
        ylabel(YLabel{iv})
        end
end




for i=1:3
    subplot(5,3,i+3)
    hold on
    for j=1:2
        for it=X
            
            inds=find(iCndTok'==j & MnkID~=pind(i) & TIB==it);
            l2 = length(inds);
            if it==1
                l1=length(inds);
                Yout=nan(l1,Nt);
                Ys=nan(l1,Nt);
            end
            Yout(1:l2,it)=y(inds);
            
            
        end
        for is=X
            Ys(:,is)=nanmean(Yout(:,is:is+flw(is)),2);
        end
        Mn=nanmean(Ys);
        SE=nanstd(Ys)./sqrt(l1);
        shadedErrorBar(X,Mn(X),SE(X),{'color',ColTok(j,:)},1)
    end
    set(gca,'tickdir','out','xlim',[X(1) X(end)],'ylim',Ylim(iv,:))
    plot(get(gca,'xlim'),[.8 .8],'linestyle','--','linewidth',1,'color','k')
    title ([Title{i} , 'Token Cnd'])
        if i==1
        ylabel(YLabel{iv})
        end
end



Ys=[];
for i=1:3
    subplot(5,3,i+6)
    hold on
    for j=1:4
        for it=X
            
            inds=find(FUS.Cnd4==j & iCndDim'>0  & FUS.PP ==2 & MnkID~=pind(i) & TIB==it & FUS.Inc);
            l2 = length(inds);
            if it==1
                l1=length(inds);
                Yout=nan(l1,Nt);
                Ys=nan(l1,Nt);
            end
            Yout(1:l2,it)=y(inds);
            
            
        end
        for is=X
            Ys(:,is)=nanmean(Yout(:,is:is+flw(is)),2);
        end
        Mn=nanmean(Ys);
        SE=nanstd(Ys)./sqrt(l1);
        shadedErrorBar(X,Mn(X),SE(X),{'color',ColTUS(j,:)},1)
    end
    set(gca,'tickdir','out','xlim',[X(1) X(end)],'ylim',Ylim(iv,:))
    plot(get(gca,'xlim'),[.8 .8],'linestyle','--','linewidth',1,'color','k')
    title ([Title{i} , 'TUS Cnd. All'])
        if i==1
        ylabel(YLabel{iv})
        end
end

Ys=[];
for i=1:3
    subplot(5,3,i+9)
    hold on
    for j=1:4
        for it=X
            
            inds=find(FUS.Cnd4==j & iCndDim'>0 & iCndTok'==1 & FUS.PP ==2 & MnkID~=pind(i) & TIB==it & FUS.Inc);
            l2 = length(inds);
            if it==1
                l1=length(inds);
                Yout=nan(l1,Nt);
                Ys=nan(l1,Nt);
            end
            Yout(1:l2,it)=y(inds);
            
            
        end
        for is=X
            Ys(:,is)=nanmean(Yout(:,is:is+flw(is)),2);
        end
        Mn=nanmean(Ys);
        SE=nanstd(Ys)./sqrt(l1);
        shadedErrorBar(X,Mn(X),SE(X),{'color',ColTUS(j,:)},1)
    end
    set(gca,'tickdir','out','xlim',[X(1) X(end)],'ylim',Ylim(iv,:))
    plot(get(gca,'xlim'),[.8 .8],'linestyle','--','linewidth',1,'color','k')
        title ([Title{i} , 'TUS Cnd. GL'])
        if i==1
        ylabel(YLabel{iv})
        end
end

Ys=[];
for i=1:3
    subplot(5,3,i+12)
    hold on
    for j=1:4
        for it=X
            
            inds=find(FUS.Cnd4==j & iCndDim'>0 & iCndTok'==2 & FUS.PP ==2 & MnkID~=pind(i) & TIB==it & FUS.Inc);
            l2 = length(inds);
            if it==1
                l1=length(inds);
                Yout=nan(l1,Nt);
                Ys=nan(l1,Nt);
            end
            Yout(1:l2,it)=y(inds);
            
            
        end
        for is=X
            Ys(:,is)=nanmean(Yout(:,is:is+flw(is)),2);
        end
        Mn=nanmean(Ys);
        SE=nanstd(Ys)./sqrt(l1);
        shadedErrorBar(X,Mn(X),SE(X),{'color',ColTUS(j,:)},1)
    end
    set(gca,'tickdir','out','xlim',[X(1) X(end)],'ylim',Ylim(iv,:))
    plot(get(gca,'xlim'),[.8 .8],'linestyle','--','linewidth',1,'color','k')
        title ([Title{i} , 'TUS Cnd. G'])
        if i==1
        ylabel(YLabel{iv})
        end
        xlabel 'Trial in block'
end
set(gcf,'position',[520    57   652   892])
end