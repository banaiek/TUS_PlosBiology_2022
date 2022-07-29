%ErrorBar_Plot_01('baseline','ptype','scatter','facealpha','xticklabel','LineWidth','LineStyle','color','edgecolor','Nonpar')

function ErrorBar_Plot_02(Y,varargin)
%defaults:
%_____________________
X=[];
col=[0.5 0.5 0.5];
faceC=[0.5 0.5 0.5];
alphap=0.5;
edgeC=[0.1 0.1 0.1];
barplt=1;
scpts=0;
plotlegend=0;
Stype=0;
Y0=[];
Xlabel = [];
Etype='SE';
ls={'-','--'};
lw=[1,1];
sz=10;
mkr='o';



%get additional input parameters (varargin)
if~isempty(find(strcmp(varargin,'xticklabel')))
    Xlabel= varargin{find(strcmp(varargin,'xticklabel'))+1};
end
if~isempty(find(strcmp(varargin,'MarkerSize')))
    sz= varargin{find(strcmp(varargin,'MarkerSize'))+1};
end
if~isempty(find(strcmp(varargin,'ErrorType')))
    Etype= varargin{find(strcmp(varargin,'ErrorType'))+1};
end
if~isempty(find(strcmp(varargin,'linewidth')))
    lw= varargin{find(strcmp(varargin,'linewidth'))+1};
end
if~isempty(find(strcmp(varargin,'linestyle')))
    ls= varargin{find(strcmp(varargin,'linestyle'))+1};
end

if ~isempty(find(strcmp(varargin,'color')))
    col = varargin{find(strcmp(varargin,'color'))+1};
end

if ~isempty(find(strcmp(varargin,'edgecolor')))
    edgeC = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if ~isempty(find(strcmp(varargin,'facealpha')))
    alphap = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if ~isempty(find(strcmp(varargin,'ptype')))
    barplt = varargin{find(strcmp(varargin,'ptype'))+1};
end
if ~isempty(find(strcmp(varargin,'scatter')))
    scpts = varargin{find(strcmp(varargin,'scatter'))+1};
end
if ~isempty(find(strcmp(varargin,'plotlegend')))
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end

if ~isempty(find(strcmp(varargin,'Nonpar')))
    Stype = varargin{find(strcmp(varargin,'Nonpar'))+1};
end

if ~isempty(find(strcmp(varargin,'baseline')))
    Y0 = varargin{find(strcmp(varargin,'baseline'))+1};
end
% 
% Ymin=min([min(Y,[],'all'),min(Y0,[],'all')]);
% Ymax=max([max(Y,[],'all'),max(Y0,[],'all')]);

if iscell(Y)==0
    Y = num2cell(Y,1);
end

L=size(Y,2);


if (~isempty(Y0) && ~iscell(Y0))
    Y0 = num2cell(Y0,1);
end

if isempty('Xlabel')
    Xlabel = [1:L];
end

x = 1:L;
Mn = cellfun(@(x) nanmean(x,1),Y);
Md = cellfun(@(x) nanmedian(x,1),Y);
Ls = cellfun(@(x) sum(~isnan(x)),Y);
SD = cellfun(@(x) nanstd(x,[],1),Y);
SE = SD./sqrt(Ls);


if ~Stype
    y=Mn;
else
    y=Md;
end

if strcmp(Etype,'SD')
    yd=y-SD;
    yu=y+SD;
else
    yd=y-SE;
    yu=y+SE;
end


if ~isempty(Y0)
    x0 = [1:L]-1/4;
    Mn0 = cellfun(@(x) nanmean(x,1),Y0);
    Md0 = cellfun(@(x) nanmedian(x,1),Y0);
    Ls0 = cellfun(@(x) sum(~isnan(x)),Y0);
    SD0 = cellfun(@(x) nanstd(x,[],1),Y0);
    SE0 = SD0./sqrt(Ls0);
    
    if ~Stype
        y0=Mn0;
    else
        y0=Md0;
    end
y(isnan(y0))=[];    
    if strcmp(Etype,'SD')
        yd0=y0-SD0;
        yu0=y0+SD0;
    else
        yd0=y0-SE0;
        yu0=y0+SE0;
    end
end
% plot Error Bars
if ~barplt
    
%    Ymin=min([yd0-SE0,yd-SE]); 
%    Ymax=max([yu0+SE0,yu+SE0]); 
    
    for i=1:L
        Y{1,i}(isnan(Y{1,i}))=[];
        if size(col,1)==1
            if scpts
                scatter(i*ones(1,Ls(i)),Y{1,i},sz/5,col,'filled',mkr)
            end
            plot([x(i) ;x(i) ],[yd(i);yu(i)],'color',col,'linewidth',lw(1),'linestyle',ls{1})
            scatter(x(i),y(i),sz,col,'filled',mkr)
        else
            if scpts
                scatter(i*ones(1,Ls(i)),Y{1,i},sz/5,col(i,:),'filled',mkr)
            end
            plot([x(i) ;x(i) ],[yd(i);yu(i)],'color',col(i,:),'linewidth',lw(1),'linestyle',ls{1})
            scatter(x(i),y(i),sz,col(i,:),'filled',mkr)
        end
    end
    
    if ~isempty(Y0)
        
        
        for i=1:L
            Y0{1,i}(isnan(Y0{1,i}))=[];
            if size(col,1)==1
                if scpts
                    scatter(x0(i)*ones(1,Ls(i)),Y0{1,i},sz/5,col,'filled',mkr)
                end
                plot([x0(i) ;x0(i) ],[yd0(i);yu0(i)],'color',col,'linewidth',lw(2),'linestyle',ls{2})
                scatter(x0(i),y0(i),sz,col,'filled',mkr)
            else
                if scpts
                    scatter(i*ones(1,Ls(i)),Y0{1,i},sz/5,col(i,:),'filled',mkr)
                end
                plot([x0(i) ;x0(i) ],[yd0(i);yu0(i)],'color',col(i,:),'linewidth',lw(1),'linestyle',ls{2})
                scatter(x0(i),y0(i),sz,col(i,:),mkr)
            end
        end
        set(gca,'xlim',[0 L+1],'xtick',[1:L]-(1/8),'xticklabel',Xlabel,'tickdir','out')
        
    else
        set(gca,'xlim',[0 L+1],'xtick',[1:L],'xticklabel',Xlabel,'tickdir','out')
    end
    
else
    Cons=0;
    if ~isempty(Y0)
        Cons=.2;
    end
    b1=bar(x+Cons,y,.3,'FaceColor','flat','FaceAlpha',alphap,'EdgeColor',edgeC);
    hold on
    b1.CData=col;
    if scpts
        for i=1:L
            Y{1,i}(isnan(Y{1,i}))=[];
            scatter((x(i)+Cons)*ones(1,Ls(i)),Y{1,i},sz,col(i,:),'filled',mkr)
            plot([(x(i)+Cons) ;(x(i)+Cons)],[y(i);yu(i)],'color',[.1 .1 .1],'linewidth',2)
        end
    end
    if ~isempty(Y0)
        b0=bar(x-Cons,y0,.3,'FaceColor','flat','FaceAlpha',alphap,'EdgeColor',edgeC);
        b0.CData=[.6 .6 .6];
        if scpts
            for i=1:L
                Y0{1,i}(isnan(Y0{1,i}))=[];
                scatter((x(i)-Cons)*ones(1,length(Y0{1,i})),Y0{1,i},sz,[.5 .5 .5],mkr)
                plot([x(i)-Cons ;x(i)-Cons ],[y0(i);yu0(i)],'color',[.1 .1 .1],'linewidth',2)
            end
        end
        
    end
    set(gca,'xlim',[0 L+1],'xtick',[1:L],'xticklabel',Xlabel,'tickdir','out')
end


end



