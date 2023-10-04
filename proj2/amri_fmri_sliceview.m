%%
% amri_fmri_sliceview: view mri/fmri dataset as slices
%
% Usage
%   [uimg,oimg]=amri_fmri_sliceview(ulay [,olay,'keyword','value'...]);
%   [uimg,oimg]=amri_fmri_sliceview(ulay,olay,'ui',1);
%
% Inputs
%   ulay: 3D or 4D array
%   olay: 3D or 4D array
%
% Keywords
%       ui: user interface [0|1]
%     mask: 3D mask applied to olay
%      new: create a new figure [0|1]
%     grid: [nrow ncol]
%    thres: [lower upper] for olay
%    oclim: [ocmin ocmax] for olay
%    uclim: [ucmin ucmax] for ulay
%    ocmap: overlay colormap {default: 'y2c'}
%    ocgap: using gap in colormap for thresholded overlay {1|0,default 0}
%    ucmap: underlay colormap {default: 'gray'}
%    alpha: control transparency for overlay
%     zoom: zoom ratio. zoomin (>1) or zoomout (<1) {default: 1}
%
% Version
%   1.17

%% History
% 1.00 - 05/27/2010 - ZMLIU - create the original file
% 1.01 - 06/10/2010 - ZMLIU - add 'grid' keyword to customize display matrix
% 1.02 - 06/25/2010 - ZMLIU - fix a minor "bug" regarding ucmin, ucmax
% 1.03 - 06/29/2010 - ZMLIU - assemble all slices into one image and display it as a whole
% 1.04 - 06/29/2010 - ZMLIU - display overlay on top of underlay
% 1.05 - 07/01/2010 - ZMLIU - make overlay transparent using 'alpha' keyword
% 1.06 - 07/29/2010 - ZMLIU - default colormap 'yellow-red-blue-cray'
% 1.07 - 07/30/2010 - ZMLIU - add colorbar if specified
% 1.08 - 08/13/2010 - ZMLIU - convert int16 to double for ucmax,ucmin etc
% 1.09 - 12/20/2010 - ZMLIU - add 'ocgap' keyword
% 1.10 - 01/06/2011 - ZMLIU - fix a minor bug 
% 1.11 - 09/12/2011 - ZMLIU - allow for a mask keyword
% 1.12 - 04/06/2012 - ZMLIU - add a new colormap "spectral" with the use of
%                           - spectral function developed by K. Worsley (in
%                           - memory of a great fMRI analysis expert)
% 1.13 - 05/24/2012 - ZMLIU - allow oclim to contain four numbers
%                           - [neg_min neg_max pos_min pos_max]
%                           - ocgap is not used anymore 
% 1.14 - 12/03/2013 - ZMLIU - remove ocgap keyword
% 1.15 - 06/09/2015 - ZMLIU - add basic GUI
% 1.16 - 07/09/2015 - ZMLIU - plot multiple time series on click
% 1.17 - 08/03/2021 - ZMLIU - return h_image from show_underlay


%%

function amri_fmri_sliceview(ulay,varargin)

if nargin<1
    eval('help amri_fmri_sliceview');
    return
end

if ndims(ulay)<2
    error('amri_fmri_sliceview(): ulay dimensions must be larger than 1');
elseif ismatrix(ulay)
    ulay = reshape(ulay,size(ulay,1),size(ulay,2),1);
elseif ndims(ulay)==4
    ulay_data = ulay;
    ulay = sqrt(sum(ulay.^2,4));
end
    
%% default settings
uval=[];
oval=[];
uimg=[];
oimg=[];

ncol=nan;          % # of slices in each row
nrow=nan;          % # of slices in each column
zoomfac=nan;       % zoom in (>1) or out (<1)
olower=0;          % lower threshold for overlay
oupper=0;          % upper threshold for overlay
ocmap='y2c';       % coloarmap for overlay (yellow-red-blue-cray)
% ocgap=0;           % use gap in colorbar or not
oalpha=1;          % transparency parameter
ucmap='gray';      % colormap for underlay
ucmapni=256;       % number of indices used in the underlay colormap
ocmapni=256;       % number of indices used in the overlay colomap
% colorbar_width=15; % colorbar width
% textbox_height=30; % textbox height
ui_width=50;       % ui width
ui_height=20;      % ui height
flag_newfig=0;     % new figure (1) or not (0)
flag_ui=1;         % display user interface

h_plot=nan;        % handle to plot figure
h_tens=nan;        % handle to tensor figure

% flag_ucolorbar=0;  % display underlay colorbar (1) or not (0)
% flag_ocolorbar=0;  % display overlay colorbar (1) or not (0)
% flag_utextbox=0;   % show text message (1) or not (0) for underlay
% flag_otextbox=0;   % show text message (1) or not (0) for overlay

% find the intensity range for ulay [ucmin ucmax]
orig_size=size(ulay);
ulay=ulay(:);
ucmin=double(min(ulay(ulay~=0)));
ucmax=double(max(ulay(ulay~=0)));
if ucmin==ucmax
    ucmin=min([ucmin 0]);
    ucmax=max([ucmax 0]);
end
ulay=reshape(ulay,orig_size);

%% keyword-value pairs
if rem(nargin,2)==1
     kwstart = 1;
%      flag_ocolorbar=0;
%      flag_otextbox=0;
else
     olay=varargin{1};
     if ndims(olay)==4
         olay_data=olay;
         olay=olay(:,:,:,1);
     end
     kwstart = 2;
end

for i=kwstart:2:size(varargin,2)
    Keyword = varargin{i};
    Value   = varargin{i+1};
    if ~ischar(Keyword) 
        printf('amri_fmri_sliceview(): keywords must be strings')
        return
    end
    if strcmpi(Keyword,'newfig') || ...
       strcmpi(Keyword,'new')
        if isnumeric(Value)
            flag_newfig=(Value(1)~=0);
        end
    elseif strcmpi(Keyword,'mask')
        mask=Value;
    elseif strcmpi(Keyword,'grid') && isnumeric(Value)
        if length(Value)==1
            nrow=Value;
        else
            nrow=Value(1);
            ncol=Value(2);
        end
    elseif (strcmpi(Keyword,'zoomfac') || ...
            strcmpi(Keyword,'zoomfactor') || ...
            strcmpi(Keyword,'zoom')) && ...
            isnumeric(Value)
        zoomfac=Value(1);
        if zoomfac<0
            zoomfac=1;
        end
    elseif strcmpi(Keyword,'ocmap')
        ocmap=Value;
    elseif strcmpi(Keyword,'ucmap')
        ucmap=Value;
    elseif strcmpi(Keyword,'alpha') || strcmpi(Keyword,'oalpha')
        oalpha=Value;
    elseif strcmpi(Keyword,'uclim')
        if isscalar(Value)
            ucmax = abs(Value);
            ucmin = 0-abs(Value);
        elseif isvector(Value)
            if length(Value)==2
                ucmin = Value(1);
                ucmax = Value(2);
            end
        end
        ucmin=double(ucmin);
        ucmax=double(ucmax);
    elseif strcmpi(Keyword,'thres') || strcmpi(Keyword,'threshold')
        if isnumeric(Value)
            if isscalar(Value)
                olower=-Value;
                oupper=Value;
            else
                olower=Value(1);
                oupper=Value(2);
            end
        end
    elseif strcmpi(Keyword,'oclim')
        if isscalar(Value)
            oclim=[-abs(Value) abs(Value)]';
        else
            oclim=sort(Value(:),1,'ascend');
        end
    elseif strcmpi(Keyword,'gui') || strcmpi(Keyword,'ui')
        if Value(1)>0
            flag_ui=1;
        end
% 2012-05-24        
%         if isscalar(Value)
%             ocmax = abs(Value);
%             ocmin = 0-abs(Value);
%         elseif isvector(Value)
%             if length(Value)==2
%                 ocmin = Value(1);
%                 ocmax = Value(2);
%             end
%         end     
%         ocmin=double(ocmin);
%         ocmax=double(ocmax);
%     elseif strcmpi(Keyword,'colorbar')
%         flag_ocolorbar=Value;
%         flag_ucolorbar=Value;
%     elseif strcmpi(Keyword,'ocgap')
%         ocgap=Value(1);
    else
        fprintf(['amri_fmri_sliceview(): ' Keyword ' is an unknown keyword']);
    end
end

%%
if ischar(ucmap)
    ucmap=eval([ucmap '(' num2str(ucmapni) ')']);
elseif isnumeric(ucmap)
    ucmapni=size(ucmap,1);
    ucmap=ucmap/max(ucmap(:));
end

[nx,ny,nz] = size(ulay);

if isempty(findobj('Type','figure'))
    flag_newfig=1;
end

if flag_newfig
    figure;
else
    clf;
end

h_figure=gcf;

junk=get(h_figure,'Position');
fig_left   = junk(1);
fig_bottom = junk(2);
fig_width  = junk(3);
fig_height = junk(4);
set(gcf,'Color','k');   % set bkcolor to 'black'

% compute [ncol nrow] that best fits the current figure size
if isnan(ncol) && isnan(nrow)
    nrow = round(sqrt((nx/ny)/(fig_width/fig_height)*nz));
    ncol = ceil(nz/nrow);
elseif isnan(ncol) && ~isnan(nrow)
    ncol = ceil(nz/nrow);
end

% calculate zoom factor, if not specified, to best fit the figure size
if isnan(zoomfac)
    a=(fig_width-flag_ui*ui_width)/(nx*ncol);
    b=(fig_height-flag_ui*ui_height)/(ny*nrow);
    zoomfac=min(a,b);
end    

img_width =nx*ncol*zoomfac;
img_height=ny*nrow*zoomfac;
% reset figure size and position
fig_width = img_width+flag_ui*ui_width;
fig_height= img_height+flag_ui*ui_height;
scrsize=get(0,'ScreenSize');
if flag_newfig==1
    set(gcf,'Position',[0 scrsize(4)-fig_height fig_width fig_height]);
else
    set(gcf,'Position',[fig_left fig_bottom fig_width fig_height]);
end

% 08/03/2021 - return h_image from show_underlay
h_image = show_underlay;    

show_overlay;


% image click callback function
set(h_image,'ButtonDownFcn',@ImageClick_Callback);

%% colorbars
% 
% if flag_ucolorbar==1
%     subplot('Position',...
%         [(img_width+1)/fig_width ...
%          ((flag_utextbox+flag_otextbox)*textbox_height+1)/fig_height ...
%          colorbar_width/fig_width ...
%          img_height/fig_height]);
%     ucb=repmat((ucmapni:-1:1)',1,colorbar_width);
%     ucb=ind2rgb(ucb,ucmap);
%     image(ucb);
%     %axis equal;
%     axis off;
% end
% 
% % colorbar for overlay
% if flag_ocolorbar==1 && exist('olay','var')
%     subplot('Position', ...
%          [(flag_ucolorbar*colorbar_width+img_width+1)/fig_width ...
%           ((flag_utextbox+flag_otextbox)*textbox_height+1)/fig_height ...
%           colorbar_width/fig_width ...
%           img_height/fig_height]);
%     ocb=repmat((ocmapni:-1:1)',1,colorbar_width);
%     ocb=ind2rgb(ocb,ocmap);
%     image(ocb);
%     %axis equal;
%     axis off;
% end
% 

%% GUI

if flag_ui>0
    junk=get(gcf,'Position');
    fig_width=junk(3);
    fig_height=junk(4);
    hinfo=uicontrol('Style','pushbutton','String','Info',...
        'units','normalized',...
        'Position',[(fig_width-ui_width)/fig_width ...
                    (fig_height-ui_height*1)/fig_height ...
                    ui_width/fig_width ui_height/fig_height],...
        'Callback',@infobutton_Callback);
    htext=uicontrol('Style','text','String','','units','normalized',...
        'Position',[0 0 1 ui_height/fig_height]);
end


    %% nested function: show data information
    function infobutton_Callback(source,eventdata) %#ok<INUSD>
        set(hinfo,'Enable','off');
        str=['underlay:\n'  ...
             ' max=' num2str(max(ulay(ulay~=0))) ';\t' ...
             ' min=' num2str(min(ulay(ulay~=0))) '\n'];
        if exist('olay','var')
            str=[str 'overlay:\n' ...
                ' max=' num2str(max(olay(olay~=0))) ';\t' ...
                ' min=' num2str(min(olay(olay~=0))) '\n'];
        end
        msgbox(sprintf(str));
        set(hinfo,'Enable','on');
    end

    %% nested function: respond to image click event
    function ImageClick_Callback(objectHandle,eventData) %#ok<INUSL>
        [nx,ny,nz]=size(ulay);
        
        % get the key option
        key_option=get(gcf,'CurrentCharacter');
        
        % get (x,y) of the click point
        x=eventData.IntersectionPoint(1);
        y=eventData.IntersectionPoint(2);
        % calculate [I,J,K] in the 3-D image space
        I=round(mod(x,nx));
        J=ny-round(mod(y,ny));
        K=floor(y/ny)*ncol+ceil(x/nx); 

        % zmliu 7/8/15
        global IJK_Click
        IJK_Click=[I J K];
        
        % display [I,J,K] and image intensity in the status bar
        str=['IJK=[' int2str(I) ',' int2str(J) ',' int2str(K) ']; '];
        str=[str 'ulay=' num2str(ulay(I,J,K)) '; '];
        if exist('olay','var')
            str=[str 'olay=' num2str(olay(I,J,K))]; 
        end
        set(htext,'String',str);
        
        % plot data or information
        if eventData.Button==1
            % press left button
            if exist('ulay_data','var')
                prep_plot_figure(key_option);
                plot_data(ulay_data,[I J K],key_option);
            end
        elseif eventData.Button==3
            % press right button
            if exist('olay_data','var')
                prep_plot_figure(key_option);
                plot_data(olay_data,[I J K],key_option);
            end
        elseif eventData.Button==2
            % press middle button
            if exist('olay_data','var') || exist('ulay_data','var')
                prep_plot_figure(key_option);
                if exist('ulay_data','var')
                    plot_data(ulay_data,[I J K],key_option);
                end
                hold on;
                if exist('olay_data','var')
                    plot_data(olay_data,[I J K],key_option);
                end
            end
        end
        
%         % plot correlation tensor
%         if eventData.Button==1
%             % press left button
%             if exist('ulay_data','var')
%                 prep_tensor_figure(key_option);
%                 plot_tensor(ulay_data,[I J K], key_option);
%             end
%         elseif eventData.Button==3
%             % press right button
%             if exist('olay_data','var')
%                 prep_tensor_figure(key_option);
%                 plot_tensor(olay_data,[I J K], key_option);
%             end
%         elseif eventData.Button==2
%             % press middle button
%             if exist('olay_data','var') || exist('ulay_data','var')
%                 prep_tensor_figure(key_option);
%                 if exist('ulay_data','var')
%                     plot_tensor(ulay_data,[I J K],key_option);
%                 end
%                 hold on;
%                 if exist('olay_data','var')
%                     plot_tensor(olay_data,[I J K],key_option);
%                 end
%             end
%         end
        
        % set the main figure as the active figure
        figure(h_figure);
    end

    %% nested function: prepare a time series plot figure
    function prep_plot_figure(key)
        if ishghandle(h_plot)
            set(0,'CurrentFigure',h_plot);
        else
            h_plot=figure;
        end
        set(gcf,'color',get(gcf,'defaultFigureColor'));
        clf;
    end
    
    %% nested function: plot time series data
    function plot_data(data, IJK, key)
        [ni,nj,~,~]=size(data);
        I=IJK(1); J=IJK(2); K=IJK(3);
        ts1=squeeze(data(I,J,K,:));
        if strcmpi(key,'m')
            
            subplot(3,3,1);
            if I-1>0 && J+1<=nj 
                plot(ts1,'r'); hold on; plot(squeeze(data(I-1,J+1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I-1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,2);
            if J+1<=nj 
                plot(ts1,'r'); hold on; plot(squeeze(data(I,J+1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,3);
            if I+1<=ni && J+1<=nj 
                plot(ts1,'r'); hold on; plot(squeeze(data(I+1,J+1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,4);
            if I-1>0 
                plot(ts1,'r'); hold on; plot(squeeze(data(I-1,J,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,5);
            plot(squeeze(data(I,J,K,:)),'r');
            title(['(' num2str(I) ',' num2str(J) ',' num2str(K) ')']);
            
            subplot(3,3,6);
            if I+1<=ni, 
                plot(ts1,'r');hold on; plot(squeeze(data(I+1,J,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,7);
            if I-1>0 && J-1>0 
                plot(ts1,'r');hold on; plot(squeeze(data(I-1,J-1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,8);
            if J-1>0 
                plot(ts1,'r');hold on; plot(squeeze(data(I,J-1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J+1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
            subplot(3,3,9);
            if I+1<=ni && J-1>0 
                plot(ts1,'r');hold on; plot(squeeze(data(I+1,J-1,K,:)));
                cc=corrcoef(ts1,squeeze(data(I+1,J-1,K,:)));cc=cc(1,2);
            end
            title(['cc=' num2str(cc) '']);
            
        else
            plot(ts1);
        end
    end

    %% prepare for the correlation tensor figure
    function prep_tensor_figure(key)
        if ishghandle(h_tens)
            set(0,'CurrentFigure',h_tens);
        else
            h_tens=figure;
        end
        set(gcf,'color',get(gcf,'defaultFigureColor'));
        clf;
    end

    %% plot correlational tensor
    
    %% nested function: plot local correlation tensor
    function plot_corr_tensor(img,IJK,key)
        % calculate local correlation
        [nx,ny,nz,~]=size(data);
        lcc=zeros(26,1);
        if I-1>0 && J-1>0 && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J-1, K-1,:)));lcc(1)=R(1,2);
        end
        if J-1>0 && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J-1, K-1,:)));lcc(2)=R(1,2);
        end
        if I+1<=nx && J-1>0 && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J-1, K-1,:)));lcc(3)=R(1,2);
        end
        if I-1>0 && K-1>0 
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J,   K-1,:)));lcc(4)=R(1,2);
        end
        if K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J,   K-1,:)));lcc(5)=R(1,2);
        end
        if I+1<=nx && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J,   K-1,:)));lcc(6)=R(1,2);
        end
        if I-1>0 && J+1<=ny && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J+1, K-1,:)));lcc(7)=R(1,2);
        end
        if J+1<=ny && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J+1, K-1,:)));lcc(8)=R(1,2);
        end
        if I+1<=nx && J+1<=ny && K-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J+1, K-1,:)));lcc(9)=R(1,2);
        end
        if I-1>0 && J-1>0 
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J-1, K,:)));lcc(10)=R(1,2);
        end
        if J-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J-1, K,:)));lcc(11)=R(1,2);
        end
        if I+1<=nx && J-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J-1, K,:)));lcc(12)=R(1,2);
        end
        if I-1>0
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J,   K,:)));lcc(13)=R(1,2);
        end
        if I+1<=nx
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J,   K,:)));lcc(14)=R(1,2);
        end
        if I-1>0 && J+1<=ny
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J+1, K,:)));lcc(15)=R(1,2);
        end
        if J+1<=ny
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J+1, K,:)));lcc(16)=R(1,2);
        end
        if I+1<=nx && J+1<=ny
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J+1, K,:)));lcc(17)=R(1,2);
        end
        if I-1>0 && J-1>0 && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J-1, K+1,:)));lcc(18)=R(1,2);
        end
        if J-1>0 && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J-1, K+1,:)));lcc(19)=R(1,2);
        end
        if I+1<=nx && J-1>0 && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J-1, K+1,:)));lcc(20)=R(1,2);
        end
        if I-1>0 && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J,   K+1,:)));lcc(21)=R(1,2);
        end
        if K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J,   K+1,:)));lcc(22)=R(1,2);
        end
        if I+1<=nx && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J,   K+1,:)));lcc(23)=R(1,2);
        end
        if I-1>0 && J+1<=ny && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I-1, J+1, K+1,:)));lcc(24)=R(1,2);
        end
        if J+1<=ny && K+1<=nz
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I,   J+1, K+1,:)));lcc(25)=R(1,2);
        end
        if I+1<=nx && J+1<=ny && K+1<=nz 
            R=corrcoef(squeeze(img(I,J,K,:)),squeeze(img(I+1, J+1, K+1,:)));lcc(26)=R(1,2);
        end
    end

    %% nested function: show the underlay image
    function h_image = show_underlay
        [nx,ny,nz] = size(ulay);
        % assemble slices into a large "image"
        % assuming data is stored in "RPI" orientation
        uimg = zeros(ny*nrow,nx*ncol);
        for irow=1:nrow %#ok<FXUP>
            for icol=1:ncol %#ok<FXUP>
                i = (irow-1)*ncol+icol;
                xrange = (icol-1)*nx+1:icol*nx;
                yrange = (irow-1)*ny+1:irow*ny;
                if i<=nz
                    uimg(yrange,xrange) = flipud(ulay(:,:,i)');
                else
                    uimg(yrange,xrange) = flipud(zeros(ny,nx));
                end
            end    
        end
        uval = uimg;

        % convert img to indexed numbers
        if ucmax>ucmin
            uimg = round((uimg-ucmin)/(ucmax-ucmin)*(ucmapni-1)+1);
        else
            uimg = zeros(size(uimg));
        end
        uimg(uimg<1)=1; uimg(uimg>ucmapni)=ucmapni;
        uimg = ind2rgb(uimg,ucmap);

        % display the image
        subplot('Position',[0 flag_ui*ui_height/fig_height ...
                           img_width/fig_width img_height/fig_height]);

        h_image=image(uimg);
        axis equal;
        axis off;   
    end

    %% nested function: show the overlay image
    function show_overlay
        if exist('olay','var')
            hold on;
            % retrieve colormap for overlay
            if ischar(ocmap)
                if strcmpi(ocmap,'y2c') || strcmpi(ocmap,'yrbc') % yellow-red-blue-cray
                    ocmap=zeros(ocmapni,3);
                    mid=fix(ocmapni/2);
                    ocmap(1:mid,1)=linspace(0,0,mid);
                    ocmap(1:mid,2)=linspace(1,0,mid);
                    ocmap(1:mid,3)=linspace(1,1,mid);
                    ocmap(mid+1:end,1)=linspace(1,1,ocmapni-mid);
                    ocmap(mid+1:end,2)=linspace(0,1,ocmapni-mid);
                    ocmap(mid+1:end,3)=linspace(0,0,ocmapni-mid);
                elseif strcmpi(ocmap,'yrbg')                    % yellow-red-blue-green
                    ocmap=zeros(ocmapni,3);
                    mid=fix(ocmapni/2);
                    ocmap(1:mid,1)=linspace(0,0,mid);
                    ocmap(1:mid,2)=linspace(1,0,mid);
                    ocmap(1:mid,3)=linspace(0,1,mid);
                    ocmap(mid+1:end,1)=linspace(1,1,ocmapni-mid);
                    ocmap(mid+1:end,2)=linspace(0,1,ocmapni-mid);
                    ocmap(mid+1:end,3)=linspace(0,0,ocmapni-mid);
                elseif strcmpi(ocmap,'spectral')
                    ocmap=spectral(ocmapni);
                else
                    ocmap=eval([ocmap '(' num2str(ocmapni) ')']);
                end
            else
                % an N-by-3 matrix can also be loaded as the colormap
                ocmapni=size(ocmap,1);
            end
            ocmap=ocmap/max(ocmap(:));

            % colored range 
            if ~exist('oclim','var')
                a=double(max(abs(olay(:))));
                oclim=[-a a]';
            end

            % mask olay if a "mask" is provided
            if exist('mask','var')
                olay=olay.*(mask>0);
            end

            % assemble overlay slices into an image 
            oimg=zeros(ny*nrow,nx*ncol);
            for irow=1:nrow
                for icol=1:ncol
                    i = (irow-1)*ncol+icol;
                    xrange = (icol-1)*nx+1:icol*nx;
                    yrange = (irow-1)*ny+1:irow*ny;
                    if i<=nz
                        oimg(yrange,xrange) = flipud(olay(:,:,i)');
                    else
                        oimg(yrange,xrange) = flipud(zeros(ny,nx));
                    end
                end
            end
            oval=oimg;

            % convert intensities in oimg to indices in ocmap
            orig_oimg = oimg;
            % 2012-05-24, ocgap is not used anymore
            [nx,ny]=size(oimg);
            oimg=oimg(:);

            % comment on 2012-06-26
        %     if length(oclim)==2 && oclim(1)<0
        %         oclim=[oclim(1) 0 0 oclim(2)];
        %     end
            oimg_index = zeros(size(oimg));
            if length(oclim)==2
                vmin=min(oclim);
                vmax=max(oclim);
                cmin=1;
                cmax=ocmapni;
                now_what=(oimg>vmin&oimg<vmax);
                oimg_index(now_what)=...
                    round((oimg(now_what)-vmin)/(vmax-vmin)*(cmax-cmin)+cmin);
                now_what=oimg<=min(oclim);
                oimg_index(now_what)=1;
                now_what=oimg>=max(oclim);
                oimg_index(now_what)=ocmapni;
            elseif length(oclim)==4
                % negative range    
                vmin=oclim(1);
                vmax=oclim(2);
                cmin=1;
                cmax=ocmapni/2;
                now_what=(oimg>vmin&oimg<vmax);
                oimg_index(now_what)=...
                    round((oimg(now_what)-vmin)/(vmax-vmin)*(cmax-cmin)+cmin);
                now_what=(oimg<=vmin);
                oimg_index(now_what)=1;
                now_what=(oimg>=vmax&oimg<=0);
                oimg_index(now_what)=cmax;
                % positive range
                vmin=oclim(3);
                vmax=oclim(4);
                cmin=ocmapni/2+1;
                cmax=ocmapni;
                now_what=(oimg>vmin&oimg<vmax);
                oimg_index(now_what)=...
                    round((oimg(now_what)-vmin)/(vmax-vmin)*(cmax-cmin)+cmin);
                now_what=oimg>=vmax;
                oimg_index(now_what)=cmax;
                now_what=(oimg>0&oimg<vmin);
                oimg_index(now_what)=cmin;
        %         now_what=oimg<=min(oclim);
        %         oimg(now_what)=1;
        %         now_what=oimg>=max(oclim);
        %         oimg(now_what)=ocmapni;
            end
            oimg=oimg_index;
            oimg=reshape(oimg,nx,ny);
            % convert oimg to a rgb img    
            oimg=ind2rgb(oimg,ocmap);

            % set overlay transparency	
            oimgalpha=ones(size(orig_oimg))*oalpha;
            for i=1:size(uimg,1)
                for j=1:size(uimg,2)
                    if orig_oimg(i,j)<=oupper && orig_oimg(i,j)>=olower
                        oimgalpha(i,j)=0;
                    end
                end
            end
            h_image=image(oimg);
            set(h_image,'AlphaData',oimgalpha);
        end
    end

end

%%
function s = spectral(m)

%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
  0.0000 0.0000 0.0000
  0.4667 0.0000 0.5333
  0.5333 0.0000 0.6000
  0.0000 0.0000 0.6667
  0.0000 0.0000 0.8667
  0.0000 0.4667 0.8667
  0.0000 0.6000 0.8667
  0.0000 0.6667 0.6667
  0.0000 0.6667 0.5333
  0.0000 0.6000 0.0000
  0.0000 0.7333 0.0000
  0.0000 0.8667 0.0000
  0.0000 1.0000 0.0000
  0.7333 1.0000 0.0000
  0.9333 0.9333 0.0000
  1.0000 0.8000 0.0000
  1.0000 0.6000 0.0000
  1.0000 0.0000 0.0000
  0.8667 0.0000 0.0000
  0.8000 0.0000 0.0000
  0.8000 0.8000 0.8000
];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

end



