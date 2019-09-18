function iThreshold_Image

  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  %
  %  %   This scripts attempts to binarize real hemispherical images using  %
  %  %   an adaptive thresholding algorithm that can track local sky        % 
  %  %   brightness across the frame. The tool has been developed for       %
  %  %   application on images taken inside a forest.                       %
  %  %                                                                      %
  %  % VERSION                                                              %
  %  %   This script is work in progress and can be adapted as needed. It   %
  %  %   has been successfully tested on images taken during sunny and      %
  %  %   during overcast conditions.
  %  %   created 2019/09 | last modified on 2019/09/18                      %
  %  %                                                                      %
  %  % CONTRIBUTING AUTHORS                                                 %
  %  %   T.Jonas(1), C.Webster(2,1), G.Mazzotti(1), J.Malle(3,1)            %
  %  %    (1) WSL Institute for Snow and Avalanche Research SLF, Davos, CH  %
  %  %    (2) University of Edinburgh, School of GeoSciences, Edinburgh, UK %
  %  %    (3) Northumbria University, Newcastle upon Tyne, UK               %
  %  %                                                                      %
  %  % CITATION (to be updated)                                             %
  %  %   T.Jonas, C.Webster, G.Mazzotti, and J.Malle; 2019; HPEval: A       %
  %  %   canopy shortwave transmission model using high resolution          %
  %  %   hemispherical images; Agricultural and Forest Meteorology,         %
  %  %   in review                                                          %
  %  %                                                                      %
  %  % TEST DATASET                                                         %
  %  %   - hemispherical images taken around Davos Switzerland
  %  %                                                                      %
  %  % OUTPUT                                                               %
  %  %   The tool creates binarized images corresponding to the input files %
  %  %                                                                      %
  %  % METHODS                                                              %
  %  %   The scripts includes the following steps:                          %
  %  %   - import blue channel of input image                               %
  %  %   - divide image into 40x40 pixel subdomains                         %
  %  %   - evaluate the brightest pixel per subdomain                       %
  %  %   - remove brightest pixels that are not sky                         %
  %  %   - create complete map of sky brightness based on the remaining     %
  %  %     brightest pixels                                                 %
  %  %   - subtract sky brightness map from original image                  %
  %  %   - binarize adjusted image using a global threshold                 %
  %  %                                                                      %
  %  % USAGE                                                                %
  %  %   define all user settings in the source code                        %
  %  %   > iThreshold_Image                                                 %
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% USER SETTINGS

  folder_in    = 'D:\Files\Matlab Software\ForestSnow\HPEval\INPUT_IMAGES_RAW';
  folder_out   = 'D:\Files\Matlab Software\ForestSnow\HPEval\INPUT_IMAGES_BINARY';
  files        = {'REFERENCE_7858.JPG','IMG_7827.JPG','IMG_7840.JPG','IMG_7846.JPG'};
  copy_profile = 0; % set to 1 in order to use brightness profile from first image in the list for all subsequent images
  aux_plots    = 0; % set to 1 in order to display auxiliary figures
  threshold    = 205; % global threshold to be applied after removing brightness profile (the lower the number the brighter the image)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PREPARATIONS
  % create global variables
  global remgrid tungrid;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% LOOP THROUGH IMAGES
  for imix = 1:length(files)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % import image
    fname = fullfile(folder_in,files{imix});
    IMG = imread(fname);
    imgdata = double(IMG(:,:,3));
    if aux_plots
      figure;
      imagesc(imgdata);
    end;
  
    if imix == 1 || copy_profile == 0
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % create empty rasters in evaluation units
      typrad  = 20; % size of evaluation units (+/- pixels of original image)
      tunrad  = 2;  % define neighborhood (+/- evaluation units)
      maxgrid = nan(length(typrad:2*typrad:size(imgdata,1)-typrad),length(typrad:2*typrad:size(imgdata,2)-typrad)); % evaluation grid (reports brightest pixel within aggregation unit)
      tungrid = nan(length(typrad:2*typrad:size(imgdata,1)-typrad),length(typrad:2*typrad:size(imgdata,2)-typrad)); % sky grid (1 where sky, 0 else)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % loop through evaluation grid and assign maximum brightness value within evaluation unit to respective evaluation grid cell
      for xix = typrad:2*typrad:size(imgdata,2)-typrad;
        for yix = typrad:2*typrad:size(imgdata,1)-typrad;
          maxgrid(round((yix-typrad)/(2*typrad)+1),(xix-typrad)/(2*typrad)+1) = max(max(imgdata(yix-(typrad-1):yix+(typrad-1),xix-(typrad-1):xix+(typrad-1))));
        end;
      end;
      if aux_plots
        figure; ph = pcolor(maxgrid); set(ph,'linestyle','none'); set(gcf,'renderer','zbuf');
      end;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % loop through evaluation grid and identifies whether or not sky is visible in evaluation unit
      som2 = size(maxgrid,2);
      som1 = size(maxgrid,1);
      for xix = 1:som2
        for yix = 1:som1
          relpix = maxgrid(max(yix-tunrad,1):min(yix+tunrad,som1),max(xix-tunrad,1):min(xix+tunrad,som2));
          if abs(maxgrid(yix,xix) - median(relpix(:))) < 10 ...
              && maxgrid(yix,xix) > 100
            tungrid(yix,xix) = 1;
          else
            tungrid(yix,xix) = 0;
          end;
        end;
      end;
      if aux_plots
        figure; ph = pcolor(tungrid); set(ph,'linestyle','none'); set(gcf,'renderer','zbuf');
      end;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % creates updated evaluation grid with all non-sky pixels masked to 0
      remgrid = maxgrid .* tungrid;
      if aux_plots
        figure; ph = pcolor(remgrid); set(ph,'linestyle','none'); set(gcf,'renderer','zbuf');
        figure; ph = surf(remgrid);
      end;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % determins best-fit ellipsoid to match non-masked data in remgrid (which consistutes the threshold grid)
      % NOTE: this procedure was optimized for HPs taken shortly before/after sunrise/sunset under clear sky conditions
      par.x2 = 0;
      par.x1 = 0;
      par.y2 = 0;
      par.y1 = 0;
      par.cb = 125/1000;
      par.ag = 0/2/pi;
      out = fminsearch(@rmse,[par.x2,par.x1,par.y2,par.y1,par.cb,par.ag]);
      opar.x2 = out(1);
      opar.x1 = out(2);
      opar.y2 = out(3);
      opar.y1 = out(4);
      opar.cb = out(5)*1000;
      opar.ag = out(6);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % replicate best fit ellipsoid in image units
      [xgrid,ygrid] = meshgrid(1:size(imgdata,2),1:size(imgdata,1));
      xgrid = (xgrid - size(imgdata,2)/2) / 2 / typrad;
      ygrid = (ygrid - size(imgdata,1)/2) / 2 / typrad;
      [tgrid,rgrid] = cart2pol(xgrid,ygrid);
      [xgrid_rot,ygrid_rot] = pol2cart(tgrid+opar.ag,rgrid);
      viggrid = (opar.x2.*xgrid_rot.^2 + opar.x1.*xgrid_rot + opar.y2.*ygrid_rot.^2 + opar.y1.*ygrid_rot)./2 + opar.cb;
      if aux_plots
        figure; imagesc(viggrid);
      end;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold image using above ellipsoid
    % NOTE: the below internal thrshold of 285 could be adjusted if necessary
    idcorr = max(min(imgdata-viggrid+285,255),0);
    if aux_plots
      figure; imagesc(idcorr);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % binarize thresholded image
    imagegc = idcorr;
    imagegc(idcorr < threshold) = 0;
    imagegc(idcorr >= threshold) = 255;
    if aux_plots
      figure; imagesc(imagegc);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create grayscale image (RGB) from binarized image
    image_c              = [];
    image_c(:,:,1)       = imagegc./255;
    image_c(:,:,2)       = imagegc./255;
    image_c(:,:,3)       = imagegc./255;
    if aux_plots
      figure; image(image_c);
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save above grayscale image
    [pt,fn,ex] = fileparts(fname);
    imwrite(image_c,fullfile(folder_out,strcat(fn,'_bin.','png')));
  
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE SUB-FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optimization function needed in the context of fitting ellipsoid
function out = rmse(pars)

  global remgrid tungrid;
  par.x2 = pars(1);
  par.x1 = pars(2);
  par.y2 = pars(3);
  par.y1 = pars(4);
  par.cb = pars(5)*1000;
  par.ag = pars(6);
  [xgrid,ygrid] = meshgrid(1:size(remgrid,2),1:size(remgrid,1));
  xgrid = xgrid - size(remgrid,2)/2;
  ygrid = ygrid - size(remgrid,1)/2;
  [tgrid,rgrid] = cart2pol(xgrid,ygrid);
  [xgrid_rot,ygrid_rot] = pol2cart(tgrid+par.ag,rgrid);
  viggrid = (par.x2.*xgrid_rot.^2 + par.x1.*xgrid_rot + par.y2.*ygrid_rot.^2 + par.y1.*ygrid_rot)./2 + par.cb;
  difp = viggrid(logical(tungrid))-remgrid(logical(tungrid));
  sdifp = sort(difp);
  out = sqrt(mean(sdifp(1:round(0.9*length(sdifp))).^2));

end
