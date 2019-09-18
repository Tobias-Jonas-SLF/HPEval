function HPEval(sfile)

  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  % GENERAL DESCRIPTION                                                  %
  %  %   This scripts evaluates real hemispherical images to estimate local %
  %  %   sky view fraction and below-canopy shortwave radiation if the      %
  %  %   image is taken inside a forest.                                    %
  %  %                                                                      %
  %  % VERSION                                                              %
  %  %   This is a first complete version of HPEval to be references as     %
  %  %   vs 1.2 | created 2018/08 | last modified on 2019/09/14             %
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
  %  % CONTRIBUTING MATERIAL                                                %
  %  %   - solar position calculations from NOAA, found online at           %
  %  %     www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html              %
  %  %     (accessed 2018/06)                                               %
  %  %   - function to split solar radiation into direct/diffuse component  %
  %  %     as in the factorial snow model FSM 1.0 by R. Essery              %
  %  %     available at www.geosci-model-dev.net/8/3867/2015                %
  %  %     (accessed 2018/06)                                               %
  %  %   - function to convert UTM into WGS84 coordinates, found online at  %
  %  %     www.mathworks.com/matlabcentral/fileexchange/44242-utm2lonlat    %
  %  %     (accessed 2018/08)                                               %
  %  %                                                                      %
  %  % TEST DATASET                                                         %
  %  %   - hemispherical images taken around Davos Switzerland and          %
  %  %     corresponding meta data; by SLF                                  %
  %  %   - shortwave radiation data measured at locations were above images %
  %  %     were taken test images; by SLF                                   %  
  %  %   - DEM at 100 m spatial resolution, derived from TopoSwiss data     %  
  %  %     available for research; by TopoSwiss                             %
  %  %                                                                      %
  %  % OUTPUT                                                               %
  %  %   The script has a number of optional outputs, which include         %
  %  %   - reclassification of image (canopy, terrain, sky, solar track)    %
  %  %   - shortwave radiation (per horizontal, inclined, sun-facing area)  %
  %  %   - sky view fraction (for flat and/or hemispherical sensor)         %
  %  %                                                                      %
  %  % METHODS                                                              %
  %  %   The scripts includes the following steps (some of them being       %
  %  %   optional):                                                         %
  %  %   - grayscale image (according to external input on threshold and    %
  %  %     thresholding method)                                             %
  %  %   - correct image distortion (using a lens-specific profile)         %
  %  %     note: this is an essential feature and should not be skipped     %
  %  %   - mask surrounding terrain and add to image (based on DEM)         %
  %  %     note: helps to avoid misclassification of snow on the ground     %
  %  %   - calculate sky view fraction                                      %
  %  %   - calculate position of sun (according to formulas by NOAA)        %
  %  %   - plot sun disk on image behind canopy / terrain elements to       %
  %  %     evaluate the transmissivity of direct radiation through the      %
  %  %     canopy                                                           %
  %  %   - split swr data measured at open site into diffuse/direct         %
  %  %     radiation (according to formulas in FSM 1.0 by R. Essery)        %
  %  %     note: potential swr is used if no measurements are available     %
  %  %   - calculate sub-canopy swr based on diffuse radiation times sky    %
  %  %     view fraction plus direct radiation times transmissivity         %
  %  %                                                                      %
  %  % NOTATION                                                             %
  %  %   - local spherical coordinates are given as [rad,tht,phi], where    % 
  %  %     rad is radial distance in [m]                                    %
  %  %     tht is zenith angle (theta) in degree, i.e.                      %
  %  %       sun @ zenith is tht = 0;                                       %
  %  %       sun @ horizon is tht = 90                                      %
  %  %     phi is azimuthal angle (phi) in degree, where                    %
  %  %       0° is North                                                    %
  %  %       90° is East                                                    %
  %  %       180° is South                                                  %
  %  %       270° is West                                                   %
  %  %       note that hemispherical images are taken upwards so that       %
  %  %       azimuthal angles increase in counterclockwise direction,       %
  %  %       i.e. South is up if West is left                               %
  %  %                                                                      %
  %  % SETTINGS                                                             %
  %  %   All settings such as paths to I/O data and model parameters are    %
  %  %   handled through an external settings file the path/filename of     %
  %  %   which is the only input argument to this script                    %
  %  %                                                                      %
  %  % USAGE                                                                %
  %  %   > HPEval('HPEval_Settings.m')                                      %
  %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% LOAD SETTINGS
  %% > read settings file
  if nargin == 0
    spath          = which('HPEval.m');
    sfile          = fullfile(spath(1:end-8),'HPEval_Settings.m');
  end;
  cdir = cd;
  [setp,setf,~]    = fileparts(sfile);
  if ~isempty(setp)
    cd(setp);
  end;
  try
    hpeset         = feval(setf);
    cd(cdir);
  catch
    cd(cdir);
    error('Error while reading HPEval settings file');
  end;
  try
    version        = hpeset.version;
  catch
    error('Error while reading HPEval settings file');
  end;
  
  %% > update settings if necessary
  if versioncompare('1.1',version) == 1
    % hpeset.in.swr used to point to a specific file, which is now specified in the info file. This is why from now on hpeset.in.swr should only contain a path 
    [setp,setf,setx] = fileparts(hpeset.in.swr);
    if ~isempty(setx)
      hpeset.in.swr = setp;
    end;
  end;
  if versioncompare('1.2',version) == 1
    % adds new user settings, not available in previous versions 
    hpeset.par.zenith_rings = 9;
    hpeset.for.dat = {'mat'};
  end;
  
  % specify earliest version under which settings file does not require update
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PREPARATIONS
  
  %% > prepare access to data and tools
  % progress info
  if hpeset.sw.report
    disp('   -------------------------------------------------------')
    disp(['   ' datestr(now,'HH:MM:SS') ' > HPEval started'])
    drawnow;
  end;
  
  % create output folder if not existent 
  if ~exist(hpeset.out.img,'dir')
    mkdir(hpeset.out.img)
  end
  if ~exist(hpeset.out.fig,'dir')
    mkdir(hpeset.out.fig)
  end
  if ~exist(hpeset.out.dat,'dir')
    mkdir(hpeset.out.dat)
  end
  
  %% > load and extract image specifications
  % load data
  try
    fid = fopen(hpeset.in.inf,'r');
    formatstr = '%s%f%f%f%f%f%f%f%s%s%f%s%f';
    fgetl(fid);
    data = textscan(fid,formatstr,'delimiter','\t');
    img.filename             = data{1};                                        %
    img.threshold            = data{2};                                        % value between 0 and 255 (assuming 8-bit RGB image). An  image is interpreted black <  threshold, and white otherwise. Use 1 in case of externally binarized images (1-bit B/W images).
    img.gs_method            = data{3};                                        % value between 1 and 3 (supported options are detailed below, c.f. 'grayscale image')
    img.image_radius         = data{4};                                        % pixels (radius of image circle that represents a 180° field of view). Note this value will be disregarded should a specific lens profile be used
    img.x                    = data{5};                                        % easting (using coordinate system specified in HPEval_Settings.m)
    img.y                    = data{6};                                        % northing (using coordinate system specified in HPEval_Settings.m)
    img.camera_height        = data{7};                                        % m above DEM
    img.rotation_error       = data{8};                                        % °. Use positive values to rotate coordinate system counterclockwise (= prepone individual sunflecks), use negative values to rotate coordinate system clockwise (= postpone individual sunflecks)
    for imix = 1:length(img.filename)
      if isempty(data{9}{imix}) || isempty(data{10}{imix})
        img.t1(imix)         = NaN;                                      
        img.t2(imix)         = NaN;
        img.dt(imix)         = NaN;
      else
        try 
          img.t1(imix)       = datenum(data{9}{imix},'dd.mm.yyyy HH:MM:SS');   % Matlab time
          img.t2(imix)       = datenum(data{10}{imix},'dd.mm.yyyy HH:MM:SS');  % Matlab time
          img.dt(imix)       = data{11}(imix)/60/24;                           % split day
        catch
          img.t1(imix)       = datenum(data{9}{imix},'dd.mm.yyyy HH:MM');      % Matlab time
          img.t2(imix)       = datenum(data{10}{imix},'dd.mm.yyyy HH:MM');     % Matlab time
          img.dt(imix)       = data{11}(imix)/60/24;                           % split day
        end;
      end;
    end;
    img.swr                  = data{12};                                       % 
    img.swr_scaling          = data{13};                                       % applied as a multiplier to observed SWR data. This is to account for eventual sensor bias (detemined by external cross-calibation)
    img.swr_scaling(isnan(img.swr_scaling)) = 1;
    img.no_ims               = length(img.filename);
    fclose(fid);
  catch
    try
      fclose(fid);
    end;
    error('Error while reading image specification file');
  end;
  
  % test integrity of data
  if img.no_ims < 1 ...
  || img.no_ims ~= length(img.threshold) ...
  || img.no_ims ~= length(img.gs_method) ...  
  || img.no_ims ~= length(img.image_radius) ...  
  || img.no_ims ~= length(img.x) ...  
  || img.no_ims ~= length(img.y) ...  
  || img.no_ims ~= length(img.camera_height) ...  
  || img.no_ims ~= length(img.rotation_error) ...  
  || img.no_ims ~= length(img.t1) ...  
  || img.no_ims ~= length(img.t2) ...  
  || img.no_ims ~= length(img.dt) ...
  || img.no_ims ~= length(img.swr) ...
  || img.no_ims ~= length(img.swr_scaling)
    error('Error while gathering data from image specification file');
  end;
  
  % identify if a refence image is available, and separate info from main set of image metadata
  ref = [];
  for imix = 1:img.no_ims
    if length(img.filename{imix}) >= 9 && strcmpi(img.filename{imix}(1:9),'reference')
      ref.filename           = img.filename{imix};
      ref.threshold          = img.threshold(imix);
      ref.gs_method          = img.gs_method(imix);
      ref.image_radius       = img.image_radius(imix);
      ref.x                  = img.x(imix); 
      ref.y                  = img.y(imix); 
      ref.camera_height      = img.camera_height(imix);
      ref.rotation_error     = img.rotation_error(imix);
      ref.t1                 = img.t1(imix);
      ref.t2                 = img.t2(imix);
      ref.dt                 = img.dt(imix);
      ref.swr                = img.swr{imix};
      ref.swr_scaling        = img.swr_scaling(imix);
      img.filename           = img.filename([1:imix-1,imix+1:end]);
      img.threshold          = img.threshold([1:imix-1,imix+1:end]);
      img.gs_method          = img.gs_method([1:imix-1,imix+1:end]);
      img.image_radius       = img.image_radius([1:imix-1,imix+1:end]);
      img.x                  = img.x([1:imix-1,imix+1:end]); 
      img.y                  = img.y([1:imix-1,imix+1:end]); 
      img.camera_height      = img.camera_height([1:imix-1,imix+1:end]); 
      img.rotation_error     = img.rotation_error([1:imix-1,imix+1:end]);
      img.t1                 = img.t1([1:imix-1,imix+1:end]); 
      img.t2                 = img.t2([1:imix-1,imix+1:end]); 
      img.dt                 = img.dt([1:imix-1,imix+1:end]);
      img.swr                = img.swr([1:imix-1,imix+1:end]);
      img.swr_scaling        = img.swr_scaling([1:imix-1,imix+1:end]);
      img.no_ims             = img.no_ims - 1;
      break
    end;
  end;
  
  % identify if a mask image is available, and separate info from main set of image metadata
  mask = [];
  for imix = 1:img.no_ims
    if length(img.filename{imix}) >= 4 && strcmpi(img.filename{imix}(1:4),'mask')
      mask.filename          = img.filename{imix};
      mask.threshold         = img.threshold(imix);
      mask.gs_method         = img.gs_method(imix);
      img.filename           = img.filename([1:imix-1,imix+1:end]);
      img.threshold          = img.threshold([1:imix-1,imix+1:end]);
      img.gs_method          = img.gs_method([1:imix-1,imix+1:end]);
      img.image_radius       = img.image_radius([1:imix-1,imix+1:end]);
      img.x                  = img.x([1:imix-1,imix+1:end]); 
      img.y                  = img.y([1:imix-1,imix+1:end]); 
      img.camera_height      = img.camera_height([1:imix-1,imix+1:end]); 
      img.rotation_error     = img.rotation_error([1:imix-1,imix+1:end]);
      img.t1                 = img.t1([1:imix-1,imix+1:end]); 
      img.t2                 = img.t2([1:imix-1,imix+1:end]); 
      img.dt                 = img.dt([1:imix-1,imix+1:end]);
      img.swr                = img.swr([1:imix-1,imix+1:end]);
      img.swr_scaling        = img.swr_scaling([1:imix-1,imix+1:end]);
      img.no_ims             = img.no_ims - 1;
      break
    end;
  end;
  
  %% > load and extract DEM data
  % load data (to be provided in ascii format)
  if (hpeset.sw.terrain || hpeset.sw.swr) && ~isempty(hpeset.in.dem)
    dem                      = load_ascii_grid(hpeset.in.dem);  
    dem.xvec                 = (dem.xllcorner:dem.cellsize:dem.xllcorner+dem.cellsize*(dem.ncols-1))+dem.cellsize/2 ; % refers to center of cell
    dem.yvec                 = (dem.yllcorner:dem.cellsize:dem.yllcorner+dem.cellsize*(dem.nrows-1))+dem.cellsize/2'; % refers to center of cell
    [dem.xmat dem.ymat]      = meshgrid(dem.xvec,dem.yvec);                % refers to center of cell
    dem.zmat                 = dem.data;                                   % refers to center of cell
    dem.zmat(abs(dem.zmat-dem.NODATA_value)<eps) = NaN;
    dem                      = rmfield(dem,'data');
    
    % calculate slope and aspect from DEM
    [fx fy]                  = gradient(dem.zmat,dem.cellsize);
    [~,rad]                  = cart2pol(fx,fy);
    dem.slp                  = atan(rad)*180/pi;                           % slp is the slope angle in degree, i.e. flat is slp = 0; a vertical face is slp = 90
    dem.asp                  = atan2(-fx,-fy)*180/pi;
    dem.asp(dem.asp<0)       = 360 + dem.asp(dem.asp<0);                   % asp is the aspect angle in degree, so that 0° is North, 90° is East, 180° is South, and 270° is East                                      
  elseif hpeset.sw.terrain && isempty(hpeset.in.dem)
    error('Error: a DEM is needed to mask surrounding terrain. Either provide DEM or set hpeset.sw.terrain = 0');
  end;
  
  %% > load and extract SWR data from reference site
  % load data (to be provided in ascii format, c.f. example file)
  if ~isempty(ref) && ~isempty(ref.swr)
    try
      fid                    = fopen(fullfile(hpeset.in.swr,ref.swr),'r');
      formatstr              = '%s%f';
      fgetl(fid);
      data                   = textscan(fid,formatstr,'delimiter','\t');
      try
        swr.time             = datenum(data{1},'dd.mm.yyyy HH:MM:SS');
      catch
        swr.time             = datenum(data{1},'dd.mm.yyyy HH:MM');
      end;
      swr.data               = data{2}*ref.swr_scaling;
      fclose(fid);
    catch
      try
        fclose(fid);
      end;
      error('Error while reading swr data file');
    end;
  else
    swr                      = [];
  end;

  %% > load and extract image mask
  if ~isempty(mask)
    % read image file
    ipath                    = fileparts(hpeset.in.img);
    mask_gs                  = double(imread(fullfile(ipath,mask.filename)));
    gs_method                = mask.gs_method;
    threshold                = mask.threshold;
    
    % grayscale image
    switch gs_method
    case 1
      % grayscale image using blue channel
      mask_gs                = mask_gs(:,:,3);
    case 2
      % grayscale image using luminance
      mask_gs                = 0.2989*mask_gs(:,:,1) + 0.58705*mask_gs(:,:,2) + 0.11405*mask_gs(:,:,3);
    case 3
      % grayscale image using mean of RGB channels / or use to handle image that only features one channel
      mask_gs                = mean(mask_gs,3);
    otherwise
      error('Error while gray-scaling image, unknown option');
    end;
   
    % binarize image
    mask_gs(mask_gs <  threshold) = 0;
    mask_gs(mask_gs >= threshold) = 1;
  end;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CALCULATIONS
  
  %% > loop through image files (start with reference image if exist)
  if isempty(ref)
    stix = 1;
  else
    stix = 0;
  end;
  for imix = stix:length(img.filename)
    
    %% > read and binarize image
    % read image file
    if imix == 0
      ipath                  = fileparts(hpeset.in.img);
      image_gs               = double(imread(fullfile(ipath,ref.filename)));
      gs_method              = ref.gs_method;
      threshold              = ref.threshold;
      camera_height          = ref.camera_height;
      rotation_error         = ref.rotation_error;
      image_radius           = ref.image_radius;      
      loc_x                  = ref.x;
      loc_y                  = ref.y;
    else
      ipath                  = fileparts(hpeset.in.img);
      image_gs               = double(imread(fullfile(ipath,img.filename{imix})));
      gs_method              = img.gs_method(imix);
      threshold              = img.threshold(imix);
      camera_height          = img.camera_height(imix);
      rotation_error         = img.rotation_error(imix);
      image_radius           = img.image_radius(imix);      
      loc_x                  = img.x(imix);
      loc_y                  = img.y(imix);
      val_data               = img.swr{imix};
      swr_scaling            = img.swr_scaling(imix);
    end;
    
    % grayscale image
    switch gs_method
    case 1
      % grayscale image using blue channel
      image_gs               = image_gs(:,:,3);
    case 2
      % grayscale image using luminance
      image_gs               = 0.2989*image_gs(:,:,1) + 0.58705*image_gs(:,:,2) + 0.11405*image_gs(:,:,3);
    case 3
      % grayscale image using mean of RGB channels / or use to handle image that only features one channel
      image_gs               = mean(image_gs,3);
    otherwise
      error('Error while gray-scaling image, unknown option');
    end;
   
    % binarize image
    image_gs(image_gs <  threshold) = 0;
    image_gs(image_gs >= threshold) = 1;
    
    %% > prepare image for adding features
    % create lens profile
    switch hpeset.par.lens_profile
    case 0
      % this is the default profile assuming a perfect equi-angular lens
      lens_profile_tht       = [0:10:90];
      lens_profile_rpix      = [0:1/9:1]; 
      image_center           = [];                                         % leave empty to use image dimensions devided by 2
    case 1
      % this is the profile for our Sigma 4.5mm F2.8 EX DC HSM Circular Fisheye (s/n 12505220)
      lens_profile_tht       = [0:10:90];
      lens_profile_rpix      = [0, 0.1257, 0.2497, 0.3729, 0.4916, 0.6072, 0.7176, 0.8217, 0.9175, 1];
      image_radius           = 1436;                                       % overrides manual input, dermined using / assuming Canon EOS 600D camera
      image_center           = [];                                         % dermined using / assuming Canon EOS 600D camera (long axis / short axis, both in pixels)
    case 2
      % this is the profile for our Yasuhara Madoka 180 7.3mm F4 Circular Fisheye (s/n 02692) - this is the newer good lens
      lens_profile_tht       = [0:10:90];
      lens_profile_rpix      = [0, 0.1671, 0.3305, 0.4833, 0.6230, 0.7450, 0.8450, 0.9198, 0.9704, 1];
      image_radius           = 1571;                                       % overrides manual input, dermined using / assuming Sony Nex-6 camera
      image_center           = [2407, 1590];                               % dermined using / assuming Sony Nex-6 camera (long axis / short axis, both in pixels)
    case 3
      % this is the profile for our Yasuhara Madoka 180 7.3mm F4 Circular Fisheye (s/n 00749) - this is the old decentered lens
      lens_profile_tht       = [0:10:90];
      lens_profile_rpix      = [0, 0.1671, 0.3305, 0.4833, 0.6230, 0.7450, 0.8450, 0.9198, 0.9704, 1];
      image_radius           = 1571;                                       % overrides manual input, dermined using / assuming Sony Nex-6 camera
      image_center           = [2435, 1547];                               % dermined using / assuming Sony Nex-6 camera (long axis / short axis, both in pixels)
    otherwise
      error('Error while preparing lens profile, unknown option');
    end;
    if isempty(image_center)
      image_center           = fliplr(size(image_gs)/2);                 
    end;
    
    % create pixel matrices for radius, phi, and theta
    [xgrid,ygrid]            = meshgrid(1:size(image_gs,2),1:size(image_gs,1));
    rgrid                    = sqrt((xgrid - image_center(1)).^2 + (ygrid - image_center(2)).^2);  % radius from images center in pixels
    tgrid                    = interp1(lens_profile_rpix*image_radius,lens_profile_tht,rgrid(:));
    tgrid                    = reshape(tgrid,size(rgrid,1),size(rgrid,2));                         % zenith angle (theta) in degree, i.e. zenith is phi = 0; horizon is phi = 90
    pgrid                    = atan2(ygrid - image_center(2),xgrid - image_center(1)); 
    pgrid                    = 90-rad2deg(pgrid);
    pgrid(pgrid<0)           = pgrid(pgrid<0) + 360;                                               % azimuth angle (phi) in degree, so that 0° is North, 90° is East, 180° is South, and 270° is West. Use image(pgrid/6) to verify, i.e. not pcolor(pgrid). Note that a hemispherical image has flipped compass headings as it is taken upwards (i.e. South is up if West is left).
    
    % mask irrelevant parts of the image
    delete_pix               = rgrid > image_radius;
    image_gs(delete_pix)     = 0.5;
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TERRAIN SHADING CALCULATIONS
    if hpeset.sw.terrain
    
      %% > calculate horizon line to acount for terrain shading
      if hpeset.sw.terrain < 3 || imix == stix || (hpeset.sw.terrain > 4 && imix == 1)
        % transfer DEM to spherical coordinates relative to camera position
        dem.xpc                        = dem.xmat - loc_x;
        dem.ypc                        = dem.ymat - loc_y;
        dem.zpc                        = dem.zmat - interp2(dem.xmat,dem.ymat,dem.zmat,loc_x,loc_y) - camera_height;
        [dem.phi,dem.tht,dem.rad]      = cart2sph(dem.xpc(:),dem.ypc(:),dem.zpc(:));
        dem.phi                        = (pi/2) - dem.phi;
        dem.phi                        = dem.phi * (180/pi);
        dem.phi                        = dem.phi + rotation_error;         % rotate terrain to acount for rotation error in orientation of image
        dem.phi(dem.phi>=360)          = dem.phi(dem.phi>=360) - 360;              
        dem.phi(dem.phi<0)             = dem.phi(dem.phi<0) + 360;         % phi now is the azimuth angle (phi) in degree, so that 0° is North, 90° is East, 180° is South, and 270° is West
        dem.tht                        = (pi/2) - dem.tht;
        dem.tht                        = dem.tht * (180/pi);
        dem.tht(dem.tht >= 90,1)       = 90;                               % tht now is the zenith angle (theta) in degree, i.e. zenith is phi = 0; horizon is phi = 90

        % crop DEM data outside of a horizontal radius of hpeset.par.terrain_perim
        keepix                         = find(dem.rad.*sin(dem.tht/360*2*pi) < hpeset.par.terrain_perim);
        dem.phi                        = dem.phi(keepix);
        dem.tht                        = dem.tht(keepix);
        dem.rad                        = dem.rad(keepix);
        
        % calculate horizon line
        hpeset.sw.horizon              = 2;                                % hardcoded switch to reflect recommended choice
        switch hpeset.sw.horizon
        case 1 % simple approach, may need some additional smoothing if applied at too high azimuthal resolution / too coarse DEM
          phibins                      = [180/hpeset.par.terrain_angs:360/hpeset.par.terrain_angs:360-180/hpeset.par.terrain_angs]';
          thtmin                       = ones(length(phibins),1).*90;                     
          for pbix = 1:length(phibins)
            fix                        = find(dem.phi >= phibins(pbix) - 180/hpeset.par.terrain_angs & dem.phi < phibins(pbix) + 180/hpeset.par.terrain_angs);
            if isempty(fix)
              thtmin(pbix)             = NaN;
            else
              thtmin(pbix)             = min(dem.tht(fix));
            end;
          end;
          horizon.phi                  = phibins;
          kix                          = find(~isnan(thtmin));         
          horizon.tht                  = min(90,interp1q([phibins(kix)-360;phibins(kix);phibins(kix)+360],[thtmin(kix);thtmin(kix);thtmin(kix)],phibins)); 
        case 2 % more accurate approach, note: could be used to additionally record distance to horizon)
          radbins                      = [dem.cellsize*4:sqrt(2)*dem.cellsize:hpeset.par.terrain_perim];
          phibins                      = [180/hpeset.par.terrain_angs:360/hpeset.par.terrain_angs:360-180/hpeset.par.terrain_angs]';
          thtmin                       = ones(length(phibins),1).*90;                     
          for rbix = length(radbins)-1:-1:1
            fix1                       = find(dem.rad >= radbins(rbix) & dem.rad < radbins(rbix+1));
            [~,fix2]                   = sort(dem.phi(fix1));
            thtmin                     = min(thtmin,interp1q([dem.phi(fix1(fix2(end)))-360;dem.phi(fix1(fix2));dem.phi(fix1(fix2(1)))+360],[dem.tht(fix1(fix2(end)));dem.tht(fix1(fix2));dem.tht(fix1(fix2(1)))],phibins));
          end;
          horizon.phi                  = phibins;
          horizon.tht                  = thtmin;
        otherwise
          error('Error while calculating horizon line, unknown option');
        end;
      end;

      %% > add terrain shading to image
      if (hpeset.sw.terrain == 4 && imix == stix) || (hpeset.sw.terrain == 6 && imix == 1)
        % calc terrain mask and plot on top of canopy elements
        if hpeset.par.terrain_angs > 360                                   % plot at native resolution
          for pix = 1:length(horizon.phi)
            fix              = pgrid >= horizon.phi(pix) - 180/hpeset.par.terrain_angs & pgrid < horizon.phi(pix) + 180/hpeset.par.terrain_angs & tgrid >= horizon.tht(pix);
            image_gs(fix)    = 0.2;
          end;
        else                                                               % resample horizon line to plot at least at 1° resolution
          horizon.rphi       = [180/360:360/360:360-180/360]';
          horizon.rtht       = interp1([horizon.phi(end)-360; horizon.phi; horizon.phi(1)+360],horizon.tht([end,1:end,1]),horizon.rphi);
          for pix = 1:360
            fix              = pgrid >= horizon.rphi(pix) - 180/360 & pgrid < horizon.rphi(pix) + 180/360 & tgrid >= horizon.rtht(pix);
            image_gs(fix)    = 0.2;
          end;
        end;
        ter.mask             = image_gs == 0.2;
        ter.size             = size(image_gs);
        
      elseif (hpeset.sw.terrain == 4 && imix > stix) || (hpeset.sw.terrain == 6 && imix > 1)
        % copy-paste terrain mask from first relevant image
        if any(abs(ter.size - size(image_gs)) > 0)
          error('Error while inheriting terrain mask between images due to different image sizes. In this case option hpeset.sw.terrain == 4 will not work, please change this setting');
        end;
        image_gs(ter.mask)   = 0.2;
        
      else
        % calc terrain mask and plot below canopy elements
        if hpeset.par.terrain_angs > 360 || hpeset.sw.terrain == 2         % plot at native resolution
          for pix = 1:length(horizon.phi)
            fix              = pgrid >= horizon.phi(pix) - 180/hpeset.par.terrain_angs & pgrid < horizon.phi(pix) + 180/hpeset.par.terrain_angs & tgrid >= horizon.tht(pix);
            image_gs(fix)    = min(image_gs(fix),0.2);
          end;
        else                                                               % resample horizon line to plot at least at 1° resolution
          horizon.rphi       = [180/360:360/360:360-180/360]';
          horizon.rtht       = interp1([horizon.phi(end)-360; horizon.phi; horizon.phi(1)+360],horizon.tht([end,1:end,1]),horizon.rphi);
          for pix = 1:360
            fix              = pgrid >= horizon.rphi(pix) - 180/360 & pgrid < horizon.rphi(pix) + 180/360 & tgrid >= horizon.rtht(pix);
            image_gs(fix)    = min(image_gs(fix),0.2);
          end;
        end;
      end;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% APPLY IMAGE MASK
    if ~isempty(mask)
      if any(abs(size(mask_gs) - size(image_gs)) > 0)
        error('Error while applying mask due to different image sizes');
      end;
      image_gs(mask_gs == 0) = 0.3;
    end;
    image_gs(delete_pix)     = 0.5;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SKYVIEW FRACTION CALCULATIONS
    if hpeset.sw.svf || hpeset.sw.swr
    
      %% > calculate sky view fraction
      % loop through zenith ring
      ring_tht                         = 0:90/hpeset.par.zenith_rings:90;                                            % zenith angle (theta) in degree, i.e. zenith is phi = 0; horizon is phi = 90
      ring_rad                         = interp1(lens_profile_tht,lens_profile_rpix*image_radius,ring_tht);          % equivalent radius in pixels that match above zenith angles

      % loop through zenith angle rings
      white_to_all_ratio               = nan(1,length(ring_rad)-1);
      white_to_all_ratio_S             = nan(1,length(ring_rad)-1);
      white_to_all_ratio_N             = nan(1,length(ring_rad)-1);
      surface_area_ratio_hemi          = nan(1,length(ring_rad)-1);
      surface_area_ratio_flat          = nan(1,length(ring_rad)-1);
      for rix = 1:length(ring_rad)-1

        % identify pixels within zenith ring
        inner_rad                      = ring_rad(rix);                                                              % inner radius of zenith ring
        outer_rad                      = ring_rad(rix+1);                                                            % outer radius of zenith ring
        relevant_pix                   = find(rgrid > inner_rad & rgrid <= outer_rad);                               % pixels within zenith ring
        relevant_pix_S                 = find(rgrid > inner_rad & rgrid <= outer_rad & pgrid > 90 & pgrid <= 270);   % pixels within zenith ring and southerly aspects  
        relevant_pix_N                 = find(rgrid > inner_rad & rgrid <= outer_rad & (pgrid <= 90 | pgrid > 270)); % pixels within zenith ring and northerly aspects

        % calculate tranmissivity and respective weights per ring 
        white_to_all_ratio(rix)        = sum(image_gs(relevant_pix) == 1) / length(relevant_pix);                    % raw transmissivity per ring
        white_to_all_ratio_S(rix)      = sum(image_gs(relevant_pix_S) == 1) / length(relevant_pix_S);                % raw transmissivity per ring for southerly aspects
        white_to_all_ratio_N(rix)      = sum(image_gs(relevant_pix_N) == 1) / length(relevant_pix_N);                % raw transmissivity per ring for northerly aspects
        surface_area_ratio_hemi(rix)   = cos(ring_tht(rix)/360*2*pi) - cos(ring_tht(rix+1)/360*2*pi);                % surface area*  of zenith ring per surface area*  of the entire hemisphere, * areas on hemisphere
        surface_area_ratio_flat(rix)   = sin(ring_tht(rix+1)/360*2*pi)^2 - sin(ring_tht(rix)/360*2*pi)^2;            % surface area** of zenith ring per surface area** of the entire hemisphere, ** areas projected on horizontal surface  
      end;

      % calculate sky-view fraction (SVF)
      if imix == 0
        ref.SVF_hemi                   = sum(white_to_all_ratio.*surface_area_ratio_hemi);                           % SVF to reflect perspective of hemipherically shaped sensor surface or plant (weights zenith rings according to their surface area on the hemisphere)
        ref.SVF_flat                   = sum(white_to_all_ratio.*surface_area_ratio_flat);                           % SVF to reflect perspective of a horizonatal flat uplooking sensor surface (weights zenith rings according to their surface area projected onto a horizonatal flat surface)
        ref.SVF_flat_S                 = sum(white_to_all_ratio_S.*surface_area_ratio_flat);                         % SVF for southerly aspects only
        ref.SVF_flat_N                 = sum(white_to_all_ratio_N.*surface_area_ratio_flat);                         % SVF for northerly aspects only
      else
        img.SVF_hemi(imix)             = sum(white_to_all_ratio.*surface_area_ratio_hemi);                           % SVF to reflect perspective of hemipherically shaped sensor surface or plant (weights zenith rings according to their surface area on the hemisphere)
        img.SVF_flat(imix)             = sum(white_to_all_ratio.*surface_area_ratio_flat);                           % SVF to reflect perspective of a horizonatal flat uplooking sensor surface (weights zenith rings according to their surface area projected onto a horizonatal flat surface)
        img.SVF_flat_S(imix)           = sum(white_to_all_ratio_S.*surface_area_ratio_flat);                         % SVF for southerly aspects only
        img.SVF_flat_N(imix)           = sum(white_to_all_ratio_N.*surface_area_ratio_flat);                         % SVF for northerly aspects only
      end;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SWR CALCULATIONS
    if hpeset.sw.swr
    
      if imix == 0 && isempty(ref.swr)
        error('Error while calculating solar track: if a reference image is listed, associated swr data must be available');
      end;
      
      if imix > 0 || ~isempty(ref.swr)
        
        if hpeset.sw.swr < 2 || imix == stix
          %% > calculate solar track
          % the below calcuations are based on formulas found at https://www.esrl.noaa.gov/gmd/grad/solcalc/calcdetails.html (accessed June 2018)

          % define time axis for calculations
          if imix == 0
            if hpeset.sw.swr < 2 
              loc_time       = swr.time;
            else
              loc_time       = (ref.t1:ref.dt:ref.t2)';
            end;
          else
            loc_time         = (img.t1(imix):img.dt(imix):img.t2(imix))';
          end;
          if hpeset.sw.swr >= 2 && imix == stix && any(isnan(loc_time))
            error('Error while calculating solar track, time axis information are missing');
          end;

          % conversion of x/y  to lat/lon
          switch hpeset.par.coor_system
          case 'CH1903'
            xd               = (loc_x - 600000)/1000000;
            yd               = (loc_y - 200000)/1000000;
            lon              = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
            lat              = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;
          case 'CH1903+'
            xd               = (loc_x - 2600000)/1000000;
            yd               = (loc_y - 1200000)/1000000;
            lon              = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
            lat              = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;
          case 'UTM'
            zone             = sscanf(hpeset.par.utm_zone,'%f%c');
            hemi             = char(zone(2));
            [lat,lon]        = utm2deg(loc_x,loc_y,zone(1),hemi);
          otherwise
            error('Error while converting x/y to lat/lon, the coordinate system is not yet supported');
          end;

          % conversion of time
          time_vec                     = datevec(loc_time);                % time stamps should represent an instantaneous sample
          day_of_year                  = (datenum(time_vec(:,1),time_vec(:,2),time_vec(:,3)));
          time_of_day                  = (datenum(00,00,00,time_vec(:,4),time_vec(:,5),time_vec(:,6)));
          julian_day                   = day_of_year + 1721058.5 + time_of_day - hpeset.par.time_zone/24;       
          julian_century               = (julian_day - 2451545) ./ 36525;

          % calculate solar elevation angle
          geom_mean_long_sun_deg       = mod(280.46646 + julian_century .* (36000.76983 + julian_century .* 0.0003032),360);
          geom_mean_anom_sun_deg       = 357.52911 + julian_century .* (35999.05029 - 0.0001537 .* julian_century);
          eccent_earth_orbit           = 0.016708634 - julian_century .* (0.000042037 + 0.0000001267 .* julian_century);
          sun_eq_of_ctr                = sin(deg2rad(geom_mean_anom_sun_deg)) .* (1.914602 - julian_century .* (0.004817 + 0.000014 .* julian_century)) + sin(deg2rad(2 .* geom_mean_anom_sun_deg)) .* ( 0.019993 - 0.000101 .* julian_century) + sin(deg2rad(3 .* geom_mean_anom_sun_deg)) .* 0.000289;
          sun_true_long_deg            = sun_eq_of_ctr + geom_mean_long_sun_deg;
          sun_app_long_deg             = sun_true_long_deg - 0.00569 - 0.00478 .* sin(deg2rad(125.04 - 1934.136 * julian_century));
          mean_obliq_ecliptic_deg      = 23 + (26 + ((21.448 - julian_century .* (46.815 + julian_century .* (0.00059 - julian_century .* 0.001813)))) / 60) / 60;
          obliq_corr_deg               = mean_obliq_ecliptic_deg + 0.00256 .* cos(deg2rad(125.04 - 1934.136 .* julian_century));
          sun_declin_deg               = rad2deg(asin(sin(deg2rad(obliq_corr_deg)) .* sin(deg2rad(sun_app_long_deg))));
          var_y                        = tan(deg2rad(obliq_corr_deg ./ 2)) .* tan(deg2rad(obliq_corr_deg ./ 2));
          eq_of_time_minutes           = 4 * rad2deg(var_y .* sin(2 .* deg2rad(geom_mean_long_sun_deg)) - 2 .* eccent_earth_orbit .* sin(deg2rad(geom_mean_anom_sun_deg)) + 4 .* eccent_earth_orbit .* var_y .* sin(deg2rad(geom_mean_anom_sun_deg)) .* cos(2 .* deg2rad(geom_mean_long_sun_deg)) - 0.5 .* var_y .* var_y .* sin(4 .* deg2rad(geom_mean_long_sun_deg)) - 1.25 .* eccent_earth_orbit .* eccent_earth_orbit .* sin(2 .* deg2rad(geom_mean_anom_sun_deg)));
          true_solar_time_min          = mod(time_of_day .* 1440 + eq_of_time_minutes + 4 .* lon - 60 .* hpeset.par.time_zone,1440);
          hour_angle_deg               = NaN(size(true_solar_time_min));
          idx1                         = true_solar_time_min ./4 < 0;
          idx2                         = true_solar_time_min ./4 >= 0;
          hour_angle_deg(idx1)         = true_solar_time_min(idx1) ./4 + 180;
          hour_angle_deg(idx2)         = true_solar_time_min(idx2) ./4 - 180;
          solar_zenith_angle_deg       = rad2deg(acos(sin(deg2rad(lat)) .* sin(deg2rad(sun_declin_deg)) + cos(deg2rad(lat)) .* cos(deg2rad(sun_declin_deg)) .* cos(deg2rad(hour_angle_deg))));
          solar_elev_angle_deg         = 90 - solar_zenith_angle_deg;

          % calculate atmospheric diffraction dependent on solar elevation angle
          approx_atm_refrac_deg        = NaN(size(solar_elev_angle_deg));
          approx_atm_refrac_deg(solar_elev_angle_deg > 85) = 0;
          idx1                         = solar_elev_angle_deg > 5 & solar_elev_angle_deg <= 85;
          approx_atm_refrac_deg(idx1)  = (58.1 ./ tan(deg2rad(solar_elev_angle_deg(idx1))) - 0.07 ./ (tan(deg2rad(solar_elev_angle_deg(idx1)))).^3 + 0.000086 ./ (tan(deg2rad(solar_elev_angle_deg(idx1)))).^5) ./ 3600; 
          idx2                         = solar_elev_angle_deg > -0.757 & solar_elev_angle_deg <= 5;
          approx_atm_refrac_deg(idx2)  = (1735 + solar_elev_angle_deg(idx2) .* (-518.2 + solar_elev_angle_deg(idx2) .* (103.4 + solar_elev_angle_deg(idx2) .* (-12.79 + solar_elev_angle_deg(idx2) .* 0.711)))) ./ 3600;
          idx3                         = solar_elev_angle_deg <= -0.757;
          approx_atm_refrac_deg(idx3)  = (-20.772 ./ tan(deg2rad(solar_elev_angle_deg(idx3)))) ./ 3600;
          solar_elev_corr_atm_ref_deg  = solar_elev_angle_deg + approx_atm_refrac_deg;

          % calculate solar azimuth angle depending on hour angle
          solar_azimuth_angle          = NaN(size(hour_angle_deg));
          idx1                         = hour_angle_deg > 0;
          solar_azimuth_angle(idx1)    = mod(round((rad2deg(acos(((sin(deg2rad(lat)) .* cos(deg2rad(solar_zenith_angle_deg(idx1)))) - sin(deg2rad(sun_declin_deg(idx1)))) ./ (cos(deg2rad(lat)) .* sin(deg2rad(solar_zenith_angle_deg(idx1)))))) + 180)*100000)/100000,360);
          idx2                         = hour_angle_deg <= 0;
          solar_azimuth_angle(idx2)    = mod(round((540 - rad2deg(acos(((sin(deg2rad(lat)) .* cos(deg2rad(solar_zenith_angle_deg(idx2)))) - sin(deg2rad(sun_declin_deg(idx2)))) ./ (cos(deg2rad(lat)) .* sin(deg2rad(solar_zenith_angle_deg(idx2)))))))*100000)/100000,360);      
          
          % convert solar position to phi/tht
          tht                          = 90 - solar_elev_corr_atm_ref_deg;          % convert to zenith angle (theta) in degree, i.e. zenith is phi = 0; horizon is phi = 90
          phi                          = solar_azimuth_angle;
          phi(phi>=360)                = phi(phi>=360) - 360;              
          phi(phi<0)                   = phi(phi<0)+360;                            % convert to azimuth angle (phi) in degree, so that 0° is North, 90° is East, 180° is South, and 270° is West
          sin_solar_elev               = sin(solar_elev_corr_atm_ref_deg/360*2*pi); % sin of solar elevation angle, where elevation angle = 90 - zenith angle
        
        end
        
        % correct solar azimuth angle for rotation error in orientation of image
        cphi                           = phi + rotation_error;                      % rotate solar position to acount for rotation error in orientation of image
        cphi(cphi>=360)                = cphi(cphi>=360) - 360;              
        cphi(cphi<0)                   = cphi(cphi<0)+360;                 

        %% > add solar track to image and determine raw canopy transmissivity
        if hpeset.sw.scatter
          % use blending of multiple solar disks (true and enlarged size) to simulate the effect of haze / forward scattering of direct radiation
          drad                         = hpeset.par.scatter_drad;
        else
          % use solar disk (true size)
          drad                         = hpeset.par.scatter_drad(1);
        end;
        trans_for                      = zeros(length(loc_time),length(drad));
        
        % add sun track, using simple square to represent sun
        if ~hpeset.sw.detsun
          prad                         = interp1(lens_profile_tht,lens_profile_rpix*image_radius,tht);
          keepix                       = find(tht <= 90);
          x                            = image_center(1) + sin(deg2rad(cphi(keepix))) .* prad(keepix);
          y                            = image_center(2) + cos(deg2rad(cphi(keepix))) .* prad(keepix);
          for six = 1:length(x)
            for dix = 1:length(drad)
              dpix                     = drad(dix)/90*image_radius*sqrt(pi)/2;      % half the side length of a square (in pixels) that has the same area than a circle of radius drad (in degree)
              xmm                      = [max(1,round(x(six)-dpix)):min(size(image_gs,2),round(x(six)+dpix))];
              ymm                      = [max(1,round(y(six)-dpix)):min(size(image_gs,1),round(y(six)+dpix))];          
              switch dix
              case 1
                image_gs(ymm,xmm)      = min(image_gs(ymm,xmm),0.6);
                nolp_1                 = length(find(image_gs(ymm,xmm) == 0.6));
                noap_1                 = numel(image_gs(ymm,xmm));
                trans_for(keepix(six),1)         = nolp_1/noap_1;
              case 2
                image_gs(ymm,xmm)      = min(image_gs(ymm,xmm),0.61);
                nolp_2                 = length(find(image_gs(ymm,xmm) >= 0.60 & image_gs(ymm,xmm) <= 0.61));
                noap_2                 = numel(image_gs(ymm,xmm));
                trans_for(keepix(six),2)         = (nolp_2-nolp_1)/(noap_2-noap_1);
              case 3
                image_gs(ymm,xmm)      = min(image_gs(ymm,xmm),0.62);
                nolp_3                 = length(find(image_gs(ymm,xmm) >= 0.60 & image_gs(ymm,xmm) <= 0.62));
                noap_3                 = numel(image_gs(ymm,xmm));
                trans_for(keepix(six),3)         = (nolp_3-nolp_2)/(noap_3-noap_2);
              end;
            end;
          end;
          
        % add sun track, using realistic lens projection of sun disc
        else
          for tix = 1:length(loc_time)
            for dix = 1:length(drad)
              % outline of sun in direction of tht
              sun.drad                 = drad(dix);
              sun.ctht                 = tht(tix);
              sun.cphi                 = cphi(tix);
              if sun.ctht + sun.drad > 90
                sun.dr                 = 2*diff(interp1([-fliplr(lens_profile_tht(2:end)) lens_profile_tht],([-fliplr(lens_profile_rpix(2:end)) lens_profile_rpix]).*image_radius,sun.ctht+[-sun.drad 0]));
                sun.dt                 = 2*pi*interp1(lens_profile_tht,lens_profile_rpix*image_radius,sun.ctht)/(sin(deg2rad(sun.ctht))*360/(2*sun.drad));
              elseif sun.ctht > 0
                sun.dr                 = diff(interp1([-fliplr(lens_profile_tht(2:end)) lens_profile_tht],([-fliplr(lens_profile_rpix(2:end)) lens_profile_rpix]).*image_radius,sun.ctht+[-sun.drad sun.drad]));
                sun.dt                 = 2*pi*interp1(lens_profile_tht,lens_profile_rpix*image_radius,sun.ctht)/(sin(deg2rad(sun.ctht))*360/(2*sun.drad));
              else
                sun.dr                 = diff(interp1([-fliplr(lens_profile_tht(2:end)) lens_profile_tht],([-fliplr(lens_profile_rpix(2:end)) lens_profile_rpix]).*image_radius,sun.ctht+[-sun.drad sun.drad]));
                sun.dt                 = sun.dr;
              end;
              Frel1                    = [-sqrt((sun.dt/2)^2-(sun.dr/2)^2),interp1(lens_profile_tht,lens_profile_rpix.*image_radius,sun.ctht)]; % position of focal point 1 of ellipe relative to image center (assuming azimuth = 0/North -> at the bottom of thed image)
              Frel2                    = [+sqrt((sun.dt/2)^2-(sun.dr/2)^2),interp1(lens_profile_tht,lens_profile_rpix.*image_radius,sun.ctht)]; % position of focal point 2 of ellipe relative to image center (assuming azimuth = 0/North -> at the bottom of thed image) 
              RMat                     = [cos(deg2rad(sun.cphi)) sin(deg2rad(-sun.cphi)); sin(deg2rad(sun.cphi)) cos(deg2rad(sun.cphi))];       % rotation matrix to position ellipse according to real azimuth
              Fabs1                    = image_center+Frel1*RMat;                                                                               % [xcoor, ycoor] absolute position of focal point 1 of ellipe
              Fabs2                    = image_center+Frel2*RMat;                                                                               % [xcoor, ycoor] absolute position of focal point 2 of ellipe
              switch dix
              case 1
                fix1                   = sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt;                         % identify points within ellipse
                fix2                   = sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt & rgrid <= image_radius; % identify points within ellipse and inside image circle
                image_gs(fix2)         = min(image_gs(fix2),0.6);
                trans_for(tix,1)       = length(find(image_gs(fix2) == 0.6))/numel(image_gs(fix1));
              case 2
                fix3                   = logical([sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt] .* ~fix1);                         % identify points within ellipse
                fix4                   = logical([sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt & rgrid <= image_radius] .* ~fix2); % identify points within ellipse and inside image circle
                image_gs(fix4)         = min(image_gs(fix4),0.61);
                trans_for(tix,2)       = length(find(image_gs(fix4) >= 0.60 & image_gs(fix4) <= 0.61))/numel(image_gs(fix3));
              case 3
                fix5                   = logical([sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt] .* ~fix3 .* ~fix1);                         % identify points within ellipse
                fix6                   = logical([sqrt((xgrid-Fabs1(1)).^2+(ygrid-Fabs1(2)).^2)+sqrt((xgrid-Fabs2(1)).^2+(ygrid-Fabs2(2)).^2) <= sun.dt & rgrid <= image_radius] .* ~fix4 .* ~fix2); % identify points within ellipse and inside image circle
                image_gs(fix6)         = min(image_gs(fix6),0.62);
                trans_for(tix,3)       = length(find(image_gs(fix6) >= 0.60 & image_gs(fix6) <= 0.62))/numel(image_gs(fix5));
              end;
            end;  
          end;
        end;
        trans_for(isnan(trans_for)) = 0;

        %% > calculate transmissivity for direct shortwave radiation
        tweight = hpeset.par.scatter_dwgh;
        switch length(drad)
        case 3
          trans_for_wgt       = (trans_for(:,1)*tweight(1) + trans_for(:,2)*tweight(2) + trans_for(:,3)*tweight(3))/sum(tweight(1:3));
        case 2
          trans_for_wgt       = (trans_for(:,1)*tweight(1) + trans_for(:,2)*tweight(2))/sum(tweight(1:2));
        case 1
          trans_for_wgt       = trans_for(:,1);
        otherwise
          error('Error while calculating transmissivity for direct swr, unknown option');
        end      

        %% > calculate shortwave radiation above canopy
        if isempty(ref) || isempty(ref.swr)
          % this is maximum potential swr, assuming an atmospheric transmissivity of 1 (which also means that all swr is direct / none is diffuse; this can be changed as needed) 
          swr.opn  = max(1367*sin_solar_elev,0);
        else
          % this is swr as measured at the open site, but only if the sun is visisble above the horizon, all other values have been set to nan   
          if imix == 0
            if hpeset.sw.swr < 2
              swr.fdata       = max(swr.data,0);
              swr.fdata(trans_for_wgt < 0.95) = NaN;
              swr.ftime       = swr.time;
            else
              swr.fdata       = interp1(swr.time,max(swr.data,0),loc_time);
              swr.fdata(trans_for_wgt < 0.95) = NaN;
              swr.ftime       = loc_time;
            end;
          else
            if ~isfield(swr,'fdata')
              swr.opn         = interp1(swr.time,max(swr.data,0),loc_time);
            else
              swr.opn         = interp1(swr.ftime,swr.fdata,loc_time);
            end;
          end;
        end;
      end;
      
      if imix > 0 && ~any(isnan(sin_solar_elev))
        
        %% > determine local slope/aspect from DEM
        if ~isempty(hpeset.in.dem)
          [minx,colix]        = min(abs(dem.xvec - loc_x));
          [miny,rowix]        = min(abs(dem.yvec - loc_y));
          if minx > dem.cellsize || miny > dem.cellsize
            loc_slp           = NaN;
            loc_asp           = NaN;
          else
            loc_slp           = dem.slp(rowix,colix);
            loc_asp           = dem.asp(rowix,colix);
          end;
        else
          loc_slp             = NaN;
          loc_asp             = NaN;
        end;
        
        %% > calculate shortwave radiation below canopy
        % split swr into diffuse/direct radiation       
        trans_atm             = swr.opn./max(1367*sin_solar_elev,0);                 % atmospheric transmissivity
        dif_frac              = ones(size(loc_time));                                % fraction of diffuse radiation
        fix                   = sin_solar_elev > 0;
        dif_frac(fix)         = 0.165;
        fix                   = sin_solar_elev > 0 & trans_atm < 0.80;
        dif_frac(fix)         = 0.9511 - 0.1604.*trans_atm(fix) + 4.388.*trans_atm(fix).^2 - 16.638.*trans_atm(fix).^3 + 12.336.*trans_atm(fix).^4;
        fix                   = sin_solar_elev > 0 & trans_atm < 0.22;
        dif_frac(fix)         = 1 - 0.09.*trans_atm(fix);
        dir_frac              = 1-dif_frac;
        
        % determine maximum direct swr (if sensor is pointed towards sun)
        switch hpeset.par.swr_maxdir
        case 0
          max_dir             = 1367.*ones(trans_atm);                               % using solar constant
        case 1
          max_dir             = 1367.*trans_atm.*dir_frac;                           % using measured swr transmissivity and modelled split into direct and diffuse component
        case 2
          trans_dirmax        = 1 - 0.165;                                           % value copied from from above fsm splitting function
          max_dir             = 1367.*exp(log(trans_dirmax)./max(sin_solar_elev,0)); % using de beer's law to acount for length of pathway in terms of zenith angle. The extinction coefficient is based on the maximum transmission of direct swr if sun is at zenith
        case 3
          trans_dirmax        = 1 - 0.165;                                        
          max_dir_1           = 1367.*exp(log(trans_dirmax)./max(sin_solar_elev,0));            
          max_dir_2           = 1367.*trans_atm.*dir_frac;                          
          max_dir             = min(max_dir_1,max_dir_2);                           % applying both above criteria 1 and 2
        otherwise
          error('Error while calculating maximum direct swr, unknown option');
        end;

        % calculate diffuse swr below canopy, neglecting the distinction of / possible differences between SVF_flat, SVF_tsun, and SVF_incl as the latter two are unavailable (unless you are using synthetic HPs)
        if ~isempty(ref)
          swr.for_dif         = dif_frac.*swr.opn./ref.SVF_flat.*img.SVF_flat(imix);
        else
          swr.for_dif         = dif_frac.*swr.opn.*img.SVF_flat(imix);
        end;
        
        % calculate direct swr below canopy
        cang_1                = max(sin((90-tht)/360*2*pi),0.001);          % cosine of angle between flat surface and position of sun
        cang_2                = max(sin((90-tht)/360*2*pi).*cos(loc_slp/360*2*pi) + cos((90-tht)/360*2*pi).*sin(loc_slp/360*2*pi).*cos((loc_asp-cphi)/360*2*pi),0); % cosine of angle between inclined surface and position of sun
        nanix                 = ones(size(swr.opn));
        nanix(isnan(swr.opn)) = NaN;                                       % this nan mask is necessary with old Matlab versions that ignored NaNs when calculating min or max
        swr.for_dir_flat      = min(dir_frac.*swr.opn        ,max_dir).*trans_for_wgt.*nanix;
        swr.for_dir_tsun      = min(dir_frac.*swr.opn./cang_1,max_dir).*trans_for_wgt.*nanix;
        swr.for_dir_incl      = min(dir_frac.*swr.opn./cang_1,max_dir).*trans_for_wgt.*cang_2.*nanix;
        
        % calculate total swr below canopy
        swr.for_all_flat      = swr.for_dif + swr.for_dir_flat;             % swr for a sensor pointing to the zenith
        swr.for_all_tsun      = swr.for_dif + swr.for_dir_tsun;             % swr for a sensor pointing towards the sun (might be useful to take the perspectiv of a plant)
        swr.for_all_incl      = swr.for_dif + swr.for_dir_incl;             % swr for a sensor normal to local aspect/slope (= aspect/slope of nearest DEM grid cell)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% RESULT OUTPUT
        
        %% > load validation data
        if ~isempty(val_data)
          % load data (to be provided in ascii format)
          try
            fid               = fopen(fullfile(hpeset.in.swr,val_data),'r');
            formatstr         = '%s%f';
            fgetl(fid);
            data              = textscan(fid,formatstr,'delimiter','\t');
            try
              swr.val_time    = datenum(data{1},'dd.mm.yyyy HH:MM:SS');
            catch
              swr.val_time    = datenum(data{1},'dd.mm.yyyy HH:MM');
            end;
            swr.val_data      = data{2}*swr_scaling;
            fclose(fid);
          catch
            try
              fclose(fid);
            end;
            error('Error while reading swr validation data file');
          end;
        else
          swr.val_time        = [];
          swr.val_data        = [];
        end;
        
        %% > plot swr data and save figure (time series)
        if ~isempty(hpeset.for.swr)
          % plot swr data
          fh = figure('visible','off');
          hold on;
          plot(loc_time,swr.opn,'b');
          plot(loc_time,swr.for_all_flat,'g');
          plot(loc_time,swr.for_dif,'k');
          plot(loc_time,swr.for_dir_flat,'y');
          if ~isempty(swr.val_data)
            plot(swr.val_time,swr.val_data,'r');
            legend({'swr open','swr for mod','swr dif for','swr dir for','swr for obs'})
          else
            legend({'swr open','swr for','swr dif for','swr dir for'})
          end;
          set(gca,'xlim',loc_time([1,end]));
          datetick('x','HH:MM','keeplimits')
          box on;
          xlabel('daytime')
          ylabel('SWR [W/m2]')
          tlabel              = ['SWR.' sprintf('%6.0f',loc_x) '.' sprintf('%6.0f',loc_y)];
          title(tlabel);

          % save figure
          if imix == 0
            [opath,ofile,~]   = fileparts(fullfile(hpeset.out.fig,ref.filename));
          else
            [opath,ofile,~]   = fileparts(fullfile(hpeset.out.fig,img.filename{imix}));
          end;
          for fix = 1:length(hpeset.for.swr)
            if strcmpi(hpeset.for.swr{fix},'fig')
              set(fh,'visible','on');
            else
              set(fh,'visible','off');
            end;
            saveas(fh,fullfile(opath,strcat(ofile,'_swr_ts')),hpeset.for.swr{fix});
          end;
          close(fh);
          drawnow;
        end;
        
        %% > plot swr data and save figure (scatter plot)
        if ~isempty(hpeset.for.swr) && ~isempty(swr.val_data) 
          % plot swr data
          fh = figure('visible','off');
          hold on;
          data_10min = nan(length(loc_time(1):1/24/6:loc_time(end)),2);
          for tix = 1:length(data_10min)
            fix = swr.val_time >= loc_time(1)+(tix-1)*(1/24/6) & swr.val_time < loc_time(1)+(tix)*(1/24/6);
            data_10min(tix,1) = mean(swr.val_data(fix));
            fix = loc_time >= loc_time(1)+(tix-1)*(1/24/6) & loc_time < loc_time(1)+(tix)*(1/24/6);
            data_10min(tix,2) = mean(swr.for_all_flat(fix));
          end;
          data_60min = nan(length(loc_time(1):1/24:loc_time(end)),2);
          for tix = 1:length(data_60min)
            fix = swr.val_time >= loc_time(1)+(tix-1)*(1/24) & swr.val_time < loc_time(1)+(tix)*(1/24);
            data_60min(tix,1) = mean(swr.val_data(fix));
            fix = loc_time >= loc_time(1)+(tix-1)*(1/24) & loc_time < loc_time(1)+(tix)*(1/24);
            data_60min(tix,2) = mean(swr.for_all_flat(fix));
          end;
          plot(data_10min(:,1),data_10min(:,2),'marker','o','markersize',2,'markerfacecolor',[0 0 1],'markeredgecolor',[0 0 1],'linestyle','none');
          plot(data_60min(:,1),data_60min(:,2),'marker','o','markersize',2,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0],'linestyle','none');
          xlim = get(gca,'xlim');
          ylim = get(gca,'ylim');
          tlim = [0,max(xlim(2),ylim(2))];
          plot(tlim,tlim,'k:');
          set(gca,'xlim',tlim,'ylim',tlim);
          legend({'10-minute data','1-hour data'},'Location','SouthEast')
          box on;
          xlabel('observed SWR [W/m2]')
          ylabel('modelled SWR [W/m2]')
          tlabel = ['SWR.' sprintf('%6.0f',loc_x) '.' sprintf('%6.0f',loc_y)];
          title(tlabel);
        
          % save figure
          if imix == 0
            [opath,ofile,~]    = fileparts(fullfile(hpeset.out.fig,ref.filename));
          else
            [opath,ofile,~]    = fileparts(fullfile(hpeset.out.fig,img.filename{imix}));
          end;
          for fix = 1:length(hpeset.for.swr)
            if strcmpi(hpeset.for.swr{fix},'fig')
              set(fh,'visible','on');
            else
              set(fh,'visible','off');
            end;
            saveas(fh,fullfile(opath,strcat(ofile,'_swr_sc')),hpeset.for.swr{fix});
          end;
          close(fh);
          drawnow;
        end;
        
        %% > assemble and save swr data
        if ~isempty(hpeset.for.dat)
          [opath,ofile,~] = fileparts(fullfile(hpeset.out.dat,img.filename{imix}));
          for fix = 1:length(hpeset.for.dat)
            switch lower(hpeset.for.dat{fix})
            case 'txt' % write swr data to text file
              try
                fid = fopen(fullfile(opath,strcat(ofile,'_swr','.txt')),'w');
                fprintf(fid,'%s\t','TIMESTAMP');            
                fprintf(fid,'%s\t','SWR_OPN_TOT');            
                fprintf(fid,'%s\t','SWR_OPN_DIR');            
                fprintf(fid,'%s\t','SWR_FOR_FLAT_TOT');           
                fprintf(fid,'%s\t','SWR_FOR_FLAT_DIR');
                fprintf(fid,'%s\t','SWR_FOR_TSUN_TOT');           
                fprintf(fid,'%s\t','SWR_FOR_TSUN_DIR');
                fprintf(fid,'%s\t','SWR_FOR_INCL_TOT');           
                fprintf(fid,'%s\t','SWR_FOR_INCL_DIR');
                fprintf(fid,'%s\n','SWR_FOR_TAU');    
                for tix = 1:length(loc_time)
                  fprintf(fid,'%s\t',datestr(loc_time(tix),'dd.mm.yyyy HH:MM:SS'));
                  fprintf(fid,'%0.1f\t%0.1f\t',[swr.opn(tix),dir_frac(tix)*swr.opn(tix)]);
                  fprintf(fid,'%0.1f\t%0.1f\t',[swr.for_all_flat(tix),swr.for_dir_flat(tix)]);
                  fprintf(fid,'%0.1f\t%0.1f\t',[swr.for_all_tsun(tix),swr.for_dir_tsun(tix)]);
                  fprintf(fid,'%0.1f\t%0.1f\t',[swr.for_all_incl(tix),swr.for_dir_incl(tix)]);
                  fprintf(fid,'%0.2f\n',trans_for_wgt(tix));
                end;
                fclose(fid);
              catch
                try
                  fclose(fid);
                end;
                error('Error while saving swr data to file');
              end;
            case 'mat' % write swr data to matlab struct array
              HPE_SWR.timestamp        = loc_time;
              HPE_SWR.swr_opn_tot      = swr.opn;
              HPE_SWR.swr_opn_dir      = dir_frac.*swr.opn;
              HPE_SWR.swr_for_flat_tot = swr.for_all_flat; 
              HPE_SWR.swr_for_flat_dir = swr.for_dir_flat; 
              HPE_SWR.swr_for_tsun_tot = swr.for_all_tsun; 
              HPE_SWR.swr_for_tsun_dir = swr.for_dir_tsun; 
              HPE_SWR.swr_for_incl_tot = swr.for_all_incl; 
              HPE_SWR.swr_for_incl_dir = swr.for_dir_incl; 
              HPE_SWR.swr_for_tau      = trans_for_wgt; 
              save(fullfile(opath,strcat(ofile,'_swr','.mat')),'HPE_SWR');
            otherwise
              error('Error while saving swr data, unknown option');
            end;            
          end;
        end;
      end;
      
      if ~isempty(hpeset.for.img)
        %% > reclassify, plot and save hemispherical image 
        % reclassify image
        fh                   = figure('visible','off');
        image_r              = ones(size(image_gs));
        image_g              = ones(size(image_gs));
        image_b              = ones(size(image_gs));
        % dark green (canopy)
        fix                  = image_gs == 0;
        image_r(fix)         = 0;
        image_g(fix)         = 0.2;
        image_b(fix)         = 0;
        % dark gray (terrain)
        fix                  = image_gs == 0.2;
        image_r(fix)         = 0.4;
        image_g(fix)         = 0.4;
        image_b(fix)         = 0.4;
        % dark brown (mask)
        fix                  = image_gs == 0.3;
        image_r(fix)         = 0.4;
        image_g(fix)         = 0.3;
        image_b(fix)         = 0.2;
        % white (mask)
        fix                  = image_gs == 0.5;
        image_r(fix)         = 1;
        image_g(fix)         = 1;
        image_b(fix)         = 1;
        if hpeset.sw.suntrack
          % yellow (sun track r1)
          fix                = image_gs == 0.6;
          image_r(fix)       = 0.94;
          image_g(fix)       = 1;
          image_b(fix)       = 0.06;
          % light orange (sun track r2)
          fix                = image_gs == 0.61;
          image_r(fix)       = 0.94;
          image_g(fix)       = 0.8;
          image_b(fix)       = 0.06;
          % dark orange (sun track r3)
          fix                = image_gs == 0.62;
          image_r(fix)       = 0.94;
          image_g(fix)       = 0.6;
          image_b(fix)       = 0.06;
        else
          % light blue (sun track r1)
          fix                = image_gs == 0.6;
          image_r(fix)       = 0.1;
          image_g(fix)       = 0.5;
          image_b(fix)       = 1;
          % light blue (sun track r2)
          fix                = image_gs == 0.61;
          image_r(fix)       = 0.1;
          image_g(fix)       = 0.5;
          image_b(fix)       = 1;
          % light blue (sun track r3)
          fix                = image_gs == 0.62;
          image_r(fix)       = 0.1;
          image_g(fix)       = 0.5;
          image_b(fix)       = 1;
        end;
        % light blue (sky)
        fix                  = image_gs == 1;
        image_r(fix)         = 0.1;
        image_g(fix)         = 0.5;
        image_b(fix)         = 1;

        % plot image
        image_c              = [];
        image_c(:,:,1)       = image_r;
        image_c(:,:,2)       = image_g;
        image_c(:,:,3)       = image_b;
        image(image_c)
        set(gca,'xlim',image_center(1)+[-1 1].*image_radius,'ylim',image_center(2)+[-1 1].*image_radius);
        title(['RHI' sprintf('%6.0f',loc_x) '.' sprintf('%6.0f',loc_y)]);
        set(gcf,'position',[100 100,640,600])

        % crop image
        image_cr             = image_c(max(image_center(2)-image_radius,1):min(image_center(2)+image_radius,size(image_c,1)),...
                                       max(image_center(1)-image_radius,1):min(image_center(1)+image_radius,size(image_c,2)),:);

        % save image
        if imix == 0
          [opath,ofile,~]    = fileparts(fullfile(hpeset.out.img,ref.filename));
        else
          [opath,ofile,~]    = fileparts(fullfile(hpeset.out.img,img.filename{imix}));
        end;
        for fix = 1:length(hpeset.for.img)
          imwrite(image_cr,fullfile(opath,strcat(ofile,'_thr.',hpeset.for.img{fix})));
        end;
        close(fh);
        drawnow;
      end;
    end;
    
    % progress info
    if hpeset.sw.report
      disp(['   ' datestr(now,'HH:MM:SS') ' > ' num2str(imix-stix+1) '/' num2str(length(img.filename)-stix+1) ' images processed'])
      drawnow;
    end;
    
  end;
  
  %% > save SVF calculations
  if hpeset.sw.svf
    try
      fid = fopen(fullfile(hpeset.out.dat,'HPE_SVF.txt'),'w');
      fprintf(fid,'%s\t','Image_file');            
      fprintf(fid,'%s\t','Easting');            
      fprintf(fid,'%s\t','Northing');            
      fprintf(fid,'%s\t','SVF_hemi');            
      fprintf(fid,'%s\t','SVF_flat');            
      fprintf(fid,'%s\t','SVF_flat_south'); 
      fprintf(fid,'%s\n','SVF_flat_north');            
      for imix = 1:img.no_ims
        fprintf(fid,'%s\t',img.filename{imix});
        fprintf(fid,'%0.0f\t%0.0f\t',[loc_x,loc_y]);
        fprintf(fid,'%0.2f\t%0.2f\t%0.2f\t%0.2f\n',[img.SVF_hemi(imix),img.SVF_flat(imix),img.SVF_flat_S(imix),img.SVF_flat_N(imix)]);
      end
      fclose(fid);
    catch
      try
        fclose(fid);
      end;
      error('Error while saving svf data to file');
    end;
  end;
  
  if hpeset.sw.report
    disp('   -------------------------------------------------------')
    drawnow;
  end;
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EMBEDDED AUXILIARY FUNCTIONS
%% > deg2rad
function angleInRadians = deg2rad(angleInDegrees)
  angleInRadians   = (pi/180) * angleInDegrees;
end

%% > rad2deg
function angleInDegrees = rad2deg(angleInRadians)
  angleInDegrees   = (180/pi) * angleInRadians;
end

%% > versioncompare
function vm = versioncompare(v1,v2)
  % returns +1 if v1 > v2; returns -1 if v1 < v2; returns 0 if v1 = v2 
  try
    vers1          = [0 0 0];
    vstr1          = v1;
    xix            = strfind(vstr1,'.');
    switch length(xix)
    case 0
      vers1(1)     = str2num(vstr1);
    case 1
      vers1(1)     = str2num(vstr1(1:xix(1)-1));
      vers1(2)     = str2num(vstr1(xix(1)+1:end));
    case 2
      vers1(1)     = str2num(vstr1(1:xix(1)-1));
      vers1(2)     = str2num(vstr1(xix(1)+1:xix(2)-1));
      vers1(3)     = str2num(vstr1(xix(2)+1:end));
    otherwise
      error(' ');
    end;
    vers2          = [0 0 0];
    vstr2          = v2;
    xix            = strfind(vstr2,'.');
    switch length(xix)
    case 0
      vers2(1)     = str2num(vstr2);
    case 1
      vers2(1)     = str2num(vstr2(1:xix(1)-1));
      vers2(2)     = str2num(vstr2(xix(1)+1:end));
    case 2
      vers2(1)     = str2num(vstr2(1:xix(1)-1));
      vers2(2)     = str2num(vstr2(xix(1)+1:xix(2)-1));
      vers2(3)     = str2num(vstr2(xix(2)+1:end));
    otherwise
      error(' ');
    end;
    if vers1(1) > vers2(1)
      vm = +1;
    elseif vers1(1) < vers2(1)
      vm = -1;
    elseif vers1(2) > vers2(2)
      vm = +1;
    elseif vers1(2) < vers2(2)
      vm = -1;
    elseif vers1(3) > vers2(3)
      vm = +1;
    elseif vers1(3) < vers2(3)
      vm = -1;
    else
      vm = 0;
    end;
  catch
    vm   = NaN;
  end;
end

%% > utm2deg
function  [Lat,Lon] = utm2deg(xx,yy,zone,hemi)
  % Credit: based on utm2lonlat found on Matlab file exchange; accessed on 2018/08/25
  % Copyright (c) 2013, Erwin N. All rights reserved. 
  x      = xx(:);
  y      = yy(:);
  sa     = 6378137.000000; sb = 6356752.314245;
  e2     = ( ( ( sa .^ 2 ) - ( sb .^ 2 ) ) .^ 0.5 ) ./ sb;
  e2squared = e2 .^ 2;
  c      = ( sa .^ 2 ) ./ sb;
  switch upper(hemi)
  case 'N'
    X    = x - 500000;
    Y    = y;
  case 'S'
    X    = x - 500000;
    Y    = y - 10000000;
  end;
  S      = ( ( zone .* 6 ) - 183 ); 
  lat    =  Y ./ ( 6366197.724 .* 0.9996 );                                    
  v      = ( c ./ ( ( 1 + ( e2squared .* ( cos(lat) ) .^ 2 ) ) ) .^ 0.5 ) .* 0.9996;
  a      = X ./ v;
  a1     = sin( 2 .* lat );
  a2     = a1 .* ( cos(lat) ) .^ 2;
  j2     = lat + ( a1 ./ 2 );
  j4     = ( ( 3 .* j2 ) + a2 ) ./ 4;
  j6     = ( ( 5 .* j4 ) + ( a2 .* ( cos(lat) ) .^ 2) ) ./ 3;
  alpha  = ( 3 ./ 4 ) .* e2squared;
  beta   = ( 5 ./ 3 ) .* alpha .^ 2;
  gamma  = ( 35 ./ 27 ) .* alpha .^ 3;
  Bm     = 0.9996 .* c .* ( lat - alpha .* j2 + beta .* j4 - gamma .* j6 );
  b      = ( Y - Bm ) ./ v;
  Epsi   = ( ( e2squared .* a.^2 ) ./ 2 ) .* ( cos(lat) ).^ 2;
  Eps    = a .* ( 1 - ( Epsi ./ 3 ) );
  nab    = ( b .* ( 1 - Epsi ) ) + lat;
  senoheps = ( exp(Eps) - exp(-Eps) ) ./ 2;
  Delt   = atan(senoheps ./ (cos(nab) ) );
  TaO    = atan(cos(Delt) .* tan(nab));
  longitude = (Delt .* (180/pi) ) + S;
  latitude = ( lat + ( 1 + e2squared .* (cos(lat).^2) - ( 3/2 ) ...
    .* e2squared .* sin(lat) .* cos(lat) .* ( TaO - lat ) ) ...
    .* ( TaO - lat ) ) .* (180/pi);
  Lat    = latitude;
  Lon    = longitude;
end

%% > load_ascii_grid
function answer = load_ascii_grid(filepath)
  answer = [];
  [path,name,ext] = fileparts(filepath);
  fid = fopen(filepath,'r');
  if fid == -1
    answer = 'File inaccessible';
    return;
  end;
  try % reading grid according to ASCII GIS format
    answer = 'Error reading grid.ncols';  
    grid.ncols = fgets(fid);
    fix = strfind(lower(grid.ncols),'ncols');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.ncols(1:fix-1) grid.ncols(fix+5:end)];
    grid.ncols = str2num(hlpstr);
    answer = 'Error reading grid.nrows';  
    grid.nrows = fgets(fid);
    fix = strfind(lower(grid.nrows),'nrows');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.nrows(1:fix-1) grid.nrows(fix+5:end)];
    grid.nrows = str2num(hlpstr);
    answer = 'Error reading grid.xllcorner';  
    grid.xllcorner = fgets(fid);
    fix = strfind(lower(grid.xllcorner),'xllcorner');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.xllcorner(1:fix-1) grid.xllcorner(fix+9:end)];
    grid.xllcorner = str2num(hlpstr);
    answer = 'Error reading grid.yllcorner';  
    grid.yllcorner = fgets(fid);
    fix = strfind(lower(grid.yllcorner),'yllcorner');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.yllcorner(1:fix-1) grid.yllcorner(fix+9:end)];
    grid.yllcorner = str2num(hlpstr);
    answer = 'Error reading grid.cellsize';  
    grid.cellsize = fgets(fid);
    fix = strfind(lower(grid.cellsize),'cellsize');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.cellsize(1:fix-1) grid.cellsize(fix+8:end)];
    grid.cellsize = str2num(hlpstr);
    answer = 'Error reading grid.NODATA_value';  
    grid.NODATA_value = fgets(fid);
    fix = strfind(lower(grid.NODATA_value),'nodata_value');
    if isempty(fix)
      error(' ');
    end;
    hlpstr = [grid.NODATA_value(1:fix-1) grid.NODATA_value(fix+12:end)];
    grid.NODATA_value = str2num(hlpstr);
    answer = 'Error reading grid.data';  
    formatstr = '';
    for cix = 1:grid.ncols
      formatstr = [formatstr '%f'];
    end;
    data = textscan(fid,formatstr,grid.nrows);
    for cix = 1:grid.ncols
      grid.data(:,cix) = data{cix};
    end;
    grid.data = flipud(grid.data); % conversion necessary in order to comply with own asciigrid standards
    clear data;
    fclose(fid);
    answer = grid;
  catch
    try
      fclose(fid);
    end;
    answer = 'File is not a standard ASCII grid';  
  end;
end