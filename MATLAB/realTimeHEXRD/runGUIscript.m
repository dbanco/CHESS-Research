%% Run gui with specified input/output dirs
filePath = '/nfs/chess/raw/2023-2/id3a/shanks-3731-a/ti-2-tension/';
center = [3026,1958];
omInd = 5;
omInd2 = 44;
startFnum = 113;
%% Ring 1 processing (outermost of these 3)
outPath = '/nfs/chess/user/dbanco/outputs/ti-2-tension_ring1b/';
r = round(610*tan(pi*7.61/180)/0.0748 - 10);
r1 = r-20;
r2 = r+20;
hexd_realTime_gui(filePath,outPath,r,r1,r2,center,omInd,omInd2,startFnum);

%% Ring 2 processing
outPath = '/nfs/chess/user/dbanco/outputs/ti-2-tension_ring2b/';
r = round(610*tan(pi*7.61/180)/0.0748 - 55);
r1 = r-12;
r2 = r+12;
hexd_realTime_gui(filePath,outPath,r,r1,r2,center,omInd,omInd2,startFnum);
    
%% Ring 3 processing (innermost)
outPath = '/nfs/chess/user/dbanco/outputs/ti-2-tension_ring3b/';
r = round(610*tan(pi*7.61/180)/0.0748 - 140);
r1 = r-15;
r2 = r+15;
hexd_realTime_gui(filePath,outPath,r,r1,r2,center,omInd,omInd2,startFnum);