run.name    = 'aggregates';
run.mode    = 0;                % 0: Deterministic
                                % 1: Probabilistic
run.mass    = 0;                % 0: Mass if specified
                                % 1: Mass computed from height, wind and duration
                                % 2: Mass computed from thermal theory - to add
                                
grid.xLim   = [-50000 200000];  % X limits (UTM, m)
grid.yLim   = [-50000 50000];   % Y limits (UTM, m)
grid.res    = [1000];             % Grid resolution (m)
grid.elev   = [0];                % Mean elevation
grid.zone   = [];               % Zone

vent.x      = 0;
vent.y      = 0;
vent.z      = 1;

wind.pth    = [];               % If not empty, then use it
wind.dir    = 90;               % Degree from North, wind direction
wind.trop   = 12000;            % Tropopause height (m asl)
wind.strat  = 15000;            % Stratosphere height (m asl)
wind.vel    = [20];     % Velocity at tropopause (m/s)
wind.model  = 1;                % 1: WM1 of Bonadonna and Phillips
                                % 2: WM2 of Bonadonna and Phillips

tgsd.pth    = [];
tgsd.med    = -1;
tgsd.std    = 2;
tgsd.agg    = [.1,.25,.5,.75,.9];
tgsd.minD   = 5;


conf.plumeHeight   = 25000;
conf.mass          = 1e12;
conf.duration      = 6; % Hours
conf.eddyConst     = 0.04;
conf.diffCoef      = 112;
conf.ftt           = 1563;
conf.lithDens      = 2600;
conf.pumDens       = 513;
conf.colStep       = 50;
conf.partStep      = 50;
conf.plumeAlpha    = 3;
conf.plumeBeta     = 2.5;


%% 
data.run    = run;
data.grid   = grid;
data.vent   = vent;
data.wind   = wind;
data.tgsd   = tgsd;
data.conf   = conf;

data = runT2(data);
