function OUT = macro_sim_core_hazard_SSA_FINAL(cfg)
% =========================================================================
% macro_sim_core_hazard_SSA_FINAL
%
% Exact event-driven direct-method SSA for the updated macro model:
%
%       V --(r*q_i)--------> V + V
%       V --(lambda_F,i)--> F
%       F --(mu)----------> empty
%
% States: 0=empty, 1=viable V, 2=ferroptotically committed F.
% q_i is the fraction of six von-Neumann neighbour sites that are empty.
%
% This solver has NO numerical time step. Between events the state and all
% propensities are constant, so the direct-method SSA is exact for the
% specified continuous-time Markov jump process.
%
% REQUIRED cfg:
%   L, R0, Tend_h, record_dt_h, r_h, mu_h, seed, lambda_mode
%
% lambda_mode:
%   'scalar'    -> lambdaF_h
%   'site_map'  -> lambdaF_site_h (LxLxL)
%   'cell_mark' -> persistent cell marks; initial map or lambdaF_h fallback
%
% Optional:
%   initial_center
%   snapshot_times_h
%   return_final_state (default true)
%   compute_cavity     (default true)
%   store_event_log    (default false)
%   verbose            (default false)
%
% For cell_mark daughters:
%   daughter_lambda_rule = 'inherit' (default) | 'sample_pool' | 'site_map'
%   sample_pool requires lambdaF_pool_h
%   site_map requires lambdaF_site_h
%
% IMPORTANT: mesoscopic FPT samples are NOT used as macro lifetime deadlines.
% The mesoscopic 0--60 min statistics have already been coarse-grained into
% lambda_F, which is used here as the continuous-time V->F hazard.
% =========================================================================

cfg = validateCfg(cfg);
rng(cfg.seed,'twister');

L = cfg.L;
nSite = L^3;

if isempty(cfg.initial_center)
    ctr = [(L+1)/2, (L+1)/2, (L+1)/2];
else
    ctr = double(cfg.initial_center(:).');
end

%% Initial state
[xg,yg,zg] = ndgrid(1:L,1:L,1:L);
radiusSqAll = (double(xg)-ctr(1)).^2 + ...
              (double(yg)-ctr(2)).^2 + ...
              (double(zg)-ctr(3)).^2;

initialMask = radiusSqAll <= cfg.R0^2;

state = zeros(nSite,1,'uint8');
initialIdx = find(initialMask(:));
state(initialIdx) = uint8(1);

lambdaCell = nan(nSite,1);

if strcmp(cfg.lambda_mode,'cell_mark')
    if ~isempty(cfg.lambdaF_initial_map_h)
        A = double(cfg.lambdaF_initial_map_h(:));
        lambdaCell(initialIdx) = A(initialIdx);
    else
        lambdaCell(initialIdx) = cfg.lambdaF_h;
    end
end

if strcmp(cfg.lambda_mode,'site_map') || ...
        (strcmp(cfg.lambda_mode,'cell_mark') && ...
         strcmp(cfg.daughter_lambda_rule,'site_map'))
    lambdaSite = double(cfg.lambdaF_site_h(:));
else
    lambdaSite = [];
end

%% Active sets
Vlist = zeros(nSite,1,'uint32');
Flist = zeros(nSite,1,'uint32');
posV = zeros(nSite,1,'int32');
posF = zeros(nSite,1,'int32');

nV = numel(initialIdx);
nF = 0;

Vlist(1:nV) = uint32(initialIdx);
for k=1:nV
    posV(initialIdx(k)) = int32(k);
end

%% Empty-neighbour counts
K = zeros(3,3,3);
K(2,2,1)=1; K(2,2,3)=1;
K(2,1,2)=1; K(2,3,2)=1;
K(1,2,2)=1; K(3,2,2)=1;

empty3 = reshape(state==0,[L L L]);
nEmptyCount = uint8(convn(double(empty3),K,'same'));
nEmptyCount = nEmptyCount(:);

%% Recording arrays
tRec = (0:cfg.record_dt_h:cfg.Tend_h).';
if isempty(tRec) || abs(tRec(end)-cfg.Tend_h)>1e-12
    tRec = [tRec; cfg.Tend_h];
end
tRec = unique(tRec,'stable');
nRec = numel(tRec);

NV = zeros(nRec,1); NF = zeros(nRec,1); NT = zeros(nRec,1);
Rmax_sites = zeros(nRec,1);
Rboundary_mean_sites = zeros(nRec,1);
roughness_sites = zeros(nRec,1);
cavity_fraction = zeros(nRec,1);
edgeOcc = false(nRec,1);

cumBirth = zeros(nRec,1);
cumFerro = zeros(nRec,1);
cumClear = zeros(nRec,1);

snapTimes = double(cfg.snapshot_times_h(:));
snapTimes = snapTimes(snapTimes>=0 & snapTimes<=cfg.Tend_h);
snapTimes = unique(snapTimes,'stable');
snapshots = repmat(struct('time_h',NaN,'state',[]),numel(snapTimes),1);

AUTC_NV_exact = 0;
AUTC_NT_exact = 0;

nBirthTotal = 0;
nFerroTotal = 0;
nClearTotal = 0;

viableExtinctionTime_h = Inf;
totalClearanceTime_h = Inf;

if cfg.store_event_log
    eventTime = zeros(0,1);
    eventType = zeros(0,1,'uint8'); % 1 birth, 2 ferro, 3 clear
    eventSite = zeros(0,1,'uint32');
    eventTarget = zeros(0,1,'uint32');
else
    eventTime = [];
    eventType = [];
    eventSite = [];
    eventTarget = [];
end

M = morphologyMetricsFromLists( ...
    state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);

NV(1)=M.NV; NF(1)=M.NF; NT(1)=M.NT;
Rmax_sites(1)=M.Rmax;
Rboundary_mean_sites(1)=M.RboundaryMean;
roughness_sites(1)=M.roughness;
cavity_fraction(1)=M.cavity;
edgeOcc(1)=M.edgeTouch;

nextRec = 2;
nextSnap = 1;

while nextSnap<=numel(snapTimes) && abs(snapTimes(nextSnap))<1e-12
    snapshots(nextSnap).time_h = 0;
    snapshots(nextSnap).state = reshape(state,[L L L]);
    nextSnap = nextSnap+1;
end

%% Exact event-driven SSA
t = 0;
nEventTotal = 0;

while t < cfg.Tend_h-1e-12

    % Current propensities
    if nV>0
        vidx = double(Vlist(1:nV));
        nEmptyV = double(nEmptyCount(vidx));

        totalDirectedVEmptyEdges = sum(nEmptyV);
        A_birth = cfg.r_h * totalDirectedVEmptyEdges / 6;

        switch cfg.lambda_mode
            case 'scalar'
                A_ferro = cfg.lambdaF_h * nV;
                ferroWeights = [];
            case 'site_map'
                ferroWeights = lambdaSite(vidx);
                A_ferro = sum(ferroWeights);
            case 'cell_mark'
                ferroWeights = lambdaCell(vidx);
                assert(all(isfinite(ferroWeights)) && all(ferroWeights>=0), ...
                    'Invalid viable-cell lambda marks encountered.');
                A_ferro = sum(ferroWeights);
            otherwise
                error('Unknown lambda_mode.');
        end
    else
        vidx = [];
        nEmptyV = [];
        totalDirectedVEmptyEdges = 0;
        A_birth = 0;
        A_ferro = 0;
        ferroWeights = [];
    end

    A_clear = cfg.mu_h * nF;
    A_total = A_birth + A_ferro + A_clear;

    % No events remain: freeze state to Tend.
    if A_total<=0
        dtRemain = cfg.Tend_h-t;
        AUTC_NV_exact = AUTC_NV_exact + nV*dtRemain;
        AUTC_NT_exact = AUTC_NT_exact + (nV+nF)*dtRemain;

        while nextSnap<=numel(snapTimes)
            snapshots(nextSnap).time_h = snapTimes(nextSnap);
            snapshots(nextSnap).state = reshape(state,[L L L]);
            nextSnap=nextSnap+1;
        end

        while nextRec<=nRec
            M = morphologyMetricsFromLists( ...
                state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);
            NV(nextRec)=M.NV; NF(nextRec)=M.NF; NT(nextRec)=M.NT;
            Rmax_sites(nextRec)=M.Rmax;
            Rboundary_mean_sites(nextRec)=M.RboundaryMean;
            roughness_sites(nextRec)=M.roughness;
            cavity_fraction(nextRec)=M.cavity;
            edgeOcc(nextRec)=M.edgeTouch;
            cumBirth(nextRec)=nBirthTotal;
            cumFerro(nextRec)=nFerroTotal;
            cumClear(nextRec)=nClearTotal;
            nextRec=nextRec+1;
        end

        t=cfg.Tend_h;
        break;
    end

    % Exact next-event time
    dtEvent = -log(rand)/A_total;
    tEvent = t+dtEvent;

    if tEvent>cfg.Tend_h
        dtRemain = cfg.Tend_h-t;
        AUTC_NV_exact = AUTC_NV_exact + nV*dtRemain;
        AUTC_NT_exact = AUTC_NT_exact + (nV+nF)*dtRemain;

        while nextSnap<=numel(snapTimes)
            snapshots(nextSnap).time_h = snapTimes(nextSnap);
            snapshots(nextSnap).state = reshape(state,[L L L]);
            nextSnap=nextSnap+1;
        end

        while nextRec<=nRec
            M = morphologyMetricsFromLists( ...
                state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);
            NV(nextRec)=M.NV; NF(nextRec)=M.NF; NT(nextRec)=M.NT;
            Rmax_sites(nextRec)=M.Rmax;
            Rboundary_mean_sites(nextRec)=M.RboundaryMean;
            roughness_sites(nextRec)=M.roughness;
            cavity_fraction(nextRec)=M.cavity;
            edgeOcc(nextRec)=M.edgeTouch;
            cumBirth(nextRec)=nBirthTotal;
            cumFerro(nextRec)=nFerroTotal;
            cumClear(nextRec)=nClearTotal;
            nextRec=nextRec+1;
        end

        t=cfg.Tend_h;
        break;
    end

    % Exact integrated burdens over the holding interval.
    AUTC_NV_exact = AUTC_NV_exact + nV*dtEvent;
    AUTC_NT_exact = AUTC_NT_exact + (nV+nF)*dtEvent;

    % Records/snapshots strictly before event
    while nextRec<=nRec && tRec(nextRec)<tEvent-1e-12
        M = morphologyMetricsFromLists( ...
            state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);
        NV(nextRec)=M.NV; NF(nextRec)=M.NF; NT(nextRec)=M.NT;
        Rmax_sites(nextRec)=M.Rmax;
        Rboundary_mean_sites(nextRec)=M.RboundaryMean;
        roughness_sites(nextRec)=M.roughness;
        cavity_fraction(nextRec)=M.cavity;
        edgeOcc(nextRec)=M.edgeTouch;
        cumBirth(nextRec)=nBirthTotal;
        cumFerro(nextRec)=nFerroTotal;
        cumClear(nextRec)=nClearTotal;
        nextRec=nextRec+1;
    end

    while nextSnap<=numel(snapTimes) && snapTimes(nextSnap)<tEvent-1e-12
        snapshots(nextSnap).time_h = snapTimes(nextSnap);
        snapshots(nextSnap).state = reshape(state,[L L L]);
        nextSnap=nextSnap+1;
    end

    % Event class
    uClass = rand*A_total;
    chosenSite = uint32(0);
    chosenTarget = uint32(0);

    if uClass < A_birth
        % ================================================================
        % Birth
        % ================================================================
        edgeTicket = randi(totalDirectedVEmptyEdges);
        cs = cumsum(nEmptyV);
        parentPos = find(cs>=edgeTicket,1,'first');
        parent = vidx(parentPos);

        nbr = emptyNeighbourIndices(parent,state,L);
        assert(~isempty(nbr),'Selected birth parent has no empty neighbour.');
        target = nbr(randi(numel(nbr)));

        targetNbr = neighbourIndices(target,L);
        childEmptyCount = sum(state(targetNbr)==0);

        state(target)=uint8(1);

        % Every viable neighbour loses one empty neighbour.
        vNbr = targetNbr(state(targetNbr)==1);
        if ~isempty(vNbr)
            vals = double(nEmptyCount(vNbr))-1;
            assert(all(vals>=0),'Negative empty-neighbour count.');
            nEmptyCount(vNbr)=uint8(vals);
        end
        nEmptyCount(target)=uint8(childEmptyCount);

        nV=nV+1;
        Vlist(nV)=uint32(target);
        posV(target)=int32(nV);

        if strcmp(cfg.lambda_mode,'cell_mark')
            switch cfg.daughter_lambda_rule
                case 'inherit'
                    lambdaCell(target)=lambdaCell(parent);
                case 'sample_pool'
                    pool=cfg.lambdaF_pool_h(:);
                    lambdaCell(target)=pool(randi(numel(pool)));
                case 'site_map'
                    lambdaCell(target)=lambdaSite(target);
                otherwise
                    error('Unknown daughter_lambda_rule.');
            end
            assert(isfinite(lambdaCell(target)) && lambdaCell(target)>=0);
        end

        nBirthTotal=nBirthTotal+1;
        chosenSite=uint32(parent);
        chosenTarget=uint32(target);
        eType=uint8(1);

    elseif uClass < A_birth+A_ferro
        % ================================================================
        % Ferroptotic commitment
        % ================================================================
        if strcmp(cfg.lambda_mode,'scalar')
            p = randi(nV);
        else
            p = weightedChoice(ferroWeights);
        end

        idx = double(Vlist(p));

        % Remove from V set by swap-with-last.
        lastIdx = double(Vlist(nV));
        Vlist(p)=uint32(lastIdx);
        posV(lastIdx)=int32(p);
        Vlist(nV)=uint32(0);
        posV(idx)=int32(0);
        nV=nV-1;

        % Add to F set.
        nF=nF+1;
        Flist(nF)=uint32(idx);
        posF(idx)=int32(nF);
        state(idx)=uint8(2);

        if strcmp(cfg.lambda_mode,'cell_mark')
            lambdaCell(idx)=NaN;
        end

        nFerroTotal=nFerroTotal+1;
        chosenSite=uint32(idx);
        eType=uint8(2);

    else
        % ================================================================
        % Clearance
        % ================================================================
        p = randi(nF);
        idx = double(Flist(p));

        lastIdx = double(Flist(nF));
        Flist(p)=uint32(lastIdx);
        posF(lastIdx)=int32(p);
        Flist(nF)=uint32(0);
        posF(idx)=int32(0);
        nF=nF-1;

        state(idx)=uint8(0);

        nbr = neighbourIndices(idx,L);
        vNbr = nbr(state(nbr)==1);
        if ~isempty(vNbr)
            vals = double(nEmptyCount(vNbr))+1;
            assert(all(vals<=6),'Empty-neighbour count exceeded 6.');
            nEmptyCount(vNbr)=uint8(vals);
        end

        nEmptyCount(idx)=uint8(0);

        if strcmp(cfg.lambda_mode,'cell_mark')
            lambdaCell(idx)=NaN;
        end

        nClearTotal=nClearTotal+1;
        chosenSite=uint32(idx);
        eType=uint8(3);
    end

    t=tEvent;
    nEventTotal=nEventTotal+1;

    if nV==0 && ~isfinite(viableExtinctionTime_h)
        viableExtinctionTime_h=t;
    end

    if (nV+nF)==0 && ~isfinite(totalClearanceTime_h)
        totalClearanceTime_h=t;
    end

    if cfg.store_event_log
        eventTime(end+1,1)=t; %#ok<AGROW>
        eventType(end+1,1)=eType; %#ok<AGROW>
        eventSite(end+1,1)=chosenSite; %#ok<AGROW>
        eventTarget(end+1,1)=chosenTarget; %#ok<AGROW>
    end

    % Records/snapshots at or just after this event.
    while nextRec<=nRec && tRec(nextRec)<=t+1e-12
        M = morphologyMetricsFromLists( ...
            state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);
        NV(nextRec)=M.NV; NF(nextRec)=M.NF; NT(nextRec)=M.NT;
        Rmax_sites(nextRec)=M.Rmax;
        Rboundary_mean_sites(nextRec)=M.RboundaryMean;
        roughness_sites(nextRec)=M.roughness;
        cavity_fraction(nextRec)=M.cavity;
        edgeOcc(nextRec)=M.edgeTouch;
        cumBirth(nextRec)=nBirthTotal;
        cumFerro(nextRec)=nFerroTotal;
        cumClear(nextRec)=nClearTotal;
        nextRec=nextRec+1;
    end

    while nextSnap<=numel(snapTimes) && snapTimes(nextSnap)<=t+1e-12
        snapshots(nextSnap).time_h=snapTimes(nextSnap);
        snapshots(nextSnap).state=reshape(state,[L L L]);
        nextSnap=nextSnap+1;
    end

    if cfg.verbose && mod(nEventTotal,500)==0
        fprintf(['event=%7d | t=%9.4f h | NV=%6d NF=%6d NT=%6d | ', ...
            'B/F/C=%d/%d/%d\n'], ...
            nEventTotal,t,nV,nF,nV+nF, ...
            nBirthTotal,nFerroTotal,nClearTotal);
    end
end

%% Safety fill
while nextRec<=nRec
    M = morphologyMetricsFromLists( ...
        state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,cfg.compute_cavity);
    NV(nextRec)=M.NV; NF(nextRec)=M.NF; NT(nextRec)=M.NT;
    Rmax_sites(nextRec)=M.Rmax;
    Rboundary_mean_sites(nextRec)=M.RboundaryMean;
    roughness_sites(nextRec)=M.roughness;
    cavity_fraction(nextRec)=M.cavity;
    edgeOcc(nextRec)=M.edgeTouch;
    cumBirth(nextRec)=nBirthTotal;
    cumFerro(nextRec)=nFerroTotal;
    cumClear(nextRec)=nClearTotal;
    nextRec=nextRec+1;
end

while nextSnap<=numel(snapTimes)
    snapshots(nextSnap).time_h=snapTimes(nextSnap);
    snapshots(nextSnap).state=reshape(state,[L L L]);
    nextSnap=nextSnap+1;
end

%% Output
OUT = struct();
OUT.t_h=tRec;

OUT.NV=NV; OUT.NF=NF; OUT.NT=NT;

OUT.Rmax_sites=Rmax_sites;
OUT.Rboundary_mean_sites=Rboundary_mean_sites;
OUT.roughness_sites=roughness_sites;
OUT.cavity_fraction=cavity_fraction;
OUT.edgeOcc=edgeOcc;

OUT.cumBirth=cumBirth;
OUT.cumFerroptosis=cumFerro;
OUT.cumClearance=cumClear;

% Exact event-time integrals.
OUT.AUTC_NV=AUTC_NV_exact;
OUT.AUTC_NT=AUTC_NT_exact;

% Record-grid versions only as diagnostics.
OUT.AUTC_NV_recordGrid=trapz(tRec,NV);
OUT.AUTC_NT_recordGrid=trapz(tRec,NT);

OUT.viableExtinctionTime_h=viableExtinctionTime_h;
OUT.totalClearanceTime_h=totalClearanceTime_h;

OUT.snapshots=snapshots;

if cfg.return_final_state
    OUT.final_state=reshape(state,[L L L]);
else
    OUT.final_state=[];
end

if strcmp(cfg.lambda_mode,'cell_mark') && cfg.return_final_state
    OUT.final_lambda_cell_h=reshape(lambdaCell,[L L L]);
else
    OUT.final_lambda_cell_h=[];
end

if cfg.store_event_log
    OUT.event_log=table(eventTime,eventType,eventSite,eventTarget, ...
        'VariableNames',{'time_h','eventType','site','target'});
else
    OUT.event_log=table();
end

OUT.cfg=cfg;

OUT.algorithm=struct();
OUT.algorithm.name='exact event-driven direct-method SSA';
OUT.algorithm.time_step='none';
OUT.algorithm.birth_hazard='r*q_i';
OUT.algorithm.ferroptosis_hazard='lambda_F,i';
OUT.algorithm.clearance_hazard='mu';
OUT.algorithm.exact_for_specified_CTMC=true;
OUT.algorithm.nEvents=nEventTotal;

N0=NV(1); B=nBirthTotal; F=nFerroTotal; C=nClearTotal;

OUT.audit=struct();
OUT.audit.NV_identity=(nV==N0+B-F);
OUT.audit.NF_identity=(nF==F-C);
OUT.audit.NT_identity=((nV+nF)==N0+B-C);
OUT.audit.sum_identity=(NT(end)==NV(end)+NF(end));
OUT.audit.bookkeepingPASS = ...
    OUT.audit.NV_identity && OUT.audit.NF_identity && ...
    OUT.audit.NT_identity && OUT.audit.sum_identity;

end

%% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function cfg = validateCfg(cfg)

required = {'L','R0','Tend_h','record_dt_h','r_h','mu_h','seed','lambda_mode'};
for k=1:numel(required)
    assert(isfield(cfg,required{k}),'Missing required cfg field: %s',required{k});
end

assert(isscalar(cfg.L) && cfg.L>=7 && mod(cfg.L,2)==1 && cfg.L==round(cfg.L), ...
    'cfg.L must be an odd integer >=7.');
assert(isscalar(cfg.R0) && cfg.R0>0 && cfg.R0<(cfg.L-1)/2,'cfg.R0 is invalid.');
assert(isscalar(cfg.Tend_h) && isfinite(cfg.Tend_h) && cfg.Tend_h>0);
assert(isscalar(cfg.record_dt_h) && isfinite(cfg.record_dt_h) && cfg.record_dt_h>0);
assert(isscalar(cfg.r_h) && isfinite(cfg.r_h) && cfg.r_h>=0);
assert(isscalar(cfg.mu_h) && isfinite(cfg.mu_h) && cfg.mu_h>=0);

cfg.lambda_mode=lower(char(string(cfg.lambda_mode)));
allowed={'scalar','site_map','cell_mark'};
assert(any(strcmp(cfg.lambda_mode,allowed)), ...
    'cfg.lambda_mode must be scalar, site_map, or cell_mark.');

if ~isfield(cfg,'lambdaF_h') || isempty(cfg.lambdaF_h), cfg.lambdaF_h=NaN; end
if ~isfield(cfg,'lambdaF_site_h'), cfg.lambdaF_site_h=[]; end
if ~isfield(cfg,'lambdaF_initial_map_h'), cfg.lambdaF_initial_map_h=[]; end
if ~isfield(cfg,'lambdaF_pool_h'), cfg.lambdaF_pool_h=[]; end
if ~isfield(cfg,'daughter_lambda_rule') || isempty(cfg.daughter_lambda_rule)
    cfg.daughter_lambda_rule='inherit';
end
cfg.daughter_lambda_rule=lower(char(string(cfg.daughter_lambda_rule)));

if strcmp(cfg.lambda_mode,'scalar')
    assert(isscalar(cfg.lambdaF_h) && isfinite(cfg.lambdaF_h) && cfg.lambdaF_h>=0, ...
        'scalar mode requires cfg.lambdaF_h >= 0.');
end

if strcmp(cfg.lambda_mode,'site_map')
    assert(isequal(size(cfg.lambdaF_site_h),[cfg.L cfg.L cfg.L]));
    assert(all(isfinite(cfg.lambdaF_site_h(:))) && all(cfg.lambdaF_site_h(:)>=0));
end

if strcmp(cfg.lambda_mode,'cell_mark')
    rules={'inherit','sample_pool','site_map'};
    assert(any(strcmp(cfg.daughter_lambda_rule,rules)),'Invalid daughter_lambda_rule.');

    if isempty(cfg.lambdaF_initial_map_h)
        assert(isscalar(cfg.lambdaF_h) && isfinite(cfg.lambdaF_h) && cfg.lambdaF_h>=0);
    else
        assert(isequal(size(cfg.lambdaF_initial_map_h),[cfg.L cfg.L cfg.L]));
        assert(all(isfinite(cfg.lambdaF_initial_map_h(:))) && ...
               all(cfg.lambdaF_initial_map_h(:)>=0));
    end

    if strcmp(cfg.daughter_lambda_rule,'sample_pool')
        assert(~isempty(cfg.lambdaF_pool_h));
        assert(all(isfinite(cfg.lambdaF_pool_h(:))) && all(cfg.lambdaF_pool_h(:)>=0));
    end

    if strcmp(cfg.daughter_lambda_rule,'site_map')
        assert(isequal(size(cfg.lambdaF_site_h),[cfg.L cfg.L cfg.L]));
        assert(all(isfinite(cfg.lambdaF_site_h(:))) && all(cfg.lambdaF_site_h(:)>=0));
    end
end

if ~isfield(cfg,'initial_center') || isempty(cfg.initial_center)
    cfg.initial_center=[];
else
    assert(numel(cfg.initial_center)==3);
end

if ~isfield(cfg,'snapshot_times_h') || isempty(cfg.snapshot_times_h)
    cfg.snapshot_times_h=[];
end

if ~isfield(cfg,'return_final_state') || isempty(cfg.return_final_state)
    cfg.return_final_state=true;
end
cfg.return_final_state=logical(cfg.return_final_state);

if ~isfield(cfg,'compute_cavity') || isempty(cfg.compute_cavity)
    cfg.compute_cavity=true;
end
cfg.compute_cavity=logical(cfg.compute_cavity);

if ~isfield(cfg,'store_event_log') || isempty(cfg.store_event_log)
    cfg.store_event_log=false;
end
cfg.store_event_log=logical(cfg.store_event_log);

if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose=false;
end
cfg.verbose=logical(cfg.verbose);

end

function nbr = neighbourIndices(idx,L)

[i,j,k]=ind2sub([L L L],idx);
nbr=zeros(6,1);
n=0;

if i>1, n=n+1; nbr(n)=sub2ind([L L L],i-1,j,k); end
if i<L, n=n+1; nbr(n)=sub2ind([L L L],i+1,j,k); end
if j>1, n=n+1; nbr(n)=sub2ind([L L L],i,j-1,k); end
if j<L, n=n+1; nbr(n)=sub2ind([L L L],i,j+1,k); end
if k>1, n=n+1; nbr(n)=sub2ind([L L L],i,j,k-1); end
if k<L, n=n+1; nbr(n)=sub2ind([L L L],i,j,k+1); end

nbr=nbr(1:n);

end

function nbr = emptyNeighbourIndices(idx,state,L)
nbr=neighbourIndices(idx,L);
nbr=nbr(state(nbr)==0);
end

function p = weightedChoice(w)
w=double(w(:));
assert(~isempty(w) && all(isfinite(w)) && all(w>=0));
s=sum(w);
assert(s>0);
u=rand*s;
cs=cumsum(w);
p=find(cs>=u,1,'first');
end

function M = morphologyMetricsFromLists( ...
    state,Vlist,nV,Flist,nF,ctr,L,radiusSqAll,computeCavity)

M=struct();
M.NV=nV; M.NF=nF; M.NT=nV+nF;

if M.NT==0
    M.Rmax=0;
    M.RboundaryMean=0;
    M.roughness=0;
    M.cavity=0;
    M.edgeTouch=false;
    return;
end

occ=zeros(M.NT,1);
if nV>0, occ(1:nV)=double(Vlist(1:nV)); end
if nF>0, occ(nV+1:end)=double(Flist(1:nF)); end

[ii,jj,kk]=ind2sub([L L L],occ);

rr=sqrt((double(ii)-ctr(1)).^2 + ...
        (double(jj)-ctr(2)).^2 + ...
        (double(kk)-ctr(3)).^2);

M.Rmax=max(rr);

isBoundary=false(M.NT,1);
for q=1:M.NT
    nbr=neighbourIndices(occ(q),L);
    isBoundary(q)=any(state(nbr)==0);
end

rb=rr(isBoundary);
if isempty(rb)
    M.RboundaryMean=M.Rmax;
    M.roughness=0;
else
    M.RboundaryMean=mean(rb);
    M.roughness=std(rb,0);
end

if computeCavity
    inside = radiusSqAll(:)<=M.Rmax^2+1e-12;
    M.cavity=nnz(inside & state==0)/max(nnz(inside),1);
else
    M.cavity=NaN;
end

M.edgeTouch=any(ii==1 | ii==L | jj==1 | jj==L | kk==1 | kk==L);

end
