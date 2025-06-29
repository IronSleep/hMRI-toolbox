function hmri_corrB1maps

    % Global hmri_def variable used across the whole toolbox
    global protocol_def

    % 1) IronSleep Protocol 7T 2025
    protocol_def.b1acq_set.tags{1}  = 'KP_AFI_7T';
    protocol_def.b1acq_set.tr{1} = [0.025, 0.125]; % TRs in seconds
    protocol_def.b1acq_set.fa{1} = 55; % flip angle in degrees
    protocol_def.b1acq_set.rfsp_angle{1} = 36; % RF spoiling angle in degrees
    protocol_def.b1acq_set.seqname{1} = 'kpafi1f';
    protocol_def.b1acq_set.p{1} = [0.000010195724382, -0.001798060630655, 1.032633895145047, 1.164429706935533]; % correction factors for 7T pTx and sTx 
    
    % 2) IronSleep Protocol 3T 2025
    protocol_def.b1acq_set.tags{2}  = 'KP_AFI_3T';
    protocol_def.b1acq_set.tr{2} = [0.05, 0.15]; % TRs in seconds
    protocol_def.b1acq_set.fa{2} = 60; % flip angle in degrees
    protocol_def.b1acq_set.rfsp_angle{2} = 129.3; % RF spoiling angle in degrees
    protocol_def.b1acq_set.seqname{2} = 'kpafi1g';
    protocol_def.b1acq_set.p{2} = [-0.000003910876391, -0.000117577121954, 1.122587693334678, -0.168746536917551]; % correction factors for 3T sTx

end