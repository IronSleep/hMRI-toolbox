function hmri_corrB1maps

    % Global hmri_def variable used across the whole toolbox
    global protocol_def

    % 1) IronSleep Protocol 7T previous
    protocol_def.b1acq_set.tags{1}  = 'KP_AFI_1f_7T';
    protocol_def.b1acq_set.tr{1} = [0.025, 0.125]; % TRs in seconds
    protocol_def.b1acq_set.fa{1} = 55; % flip angle in degrees
    protocol_def.b1acq_set.rfsp_angle{1} = 36; % RF spoiling angle in degrees
    protocol_def.b1acq_set.magnetic_field{1} = 7; % magnetic field strength in Tesla
    protocol_def.b1acq_set.seqname{1} = 'kpafi1f';
    protocol_def.b1acq_set.p{1} = [0.000010357851256, -0.001837088945544, 1.037456588545633, 1.159931060627856]; % correction factors for 7T pTx and sTx 
    
    % 2) IronSleep Protocol 7T 2025
    protocol_def.b1acq_set.tags{2}  = 'KP_AFI_1g_7T';
    protocol_def.b1acq_set.tr{2} = [0.025, 0.125]; % TRs in seconds
    protocol_def.b1acq_set.fa{2} = 55; % flip angle in degrees
    protocol_def.b1acq_set.rfsp_angle{2} = 36; % RF spoiling angle in degrees
    protocol_def.b1acq_set.magnetic_field{2} = 7; % magnetic field strength in Tesla
    protocol_def.b1acq_set.seqname{2} = 'kpafi1g';
    protocol_def.b1acq_set.p{2} = [0.000010195724382, -0.001798060630655, 1.032633895145047, 1.164429706935533]; % correction factors for 7T pTx and sTx 
    

    % 3) IronSleep Protocol 3T 2025
    protocol_def.b1acq_set.tags{3}  = 'KP_AFI_3T';
    protocol_def.b1acq_set.tr{3} = [0.05, 0.15]; % TRs in seconds
    protocol_def.b1acq_set.fa{3} = 60; % flip angle in degrees
    protocol_def.b1acq_set.rfsp_angle{3} = 129.3; % RF spoiling angle in degrees
    protocol_def.b1acq_set.magnetic_field{3} = 3; % magnetic field strength in Tesla
    protocol_def.b1acq_set.seqname{3} = 'kpafi1g';
    protocol_def.b1acq_set.p{3} = [-0.000003910876391, -0.000117577121954, 1.122587693334678, -0.168746536917551]; % correction factors for 3T sTx

end