function [ measured ] = WFObs_s_loadmeasurements( strucObs, time )
% WFOBS_S_LOADMEASUREMENTS  Loads measurement data for estimation
%
%   SUMMARY
%    This script loads SOWFA data files and formats them to the correct
%    format for usage with WFObs. Inputs and outputs are:
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - strucObs: this struct contains all the observer settings and
%       (temporary) files used for updates, such as covariance matrices,
%       ensemble/sigma point sets, measurement noise, etc.
%
%     - time: current time instant in s (belonging to datafile number).
%
%     - measured: measurement data struct with the following entries
%        *measured.power  output power of turbines, undisturbed
%        *measured.sol    velocities in format of 'x' in 'Ax=b', disturbed
%        *measured.solq   velocities in format of 'x' in 'Ax=b', undisturbed
%        *measured.u      longitudinal velocities, disturbed
%        *measured.uq     longitudinal velocities, undisturbed
%        *measured.v      lateral velocities, disturbed
%        *measured.vq     lateral velocities, undisturbed
%

% Import variables
sourcepath   = strucObs.measurementsPath;
datanroffset = strucObs.measurementsOffset;

% Load the file and nullify NaN entries
measured                        = load([sourcepath '/' num2str(datanroffset+time) '.mat']);
measured.uq(isnan(measured.uq)) = 0; % nullify any out-of-range measurements
measured.vq(isnan(measured.vq)) = 0; % nullify any out-of-range measurements

% Add noise to measurements that are to be fed into the observer
measured.u = measured.uq + strucObs.noise_obs*randn(size(measured.uq));
measured.v = measured.vq + strucObs.noise_obs*randn(size(measured.vq));

% Sort them into a vector of correct size (for Ax=b)
measured.solq = [vec(measured.uq(3:end-1,2:end-1)'); vec(measured.vq(2:end-1,3:end-1)')];
measured.sol  = [vec(measured.u(3:end-1,2:end-1)') ; vec(measured.v(2:end-1,3:end-1)') ];
end