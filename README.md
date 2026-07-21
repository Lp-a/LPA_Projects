Code samples from my PhD research in quantum sensing with diamond NV centers (Backlund Lab, UIUC), shared here as representative examples of instrumentation control, data analysis, and simulation work.

Repository structure
Example_Instrumentation_Control/   Hardware control software
Example_data_analysis/              Experimental data processing and fitting
Example_simulations/                Physics simulations used to guide experiment design

Example_Instrumentation_Control
Cobolt06MLD_Lakshmy.py — Python driver for serial communication and control of a Cobolt 06-MLD laser, including port/serial-number discovery, connection handling, and device command interfacing. Used to integrate the laser into our optical setup for automated experiment control.

Example_data_analysis
MATLAB scripts for processing raw experimental data from our NV-center measurements.

DEER_data_analysis.m — Loads signal/reference data across repeated iterations, computes contrast, and fits a decaying-sinusoid model to extract relaxation and oscillation parameters.
DEER_T1_2026_April8.m — T1 relaxometry analysis pipeline: aggregates signal/reference data across iterations and time points to extract T1 decay behavior.
DEER_Rabi_Analysis.m — Processes DEER Rabi oscillation data, computing contrast from signal/reference traces as a function of pulse sweep time.
RTIL_bulk_EPR_experimental_data_fit.m — Loads and normalizes bulk EPR experimental spectra and fits them against a simulated spin Hamiltonian model (using EasySpin) to extract spin system parameters (g-tensor, hyperfine coupling).

Example_simulations
MATLAB scripts simulating expected experimental signals to guide measurement design and interpretation.

Bulk_EPR_simulation_at_two_diff_Tc_with_easyspin.m — Simulates bulk EPR spectra for a nitroxide radical spin system at different tumbling correlation times using EasySpin, based on literature spin-Hamiltonian parameters.
T1_cross_relaxation_withmagfield_simulation.m — Models NV T1 relaxation via cross-relaxation with a nearby electron spin bath as a function of magnetic field and spin density.
calculate_xyySushkovCorr1pulse_GaussianfromBlochEqns_w13Cbg.m — Simulates single-pulse DEER Rabi signal from the Bloch equations for a driven spin system with T1/T2 relaxation.
simulate_Mamin_2012_fig3.m — Reproduces the surface-spin-correlation analysis from Mamin et al. (2012), modeling NV relaxometry sensitivity to surface spin baths for different T1/field conditions.
