This project investigates the decay of ortho-positronium—a bound state of an electron and a positron—into three photons. The focus is on analyzing the kinematics of this three-photon decay using theoretical predictions from special relativity and comparing them with experimental data.
Experiments were conducted using NaI scintillators and Na-22 radioactive sources to detect the emitted photons. Various decay configurations were explored to determine their consistency with relativistic predictions.
As a final outcome, a Dalitz plot was produced. This is a 2D representation of the phase space, showing the probability distribution of different decay configurations based on two of the three angles between the emitted photons.

## Experimental Set Up

To detect the 3 photons resulting from the decay, 3 NaI scintillators were placed in the same plane around 2 Na-22 sources, from which ortho-positronium is created. 
The 3 detectors are connected to a Faster acquisition board, which allows recording the data from the 3 detectors when they are in coincidence—meaning all three detect a signal within a 400 ns time window. 
The angles between each detector vary depending on the studied configurations, and data is collected for at least 15 hours per configuration.
Lead shielding was added around each detector to reduce the detection of photons resulting from Compton scattering occurring inside another detector.

## Configurations of the decay from a probabilist view : the Dalitz Plot

The O-Ps decays in 3 photons, which trajectories are all in a single plane because of momentum conservation. The 3-body decay can have multiple possible configurations, and the goal was to determin which ones were possible.

The data measured are energies and time, in bunches of 3 because they are coincidence measurements with 3 detectors. It is easy to find the corresponding angle configuration of the decay from these energies (energy_to_angle.py) in order to
plot a Dalitz plot of the decay, which is a probabilistic vizualization of the decay configurations. The theoretical Dalitz plot is built with a high number of possible configurations of the decay (computed based on energy and momentum conservation)
in Dadlitz_plot_angles.py.






