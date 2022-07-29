# Gooding_IOD
Original Gooding angles-only initial orbit determination [1-4]

# Fortran Files:
- extra/Gooding_master.f = The first version of the code that has many print statements for debugging purposes
- extra/Gooding_legacy.f =  A legacy version of the code that maintains some Fortran77 standards that are now obsolete
- Gooding_runtime.f = The most up-to-date version of the code that reproduces the results of Gooding's original 1993 paper [3]

# Other Files:
- in.data = The input file for the fortran code
- gooding_results.xlsx = Excel spreadsheet showing the comparison of Gooding's original 1993 paper [3] and my current code results

# Notes on code:
- Originally written to input center of force -> observing site vector and the observed direction affine unit vector.
- Also requires time between observations as inputs
- Outputs the range estimates for observation 1 and 2
- inputs and outputs are unit-agnostic


# TO DO:
[ ] Create working python version of gooding_runtime

references:<br>
[1] R. H. Gooding, “On Universal Elements, and Conversion Procedures to and from Position and Velocity”, _Celestial Mechanics_ and Dynamical Astronomy **44** (1988)<br>
[2] R. H. Gooding, “A Procedure for the Solution of Lambert’s Orbital Boundary-Value Problem”, _Celestial Mechanics and Dynamical Astronomy_ **48** (1990)<br>
[3] R. H. Gooding, “A New Procedure for Orbit Determination Based on Three Lines of Sight (Angles Only)”, _Defense Research Agency_, Technical Report **93004** (1993)<br>
[4] R. H. Gooding, “A New Procedure for the Solution of the Classical Problem of Minimal Orbit Determination from Three Lines of Sight”, _Celestial Mechanics and Dynamical Astronomy_ **66** (1997)<br>
