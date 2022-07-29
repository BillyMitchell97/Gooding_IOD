# Gooding_IOD
Original Gooding angles-only initial orbit determination

# Fortran Files:
- extra/Gooding_master.f = The first version of the code that has many print statements for debugging purposes
- extra/Gooding_legacy.f =  A legacy version of the code that maintains some Fortran77 standards that are now obsolete
- Gooding_runtime.f = The most up-to-date version of the code that reproduces the results of Gooding's original 1993 paper
- in.data = The input file for the fortran code

# Notes on code:
- Originally written to input center of force -> observing site vector and the observed direction affine unit vector.
- Also requires time between observations as inputs
- Outputs the range estimates for observation 1 and 2
- inputs and outputs are unit-agnostic


TO DO:
[ ] Create working python version of gooding_runtime
