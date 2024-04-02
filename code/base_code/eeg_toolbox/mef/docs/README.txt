
Help on MEX functions:
	MEX functions cannot have help like normal MATLAB functions. However, by providing a MATLAB function header with the same name as the MEX function and the desired help information, help can be provided. This is done for the MEF-related MEX functions. Just type `help <function_name>` in MATLAB to access the help.


Changes to mexopts.sh necessary for compilation:
	- remove -ansi option from CXX flags for the machine type in question (in the case of this server, glnxa64). This means that the C compiler will allow the use of C++-style comments and the inline keyword, both of which are used in the code. See the file "mexopts.sh.example" in this folder for a working (as of 2011.07.05) example.
