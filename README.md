# Continuability

We provide here additional material related to the paper: Do moldable applications perform better on failure-prone HPC platforms? by Valentin Le Fèvre, George Bosilca, Aurélien Bouteiller, Thomas Hérault, Atsushi Hori, Yves Robert, Jack Dongarra. This is joint-work by INRIA (FR), ICL/UTK (TN,USA) and RIKEN (JP).

The source code is composed of one C++ file where several functions compute the expected Time/Work/Yield for the 3 types of applications: Rigid, Moldable, GridShaped. An application using no spares can be simulated by setting the number of allowed failures to 0 in one the 3 model. The functions used for collecting the data for the different figures in the paper are also present.

### Compilation

g++ equations.cpp -o [out]

### Running

./[out] [[-nmwc value ]] -M [model] -[PCWYLST]

The different available parameters are:
  -n val : set the number of processors to val. By default: 22,500 (150x150).
  -c val : set the checkpointing time to val (in seconds). By default: 2 min.
  -w val : set the wait time to val (in seconds). By default: 12 hours.
  -m val : set the MTBF to val (in seconds). By default: 20 years.
  -M model : set the model to the chosen one (rigid|moldable|grid). By default: Rigid.
  -P : evaluate the models while varying the number of processors.
  -C : evaluate the models while varying the checkpoint cost.
  -W : evaluate the models while varying the wait time.
  -Y : evaluate the models while varying the MTBF.
  -L : plot a table of values related to the yield as a function of the number of allowed faults (for a given model).
  -S : make all parameters vary to describe a spectrum of results.
  -T : used for testing purposes.
  
  
