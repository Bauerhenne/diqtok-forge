# Compiler settings
CXX = g++
CXXFLAGS = -O2
MPICXX = mpicxx
MPICXXFLAGS = -D MPIusage -O2


# Rules
default : integrate.o qtokentools.o analyse_pf2mopt \
	pb_pf_0_1opt_3DIT_123ML pf_2mBayesian pf_3mBayesian \
	pf_2mfix pf_3mfix pf_2mopt


integrate.o : integrate.cpp
	$(CXX) $(CXXFLAGS) -c integrate.cpp

qtokentools.o : qtokentools.cpp
	$(CXX) $(CXXFLAGS) -c qtokentools.cpp

analyse_pf2mopt : analyse_pf2mopt.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o analyse_pf2mopt integrate.o qtokentools.o analyse_pf2mopt.cpp

pb_pf_0_1opt_3DIT_123ML : pb_pf_0_1opt_3DIT_123ML.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o pb_pf_0_1opt_3DIT_123ML integrate.o qtokentools.o pb_pf_0_1opt_3DIT_123ML.cpp

pf_2mBayesian : pf_2mBayesian.cpp integrate.o qtokentools.o
	 $(CXX) $(CXXFLAGS) -o pf_2mBayesian integrate.o qtokentools.o pf_2mBayesian.cpp

pf_3mBayesian : pf_3mBayesian.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o pf_3mBayesian integrate.o qtokentools.o pf_3mBayesian.cpp

pf_2mfix : pf_2mfix.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o pf_2mfix integrate.o qtokentools.o pf_2mfix.cpp

pf_3mfix : pf_3mfix.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o pf_3mfix integrate.o qtokentools.o pf_3mfix.cpp

pf_2mopt: pf_2mopt.cpp integrate.o qtokentools.o
	$(CXX) $(CXXFLAGS) -o pf_2mopt integrate.o qtokentools.o pf_2mopt.cpp


mpi : integrate.o qtokentools.o analyse_pf2mopt \
	pb_pf_0_1opt_3DIT_123ML_MPI pf_2mBayesian_MPI pf_3mBayesian_MPI \
	pf_2mfix_MPI pf_3mfix_MPI pf_2mopt_MPI


pb_pf_0_1opt_3DIT_123ML_MPI : pb_pf_0_1opt_3DIT_123ML.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pb_pf_0_1opt_3DIT_123ML integrate.o qtokentools.o pb_pf_0_1opt_3DIT_123ML.cpp

pf_2mBayesian_MPI : pf_2mBayesian.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pf_2mBayesian integrate.o qtokentools.o pf_2mBayesian.cpp

pf_3mBayesian_MPI : pf_3mBayesian.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pf_3mBayesian integrate.o qtokentools.o pf_3mBayesian.cpp

pf_2mfix_MPI : pf_2mfix.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pf_2mfix integrate.o qtokentools.o pf_2mfix.cpp

pf_3mfix_MPI : pf_3mfix.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pf_3mfix integrate.o qtokentools.o pf_3mfix.cpp

pf_2mopt_MPI: pf_2mopt.cpp integrate.o qtokentools.o
	$(MPICXX) $(MPICXXFLAGS) -o pf_2mopt integrate.o qtokentools.o pf_2mopt.cpp


clean :
	rm *.o analyse_pf2mopt pb_pf_0_1opt_3DIT_123ML pf_2mBayesian pf_3mBayesian pf_2mfix pf_3mfix pf_2mopt
