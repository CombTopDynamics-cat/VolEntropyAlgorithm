# VolEntropyAlgorithm
This repository contains the MAPLE procedure to compute the Volume Entropy
of geometric presentations of surface groups trough a dynamical approach
developed in [Alsedà et al. 2025].

The dynamical aproach consists in computing the topological entropy of the
Bowen-Series-like maps (a family of discontinuous maps on the circle)  
using the Milnor-Thurston theory for one dimensional maps.

The topological entropy of the Bowen-Series-like maps coincides with
the volume entropy of the group as proved in [Alsedà et al. 2025]

[Alsedà et al. 2025]
ENTROPY STABILITY AND MILNOR-THURSTON INVARIANTS FOR BOWEN-SERIES-LIKE MAPS
Lluís Alsedà, David Juher, Jérôme Los and Francesc Mañosas
Ergodic Theory and Dynamical Systems, 2025.

Next we show the code of the main program, as a summary and description of the algorithm

VolumeEntropy := proc(R)

    # Computes the volume entropy of the presentation R if R is geometric
    
    # or reports that R is not geometric otherwise
    
    # USES: LinearAlgebra library
    
    # USES: CheckRelations,CyclicOrdering,MinimalBigons,KneadingMatrix procedures
    
    local co,bigons,i,A:
    
    if CheckRelations(R)=false then
    
       print(‘The relations do not satisfy the syntax conventions‘): return:
       
    end if:
    
    co:=CyclicOrdering(R):
    if co=false then print(‘The presentation is not geometric‘): return: end if:
    bigons:=MinimalBigons(R,co):
    A:=KneadingMatrix(co,bigons,R):
    i:=Determinant(DeleteColumn(A,1)):
    return 1/min(fsolve(i,t=0..1));
end proc:
