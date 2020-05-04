---
title: 'LB Colloids: : Software package for simulating 2-dimensional 
microscale colloid transport with Lattice Boltzmann'  
tags:  
   - Python
   - FORTRAN
   - Lattice Boltzmann
   - Colloid Transport
   - Fluid Dynamics    

authors:
   - name: Joshua D. Larsen  
   orcid: 0000-0002-1218-800X  
   affiliation: 1  
   
affiliations:
   - name: The University of Arizona  
   index: 1  

date: 28 April 2020  
bibliography: paper.bib
---

# Summary

Knowledge of colloid mobility in the soil environment is important for 
understanding nutrient cycling, the transport of some contaminants, and for 
developing environmental remediation systems such as geologic filters. The 
interaction forces between colloids and soil materials are central to colloid 
transport and immobilization. These forces act at the microscale (nanometers 
to microns) and include: fluid drag (friction), Brownian motion, gravity and 
buoyancy, and fluid chemical forces (including DLVO and van der Waals 
mechanisms). Most vadose zone studies, however, consider colloids at the 
continuum scale in terms of solute transport mechanisms using parametrized 
forms of the advection-dispersion equation and absorption isotherms. 

A small number of pore scale models have also been developed to track colloid
transport in porous media (Redman and others, 2004; Gao and others, 2010; Qui and
others, 2011). These models either use Lagrangian mechanics (which are computationally
inefficient for large numbers of colloids) or exist as novel approaches to modeling micro
scale colloid-surface interactions. Furthermore, none of approaches have been released as
open source tools. The limitations of these systems leave the interdisciplinary researcher
without a practical option to gain additional insight into the mechanisms driving the
physiochemical dynamics of colloid transport within their system.

LB-Colloids (Lattice Boltzmann Colloids) is an open-source, 2-dimensional modeling
system, designed to simulate colloid transport under steady state fluid flow
condions. Because colloids are very small particles, generally 1 micron to 1 
nanometer in diameter, it is difficult and labor intensive to observe colloid
tranport processes and colloid immobilization in porous media (Torkzaban and others, 2008).
Pore-scale modeling of the fundamental forces acting upon a colloids not only 
assists in the interpretation of experimental studies regarding colloid mobility,
but it can also provide insight into the physiochemical processes controlling
transport at the macroscale (Larsen, 2018). Gao and others (2010) presented a
colloid modeling approach which was updated by Qiu and others (2011) that 
included physiochemical transport mechanisms. LB-Colloids is novel because unlike
previous approaches chemical parameters such as acid base exchange parameters 
and van der Waals surface tension measurements do not need to be parameterized.
Instead a simplified formulation is used that avoids uncommonly collected field 
and laboratory measurements.

The overall process for performing a colloid particle tracking model with 
LB-Colloids consists of three steps. A pore structure at a resolution of 1µm 
or less must be obtained by imaging techniques or by generating synthetic 
media, such as with the Penetrable-Sphere (`PSphere`) module. A saturated flow model 
is then performed with the included LB code or is imported from a suitable 
CFD simulation. Finally, the transport of individual colloids is simulated in 
the pore structure by calculating the force balance and tracking each colloid 
over time. Since some chemical and physical forces act over nanometer scale 
distances, grid refinement techniques are included and can be used to more 
accurately represent colloid transport near fluid-solid interfaces. 

### Penetrable-Sphere
Penetrable-Sphere (`PSphere`) is an algorithm that was developed to generate
synthetic porous media to perform computational fluid dynamic simulations in situations
where segmented CT imagery of the porous media is unavailable. The artificially
generated porous media resembles a natural granular medium while having a
porosity and surface area that can be specified a-priori. `Psphere` creates
synthetic porous media by generating spheres and placing them within a user 
defined image domain using gaussian fields (Larsen, 2018). The artificial porous media is 
automatically checked to assure that the porosity and surface area are within
user defined tolerances and that the porous media percolates. The media is 
automatically regenerated and only returned to the user when the artificial 
porous media meets these criteria. 

### Lattice Boltzmann 2D
Lattice Boltzmann CFD methods are able to represent complex geometries
accurately through the discretization process and the application of simple bounce back
rules. By applying either a body force or pressure boundary conditions (Zou and He, 1997) 
fluid flow can be simulated with the system. Simulations have been 
successfully used to represent fluid flow in saturated systems 
(Blunt and others, 2013), unsaturated systems (Porter and others,
2009), heat transport (He and others, 1998), and macropore fluid flow (Sukop and others,
2013). The `LB2DModel` solves the 2-dimensional, nine fluid node (D2Q9) lattice 
Boltzmann equation. Fluid particle streaming and collision is achieved by 
using the Bhantager, Gross, and Krook (BGK) (1954) solution to the lattice 
Boltzmann equation (Qian and others, 1992).

### Colloids 

Colloid transport through the soil environment in an import contributor
to soil development processes through the aggregation of clays (Bronick and Lal 2005),
contaminant transport (Saiers and Hornberger, 1996, Jaisi and others, 2008), filtration
and transport of bio-colloids (Harter and Wagner, 2000, Redman and others, 2004,
Foppen and others, 2005), and soil nutrient dynamics (Bradford and others 2008).

Most vadose zone studies, however, consider colloids at the continuum
scale in terms of solute transport mechanisms using parametrized forms of the 
advection-dispersion equation and absorption isotherms. A comprehensive, 
well-documented and publicly available framework for simulating colloids at the
microscale is still lacking.
 
 
   
