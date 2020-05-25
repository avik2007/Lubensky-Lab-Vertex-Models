# Lubensky-Lab-Vertex-Models

2D Vertex Model code written by members of David Lubensky's group at the University of Michigan.

This code is written in MATLAB and C++. MATLAB is used to create input files, which 
describe the cell packing and the physical parameters describing the model. The 
C++ code is used to simulate the dynamics of the packings. The CPP Code Base contains
the base classes for the C++ code. For any of the multiple variants developed by
the lab, you can also pull the other folders which contain subclass files.
Descriptions of these subclasses are attached are included as a textfile. The subclass 
files do not contain the parent classes, so you need to pull CPP Code Base AND whichever 
subclass you are interested in implementing. 

Base Classes:

Vertices: contains methods that create vertices and calculate important vertex characteristics

Edges: contains methods that create edges and calculate important edge characteristics

Cells: contains methods that create cells and calculate their characteristic quantities 

Tissue: contains methods that manipulate the cell, edge, and vertex objects, 
	i.e. contains numerical integrator for simulation, contains mitosis and apoptosis 
	
Vect: contains methods that allow for use of vectors and vector algebra

Works that Use this Code Base:
1. Spencer, Meryl A., Jabeen, Zahera, and Lubensky, David K.; Vertex stability and topological
transitions in vertex models of foams and epithelia," Eur. Phys. J. E, vol. 40, no. 1, p. 2, 2017.

2. Nunley, H., Nagashima, M., Martin, K., Gonzalez, A. L., Suzuki, S. C., Norton, D., Wong,
R. O. L., Raymond, P. A., Lubensky, D. K.; Defect patterns on the curved surface of fish retinae suggest mechanism of cone mosaic formation. bioRxiv, 2019.

Works Cited:
1. A. G. Fletcher, J. M. Osborne, P. K. Maini, and D. J. Gavaghan; Implementing vertex dynamics
models of cell populations in biology within a consistent computational framework," Progress in
Biophysics and Molecular Biology, vol. 113, no. 2, pp. 299 { 326, 2013.

2. S. Alt, P. Ganguly, and G. Salbreux; "Vertex models: from cell mechanics to tissue morphogenesis,"
Philosophical Transactions Of The Royal Society B-Biological Sciences, vol. 372, May 2017.

3. D. L. Barton, S. Henkes, C. J. Weijer, and R. Sknepnek; "Active Vertex Model for cell-resolution
description of epithelial tissue mechanics," PLOS Computational Biology, vol. 13, Jun 2017.

4. S. Curran, C. Strandkvist, J. Bathmann, M. de Gennes, A. Kabla, G. Salbreux, and B. Baum;
"Myosin II Controls Junction Fluctuations to Guide Epithelial Tissue Ordering," Developmental
Cell, vol. 43, pp. 480+, Nov 2017.

5. D. Bi, X. Yang, M. C. Marchetti, and M. L. Manning; Motility-Driven Glass and Jamming
Transitions in Biological Tissues," Phys. Rev. X, vol. 6, Apr 2016.

6. D. M. Sussman; cellGPU: Massively parallel simulations of dynamic vertex models," arXiv, 2017.

7. Noll, N., Mani, M., Hemmskerk, I., Streichan, S. J., & Shraiman, B. I.;"Active Tension Network model
suggests an exotic mechanical state realized in epithelial tissues." Nature physics, 13(12, 1221-1226. https://doi.org/10.1038/nphys4219. (2017)

8. Honda, H., Tanemaru, M., Nagai, T.; "A three-dimensional vertex dynamics cell model of space-filling polyhedra simulating cell behavior in a cell aggregate." Journal of Theoretical Biology,
vol 226, issue 4. February 21, 2004.

9. Honda, H.; Geometrical models for cells in tissues." Int. Rev. Cytol. 81, 191 (1983).

10. Bessel function code is adapted from https://www.atnf.csiro.au/computing/software/gipsy/sub/bessel.c. Author: M. G. R. Vogelaar at the Kapteyn Institute Groningen. 

Subclass Notes:

AGVertex - "Anisotropic Growth Vertex." AGVertex contains subclasses that allow the vertex modelto handle
on-average anisotropic growth. These subclasses are in the process of introducing a new NPT Browian dynamics model that evolves the periodic box dimensions
more correctly than the integrator in the base classes. 

MegaCellSimulation - This code is for simulation an annular packing of cells. There is a large cell in the interior (with different mechanical properties). There is a large cell on the outside of the annulus. This large cell on the outside is weird because it is extends from the edge of the annulus to each of the periodic boundaries (i.e., +x, -x, +y, -y).

PCPSimulation - This code is for simulating the dynamics of PCP proteins that live on the edges of cells. Solving the dynamics of PCP requires an implicit-explicit solver (because the equations are "stiff").

PCPMegaCellSimulation - This code is for simulating annular cell packing with PCP proteins that live on the edges of cells. (This effectively is a combination of MegaCellSimulation and PCPSimulation)

Fiber - This code allows fibers to nucleate at TCJs and move across the cell until the fiber breaks. In these simulations, the cells themselves are stationary. 
