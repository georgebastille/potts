# q-state Potts Model Simulation

A Java™ implementation of the q-state Potts Model on a 2D square lattice is used, which confirmed analytically known results for q = 3 and q = 5. The Metropolis single spin flip and Wolff cluster flip algorithms are utilized to find phase changes, investigate Finite Size Scaling effects and calculate the Heat Capacity Critical Exponent.

For details see the [report](potts_report.doc)

To run the simulation:
```java
cd code/Graphical
javac LatticeSimulator.java
java LatticeSimulator
```
