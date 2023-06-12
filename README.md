This repo contains two simulation files: 
  - reciprocatingPistonCompressor_t.m
  - galinCompressor_t.m

The goal of this work is to compare the theoretical performance of a reciprocating piston compression to that of the Galin engine compression stroke. We are doing this 
comparison to figure out how the Galin architecture stacks-up relative to the reciprocating piston based compression. 

Compression parameters:
- polytropic index: 1.3
- suction pressure: 100 kPa
- suction temperature: 295 K
- discharge pressure: 689.6 kPa (100 psi)
- intake volume: 900 mL

Common material properties: 
- steel density: 7800 kg/m^2

Reciprocating piston geometry:
- bore: 94 mm
- stroke: 101 mm
- rod: 202 mm
- piston area: 6893 mm^2
- piston mass: 9.7 kg

Galin compressor geometry:
- vane radius: 125 mm
- shaft radius: 42 mm
- vane area: 6893 mm^2 
- vane mass: 9.7 kg

Running the simulations with the above parameters give us the following outputs:

**Common output parameters**
* compression ratio []: 4.4
* work done on the gases during compression [J]: 168.4
* pressure difference, start to end [kPa]: 100.000000 -> 689.5 [14.5 -> 100.0 (psi)]
* volume difference, start to end [mL]: 900.0 -> 203.8
* temperature difference, start to end [K]: 295.0 -> 460.6

Reciprocating piston compression stroke performance:
* stroke time [ms]: 56.4
* torque requested from electrical machine as 0.8 of max [Nm]: 69.4
* based on the torque load on EM machine * 0.8 we calculate the power of the EM [kW]: 3.9
* CFM: 16.9
* CFM/HP: 4.4

Galin compression stroke performance:
* stroke time [ms]: 71.3
* torque requested from electrical machines [Nm]: 138.6
* electrical machine power req. (each) [kW]: 3.1
* CFM: 53.5
* CFM/HP: 8.8

Notes: 
* in the case of the reciprocating piston simulation, we calculate the force required to compress the gases as: work required to compress gases / stroke,
* we perform a sanity check to make sure that in both simulations the total work done by the piston and vanes on the gases is equal to the required work to compress the gases,
 * we provide a plot versus time under pics/torque_vs_time.png. This plot shows that the requested torque (while higher in magnitude) from the electrical machines is constant in the gase of Galin compressor, whereas in the reciprocating piston case it varies wildly as the lever arm changes.

Conclusions:
* the simulated CFM/HP for the reciprocating piston engine is close to the often quoted figure of 3.5 (at 100psi) in industry, so our simulation isn't too unrealistic,
* the compression stroke is 25% quicker in the case of the reciprocating piston geometry - the advantages of that are closer to adiabatic conditions with less heat lost externally,
* however, the CFM is three times greater in the case of the Galin engine which is a huge improvement,
* and CFM/HP is double even when taking into consideration that in the Galin architecture we require two electrical machines.
