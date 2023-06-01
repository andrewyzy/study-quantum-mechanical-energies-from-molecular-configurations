# Study: quantum mechanical energies from molecular configurations
Just some simple simulations about quantum mechanical energies from molecular configurations. 

The code is divided into two parts: a Molecule class and a System class.

The Molecule class represents a single molecule. When we create a new molecule, we specify the maximum energy level it can have (max_energy_level), and the position of the molecule (position).

The Molecule class has two methods:

absorb_light(): When a molecule absorbs light, it gains energy. The amount of energy gained is passed as an argument to this function. If the added energy would take the molecule beyond its maximum energy level, it sets the molecule's energy level to the maximum and returns the excess energy.
transfer_energy(): This method is used to add energy to a molecule. This can happen when another molecule has excess energy to get rid of. If the added energy would take the molecule beyond its maximum energy level, it prints a message saying that it can't accept the extra energy.

The System class represents a system of multiple molecules. It contains a list of Molecule objects (molecules).
The System class also has two methods:

add_molecule(): This function is used to add a new molecule to the system.
pulse_light(): This function simulates a pulse of light hitting the entire system. It goes through each molecule in the system, makes it absorb the light, and if there's excess energy, transfers it to the next molecule.
When you create a System and add Molecules to it, you can simulate how the molecules absorb and transfer energy in response to pulses of light. If a molecule receives too much energy from a light pulse, it transfers the excess energy to the next molecule in the system. 
