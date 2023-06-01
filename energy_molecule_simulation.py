import math
import random

class Molecule:
    def __init__(self, max_energy_level, position, orientation):
        self.max_energy_level = max_energy_level
        self.current_energy_level = 0
        self.position = position  # position as a tuple (x, y, z)
        self.orientation = orientation  # orientation as a unit vector (dx, dy, dz)

    def absorb_light(self, energy):
        potential_energy_level = self.current_energy_level + energy
        if potential_energy_level <= self.max_energy_level:
            self.current_energy_level = potential_energy_level
            print(f"Molecule at position {self.position} has moved to energy level {self.current_energy_level}")
            return 0  # no excess energy
        else:
            excess_energy = potential_energy_level - self.max_energy_level
            self.current_energy_level = self.max_energy_level
            print(f"Molecule at position {self.position} has moved to max energy level {self.max_energy_level}")
            return excess_energy  # excess energy to transfer

    def transfer_energy(self, energy):
        potential_energy_level = self.current_energy_level + energy
        if potential_energy_level <= self.max_energy_level:
            self.current_energy_level = potential_energy_level
            print(f"Molecule at position {self.position} has received energy, now at level {self.current_energy_level}")
        else:
            print("The molecule received too much energy! It cannot accept it.")

class System:
    def __init__(self, temperature):
        self.molecules = []
        self.temperature = temperature  # Temperature of the system, which influences the energy levels

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    def pulse_light(self, energy):
        for molecule in self.molecules:
            excess_energy = molecule.absorb_light(energy + self.temperature)
            if excess_energy > 0:
                # Find the nearest molecule to transfer energy to
                distances = [self.distance(molecule.position, other_molecule.position) for other_molecule in self.molecules if other_molecule != molecule]
                min_distance_index = distances.index(min(distances))
                # Add in an orientation factor which reduces the energy transferred
                orientation_factor = random.uniform(0.1, 1.0)
                self.molecules[min_distance_index].transfer_energy(excess_energy * orientation_factor)


    @staticmethod
    def distance(position1, position2):
        x1, y1, z1 = position1
        x2, y2, z2 = position2
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


system = System(temperature=2)
system.add_molecule(Molecule(5, (0, 0, 0), (1, 0, 0)))  # Molecule with max energy level 5 at position (0, 0, 0)
system.add_molecule(Molecule(5, (1, 1, 1), (0, 1, 0)))  # Molecule with max energy level 5 at position (1, 1, 1)
system.add_molecule(Molecule(5, (2, 2, 2), (0, 0, 1)))  # Molecule with max energy level 5 at position (2, 2, 2)

system.pulse_light(2)  # Gives 2 units of energy to all molecules
system.pulse_light(4)  # Gives 4 units of energy to all molecules