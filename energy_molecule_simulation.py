import math
import random
import numpy as np

class Molecule:
    def __init__(self, max_energy_level, position, orientation, molecule_type):
        self.max_energy_level = max_energy_level
        self.current_energy_level = 0
        self.position = position  # position as a tuple (x, y, z)
        self.orientation = np.array(orientation)  # orientation as a vector [dx, dy, dz]
        self.molecule_type = molecule_type  # the type of molecule (e.g., 'methane', 'ethane', 'propane')
        self.h = 0.1  # A chosen energy unit. (akin to Planck's constant), and any energy gain will be rounded to the nearest multiple of h. Note that this is a gross oversimplification and doesn't fully capture quantum mechanical behavior

        # Depending on the type of molecule, set different properties
        if molecule_type == 'CH4':  # methane
            self.activation_energy = 50000
            self.molar_mass = 16.043  # in g/mol
            self.heat_capacity = 35.7  # in J/(K*mol)
            self.density = 0.657  # in kg/m3
        elif molecule_type == 'C2H6':  # ethane
            self.activation_energy = 60000
            self.molar_mass = 30.07
            self.heat_capacity = 52.7
            self.density = 1.356
        elif molecule_type == 'C3H8':  # propane
            self.activation_energy = 70000
            self.molar_mass = 44.1
            self.heat_capacity = 73.6
            self.density = 2.009
        else:
            raise ValueError('Invalid molecule type')

    def absorb_light(self, energy, temperature):
        # Apply the Arrhenius equation to adjust the max energy level based on the temperature
        R = 8.314  # Universal gas constant in J/(mol*K)
        Ea = self.activation_energy  # Activation energy in J/mol, specific to the type of molecule
        A = self.max_energy_level  # We'll just use the max_energy_level as the pre-exponential factor for simplicity
        adjusted_max_energy_level = A * math.exp(-Ea / (R * temperature))

        # Adjust the absorbed energy based on the heat capacity
        absorbed_energy = energy / self.heat_capacity  # This is a simplification - in reality the relationship may not be linear

        potential_energy_level = self.current_energy_level + absorbed_energy
        if potential_energy_level <= adjusted_max_energy_level:
            self.current_energy_level = potential_energy_level
            print(f"Molecule at position {self.position} has moved to energy level {self.current_energy_level}")
            return 0  # no excess energy
        else:
            excess_energy = potential_energy_level - self.max_energy_level
            self.current_energy_level = self.max_energy_level
            print(f"Molecule at position {self.position} has moved to max energy level {self.max_energy_level}")
            # Adjust the excess energy based on the number of molecules in a given volume
            volume = 1  # For simplicity, let's assume a volume of 1 m3
            number_of_molecules = self.density * volume / self.molar_mass  # This is a rough estimate - in reality we'd need to use the Ideal Gas Law or another model to calculate this
            adjusted_excess_energy = excess_energy / number_of_molecules  # Again, this is a simplification - in reality the relationship may not be linear
            return adjusted_excess_energy  # return the adjusted excess energy


    def absorb_light(self, energy, temperature):
        # Apply the Arrhenius equation to adjust the max energy level based on the temperature
        R = 8.314  # Universal gas constant in J/(mol*K)
        Ea = self.activation_energy  # Activation energy in J/mol, specific to the type of molecule
        A = self.max_energy_level  # We'll just use the max_energy_level as the pre-exponential factor for simplicity
        adjusted_max_energy_level = A * math.exp(-Ea / (R * temperature))

        # Adjust the absorbed energy based on the heat capacity
               # Adjust the absorbed energy based on the heat capacity
        absorbed_energy = energy / self.heat_capacity
        print(f"Absorbed energy: {absorbed_energy}")  # Add this line

        potential_energy_level = self.current_energy_level + absorbed_energy
        print(f"Potential energy level before rounding: {potential_energy_level}")  # Add this line

        # Round to the nearest multiple of h
        potential_energy_level = round(potential_energy_level / self.h) * self.h

        if potential_energy_level <= adjusted_max_energy_level:
            self.current_energy_level = potential_energy_level
            print(f"Molecule at position {self.position} has moved to energy level {self.current_energy_level}")
            return 0  # no excess energy
        else:
            excess_energy = potential_energy_level - self.max_energy_level
            self.current_energy_level = self.max_energy_level
            print(f"Molecule at position {self.position} has moved to max energy level {self.max_energy_level}")
            volume = 1  # For simplicity, let's assume a volume of 1 m3
            number_of_molecules = self.density * volume / self.molar_mass  # This is a rough estimate - in reality we'd need to use the Ideal Gas Law or another model to calculate this
            adjusted_excess_energy = excess_energy / number_of_molecules  # Again, this is a simplification - in reality the relationship may not be linear
            return adjusted_excess_energy  # return the adjusted excess energy


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
        self.temperature = temperature  # Temperature of the system, which influences the energy slevels

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    def pulse_light(self, energy):
        for molecule in self.molecules:
            excess_energy = molecule.absorb_light(energy + self.temperature, self.temperature)
            if excess_energy > 0:
                # Find the nearest molecule to transfer energy to
                distances = [self.distance(molecule.position, other_molecule.position) for other_molecule in self.molecules if other_molecule != molecule]
                min_distance_index = distances.index(min(distances))
                # Add in an orientation factor which reduces the energy transferred
                orientation_factor = self.calculate_orientation_factor(molecule.orientation, self.molecules[min_distance_index].orientation)
                self.molecules[min_distance_index].transfer_energy(excess_energy * orientation_factor)
    
    @staticmethod
    def calculate_orientation_factor(orientation1, orientation2):
        cosine_angle = np.dot(orientation1, orientation2) / (np.linalg.norm(orientation1) * np.linalg.norm(orientation2))
        return abs(cosine_angle)  # return the absolute value to avoid negative energy transfer

    @staticmethod
    def distance(position1, position2):
        x1, y1, z1 = position1
        x2, y2, z2 = position2
        return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

system = System(temperature=2)
system.add_molecule(Molecule(5, (0, 0, 0), [1, 0, 0], 'CH4'))  # A methane molecule
system.add_molecule(Molecule(5, (1, 1, 1), [0, 1, 0], 'C2H6'))  # An ethane molecule
system.add_molecule(Molecule(5, (2, 2, 2), [0, 0, 1], 'C3H8'))  # A propane molecule


system.pulse_light(2)  # Gives 2 units of energy to all molecules
system.pulse_light(4)  # Gives 4 units of energy to all molecules