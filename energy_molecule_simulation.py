import math
import random
import numpy as np


class Molecule:
    def __init__(self, max_energy_level, position, orientation, molecule_type, spectral_absorption_coefficients):
        self.max_energy_level = max_energy_level
        self.current_energy_level = 0
        self.position = position  # position as a tuple (x, y, z)
        self.orientation = np.array(orientation)  # orientation as a vector [dx, dy, dz]
        self.molecule_type = molecule_type  # the type of molecule (e.g., 'methane', 'ethane', 'propane')
        self.h = 0.1  # A chosen energy unit. (akin to Planck's constant)
        self.spectral_absorption_coefficients = spectral_absorption_coefficients

        self.n = 0  # Principle quantum number

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

    def absorb_light(self, energy, wavelength, temperature):
        R = 8.314
        Ea = self.activation_energy
        A = self.max_energy_level
        adjusted_max_energy_level = A * math.exp(-Ea / (100 * R * temperature))
        print("adjusted_max_energy_level = ", adjusted_max_energy_level)

        if wavelength in self.spectral_absorption_coefficients:
            energy *= self.spectral_absorption_coefficients[wavelength]
        
        # Add a random factor to energy absorption
        absorbed_energy = random.uniform(0, energy / self.heat_capacity)
        
        # Quantum mechanical behavior
        energy_level = round(absorbed_energy / self.h)
        if energy_level - self.n == 1:  # Check the selection rule Δn = 1
            self.n = energy_level
            absorbed_energy = energy_level * self.h
        else:
            absorbed_energy = 0  # If the energy doesn't match a quantum leap, it's not absorbed

        print(f"Absorbed energy: {absorbed_energy}")

        potential_energy_level = self.current_energy_level + absorbed_energy
        print(f"Potential energy level before rounding: {potential_energy_level}")

        # Add a condition to check if the potential energy level is too high
        if potential_energy_level > adjusted_max_energy_level:
            absorbed_energy = adjusted_max_energy_level - self.current_energy_level
            potential_energy_level = self.current_energy_level + absorbed_energy
            return 0
        else:
            excess_energy = potential_energy_level - adjusted_max_energy_level
            self.current_energy_level = adjusted_max_energy_level
            print(f"Molecule at position {self.position} has moved to max energy level {adjusted_max_energy_level}")

            volume = 1
            number_of_molecules = self.density * volume / self.molar_mass
            adjusted_excess_energy = excess_energy / number_of_molecules
            return adjusted_excess_energy



    def transfer_energy(self, energy):
        potential_energy_level = self.current_energy_level + energy
        print("potential_energy_level = ", potential_energy_level)

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

    def pulse_light(self, energy, wavelength):
        for molecule in self.molecules:
            # Distributed energy is calculated based on quantum energy levels
            distributed_energy = round(random.uniform(0, energy) / molecule.h) * molecule.h
            excess_energy = molecule.absorb_light(distributed_energy + self.temperature, wavelength, self.temperature)
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
    

system = System(temperature=200)
spectral_absorption_coefficients = {400: 0.8, 500: 0.9, 600: 1.0, 700: 0.9, 800: 0.8}
system.add_molecule(Molecule(5, (0, 0, 0), [1, 0, 0], 'CH4', spectral_absorption_coefficients))
system.add_molecule(Molecule(5, (1, 1, 1), [0, 1, 0], 'C2H6', spectral_absorption_coefficients))
system.add_molecule(Molecule(5, (2, 2, 2), [0, 0, 1], 'C3H8', spectral_absorption_coefficients))

# Now you also need to specify the wavelength when pulsing light
system.pulse_light(2, 500)
system.pulse_light(4, 600)