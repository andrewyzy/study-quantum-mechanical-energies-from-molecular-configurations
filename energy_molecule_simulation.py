class Molecule:
    def __init__(self, max_energy_level, position):
        self.max_energy_level = max_energy_level
        self.current_energy_level = 0
        self.position = position  # position as a simple integer for now

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
    def __init__(self):
        self.molecules = []

    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    def pulse_light(self, energy):
        for molecule in self.molecules:
            excess_energy = molecule.absorb_light(energy)
            if excess_energy > 0:
                # Find a nearby molecule to transfer energy to
                # For simplicity, just find the next molecule in the list
                next_molecule_index = (self.molecules.index(molecule) + 1) % len(self.molecules)
                self.molecules[next_molecule_index].transfer_energy(excess_energy)


my_molecule = Molecule(5, 0)  # A molecule that can have energy levels 0 through 5 and is at position 0

# Molecule absorbs light and increases its energy level
my_molecule.absorb_light(2)
my_molecule.absorb_light(3)
my_molecule.absorb_light(1)

system = System()
system.add_molecule(Molecule(5, 0))  # Molecule with max energy level 5 at position 0
system.add_molecule(Molecule(5, 1))  # Molecule with max energy level 5 at position 1
system.add_molecule(Molecule(5, 2))  # Molecule with max energy level 5 at position 2

system.pulse_light(2)  # Gives 2 units of energy to all molecules
system.pulse_light(4)  # Gives 4 units of energy to all molecules
