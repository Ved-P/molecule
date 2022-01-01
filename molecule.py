# Molecule
#
# This program takes in a molecular formula and creates a Lewis diagram and a 3D
# model of the molecule as the output.
#
# Author: Ved Pradhan
# Since:  December 31, 2021

import json
import matplotlib

# Opens the JSON file for use.
with open("elements.json", "r", encoding="utf8") as file:
    data = json.load(file)

# Gets the formula from the user.
formula = input("\n\n\nWelcome to Molecule! Please enter a molecular formula"
+ "(case sensitive): ")

# Class to represent each individual atom in the molecule.
class Atom:
    def __init__(self, symbol):
        self.symbol = symbol
        self.loose_ve = -1
        self.sigma_bonds = -1
        self.pi_bonds = -1
        self.formal_charge = -1
        self.total_ve = -1
        self.element = self.get_element()
        self.enegativity = self.element["electronegativity_pauling"]
        self.expected_ve = self.get_valence_electrons()

    # Returns the number of valence electrons the atom is expected to have.
    def get_valence_electrons(self):
        if self.symbol == "He":
            return 2
        elif 9 <= self.element["ypos"] <= 10:
            return 2
        elif 2 <= self.element["xpos"] <= 12:
            return 2
        else:
            return self.element["xpos"] % 10

    # Updates the formal charge of the atom.
    def update_formal_charge(self):
        self.formal_charge = self.expected_ve - self.loose_ve - self.sigma_bonds - self.pi_bonds

    # Updates the total number of valence electrons, including shared ones.
    def update_total_ve(self):
        self.total_ve = self.loose_ve + 2 * (self.sigma_bonds + self.pi_bonds)

# Retrieves the element corresponding to the given symbol.
def get_element(symbol):
    for element in data["elements"]:
        if element["symbol"] == symbol:
            return element
