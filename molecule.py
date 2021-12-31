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

# Given the symbol of the element, return the number of its valence electrons.
def get_valence_electrons(symbol):
    for element in data["elements"]:
        if element["symbol"] == symbol:
            if symbol == "He":
                return 2
            elif 9 <= element["ypos"] <= 10:
                return 2
            elif 2 <= element["xpos"] <= 12:
                return 2
            else:
                return element["xpos"] % 10

# Gets the formula from the user.
formula = input("\n\n\nWelcome to Molecule! Please enter a molecular formula"
    + "(case sensitive): ")
