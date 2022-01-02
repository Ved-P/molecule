# Molecule
#
# This program takes in a molecular formula and creates a Lewis diagram and a 3D
# model of the molecule as the output.
#
# Author: Ved Pradhan
# Since:  December 31, 2021

import json
import matplotlib.pyplot as plt
import sys
import math

# Opens the JSON file for use.
with open("elements.json", "r", encoding="utf8") as file:
    data = json.load(file)

# Gets the formula and charge from the user.
formula = input("\n\n\nWelcome to Molecule! Please enter a molecular formula "
+ "(case sensitive): ")
temp = input("What is the charge of the molecule? Enter an integer (0 for no "
+ "charge): ")
try:
    charge = int(temp)
except ValueError:
    print("Error: '" + temp + "' is not a valid charge.\n\n\n")
    sys.exit()

# A list to store each individual atom in the molecule.
atoms = []

# A dictionary to store each type of element and its frequency.
element_frequency = {}

# A list to store the bonds between Atom objects.
bonds = []

# Class to represent each individual atom in the molecule.
class Atom:
    def __init__(self, symbol):
        self.symbol = symbol
        self.element = get_element(symbol)
        if self.element != False:
            self.enegativity = self.element["electronegativity_pauling"]
            self.expected_ve = self.get_valence_electrons()
        self.loose_ve = 0
        self.sigma_bonds = 0
        self.pi_bonds = 0
        self.formal_charge = 0
        self.total_ve = 0
        self.lewis_x = 0
        self.lewis_y = 0

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

    # Returns essential information about the atom as a string.
    def __str__(self):
        return (self.element["name"] + ": " + str(self.loose_ve) + " loose, "
        + str(self.sigma_bonds) + " sigma, " + str(self.pi_bonds) + " pi")

# Retrieves the element corresponding to the given symbol.
def get_element(symbol):
    for element in data["elements"]:
        if element["symbol"] == symbol:
            return element
    print("Error: Element '" + symbol + "' not found.\n\n\n")
    return False

# Parses through the inputted formula, splitting it into elements and frequencies.
def parse(form):
    i = 1
    while i < len(form) and not(ord('A') <= ord(form[i]) <= ord('Z')):
        i += 1
    j = i - 1
    while j >= 0 and ord('0') <= ord(form[j]) <= ord('9'):
        j -= 1
    if j < 0:
        print("Error: The formula cannot start with a number.\n\n\n")
        sys.exit()
    symbol_part = form[:j+1]
    number_part = form[j+1:i]
    rest = form[i:]
    ele = get_element(symbol_part)
    if number_part == "":
        number = 1
    else:
        number = int(number_part)
    element_frequency[symbol_part] = number
    for i in range(number):
        atoms.append(Atom(symbol_part))
    if len(rest) > 0:
        parse(rest)

# Prints a "not supported" message and quits the program.
def noSupport():
    print("Sorry, this molecule is not supported yet.\n\n\n")
    sys.exit()

# Checks if the molecule is supported.
def check():
    if len(element_frequency) != 2:
        noSupport()
    symb1 = list(element_frequency)[0]
    symb2 = list(element_frequency)[1]
    global center
    global outer
    if symb1 == "H":
        center = symb2
        outer = symb1
    elif symb2 == "H":
        center = symb1
        outer = symb2
    elif get_element(symb1)["electronegativity_pauling"] < get_element(symb2)["electronegativity_pauling"]:
        center = symb1
        outer = symb2
    elif get_element(symb1)["electronegativity_pauling"] > get_element(symb2)["electronegativity_pauling"]:
        center = symb2
        outer = symb1
    else:
        noSupport()
    if element_frequency[center] != 1:
        noSupport()

# Bonds two atoms together; updates in the object and the data structure.
def bond(atom1, atom2, type):
    bonds.append((atom1, atom2, type))
    if (type == "sigma"):
        atom1.sigma_bonds += 1
        atom2.sigma_bonds += 1
    if (type == "pi"):
        atom1.pi_bonds += 1
        atom2.pi_bonds += 1

# Distributes the valence electrons as loose ones or through bonds.
def distribute():
    total_ve = 0
    for a in atoms:
        total_ve += a.expected_ve
    total_ve -= charge
    left_ve = total_ve
    global centerAtom
    centerAtom = -1
    global outerAtoms
    outerAtoms = []
    for a in atoms:
        if a.symbol == center:
            centerAtom = a
        elif a.symbol == outer:
            outerAtoms.append(a)
    for o in outerAtoms:
        bond(centerAtom, o, "sigma")
        left_ve -= 2
    want_ve = -1
    if outer == "H" or outer == "He":
        want_ve = 0
    else:
        want_ve = 6
    if left_ve // len(outerAtoms) >= want_ve:
        for o in outerAtoms:
            o.loose_ve += want_ve
            left_ve -= want_ve
        if left_ve >= 0:
            centerAtom.loose_ve += left_ve
    else:
        noSupport()

# Draws the lewis diagram using matplotlib.
def draw_lewis():
    centerAtom.lewis_x = 0
    centerAtom.lewis_y = 0
    plt.style.use('_mpl-gallery')
    fig, ax = plt.subplots()
    fig.suptitle(formula, fontsize=14, fontweight='bold')
    ax.text(0, 0, centerAtom.symbol, verticalalignment='center', horizontalalignment='center')
    for i in range(len(outerAtoms)):
        o = outerAtoms[i]
        o.lewis_x = math.cos(2 * i * math.pi / len(outerAtoms))
        o.lewis_y = math.sin(2 * i * math.pi / len(outerAtoms))
        ax.text(o.lewis_x, o.lewis_y, o.symbol, verticalalignment='center', horizontalalignment='center')
    for b in bonds:
        x1 = (2 * b[0].lewis_x + b[1].lewis_x) / 3
        x2 = (b[0].lewis_x + 2 * b[1].lewis_x) / 3
        y1 = (2 * b[0].lewis_y + b[1].lewis_y) / 3
        y2 = (b[0].lewis_y + 2 * b[1].lewis_y) / 3
        plt.plot([x1, x2], [y1, y2], color='gray')
    for a in atoms:
        x_shift = 0
        y_shift = 0
        for i in range(a.loose_ve):
            if 0 <= i <= 1:
                x_shift = -0.2
            elif 2 <= i <= 3:
                y_shift = -0.2
            elif 4 <= i <= 5:
                x_shift = 0.2
            elif 6 <= i <= 7:
                y_shift = 0.2
            if i == 0 or i == 5:
                y_shift = 0.05
            elif i == 1 or i == 4:
                y_shift = -0.05
            elif i == 2 or i == 7:
                x_shift = -0.05
            elif i == 3 or i == 6:
                x_shift = 0.05
            ax.scatter(x = a.lewis_x + x_shift, y = a.lewis_y + y_shift + 0.03,
            s = 4, color='black')
    axes = plt.gca()
    axes.set_aspect(1)
    plt.xlim([-1.75, 1.75])
    plt.ylim([-1.7, 1.8])
    axes.axes.xaxis.set_visible(False)
    axes.axes.yaxis.set_visible(False)
    plt.show()

parse(formula)
check()
distribute()
print(element_frequency)
for a in atoms:
    print(a)
draw_lewis()
print("\n\n\n")
