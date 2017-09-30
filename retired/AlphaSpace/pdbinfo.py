#! /usr/bin/python

strip_list = []


##strip_list = ["Na+","Cl-"]

def isATOM(line):
    """
	Check if line is ATOM, excluding HETATM 
	"""
    if line[0:6] == "ATOM  ":
        return 1
    else:
        return 0


def isHETATM(line):
    """
        Check if line is HETATM, excluding ATOM
        """
    if line[0:6] == "HETATM":
        return 1
    else:
        return 0


def isAtom(line):
    """
	Check the line is atom line or not
	"""
    if line[0:6] == "ATOM  " or line[0:6] == "HETATM":
        return 1
    else:
        return 0


def atmi(line):
    """
	Return atom id
	"""
    return line[6:11]


def atmn(line):
    """
	Return atom name
	"""
    return line[12:16]


def roti(line):
    """
	Return rotimer ID
	"""
    return line[16:17]


def resn(line):
    """
	Return residue name
	"""
    return line[17:20]


def chid(line):
    """
	Return chain id
	"""
    return line[21:22]


def resi(line):
    """
	Return residue id
	"""
    return line[22:26]


def seqi(line):
    """"
	Combine chain id and residue id, in case of residue id is larger than 10000 
	which make residue id not unique
	"""
    return line[21:26]


def coord(line):
    """
	Return atom atom_xyz in format [float(x), float(y), float(z)]
	"""
    crd = [float(line[30 + 8 * i: 38 + 8 * i]) for i in range(3)]
    return crd


def isHydrogen(line):
    """
	Check is hydrogen atom or not
	"""
    if line[12] == 'H' or line[13] == 'H':
        return 1
    else:
        return 0


def isProtein(line):
    """
	Return true if not water
	"""
    if resn(line) not in strip_list and resn(line) != "WAT":
        return 1
    else:
        return 0


def isWater(line):
    """
	Return true if is water
	"""
    if resn(line) == "WAT" or resn(line) == "HOH":
        return 1
    else:
        return 0


def isNa(line):
    """
    Return true if is sodium ion
    """
    if resn(line) == "Na+" or resn(line) == "NA+" or resn(line) == "NA" or resn(line) == "Na":
        return 1
    else:
        return 0


def isCl(line):
    """
    Return true if is chlorine ion
    """
    if resn(line) == "Cl-" or resn(line) == "CL-" or resn(line) == "CL" or resn(line) == "Cl":
        return 1
    else:
        return 0


if __name__ == "__main__":
    print("  This module is for abstracting desired information from PDB file. ")
    print("  NO calculation in this module.")
