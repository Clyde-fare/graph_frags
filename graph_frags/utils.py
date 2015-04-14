__author__ = 'clyde'

unwind = lambda l: (e for sl in l for e in sl)


def get_unique(frags):
    """Get unique fragments"""
    atoms = lambda e: set(unwind(e))
    list_atoms = []
    for f in frags:
        list_atoms.append(atoms(f))

    unique_atoms = []
    for a in list_atoms:
        if a not in unique_atoms:
            unique_atoms.append(a)

    return unique_atoms

