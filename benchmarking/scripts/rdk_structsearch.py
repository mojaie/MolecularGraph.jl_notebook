from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def drawsvg(mol):
    drawer = rdMolDraw2D.MolDraw2DSVG(400,400)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:','')
    ISVG(svg)


molsiter = Chem.SDMolSupplier("../_resources/open structures.sdf")

for mol in molsiter:
    if mol is None:
        continue
    drawsvg(mol)