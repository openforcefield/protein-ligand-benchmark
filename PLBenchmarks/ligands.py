"""
ligands.py
Functions and classes for handling the ligand data.
"""

from PLBenchmarks import utils, targets

import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
from openforcefield.topology import Molecule

from PIL import Image

import io

import yaml

try:
    from importlib.resources import open_text
except ImportError:
    # Python 2.x backport
    from importlib_resources import open_text


class ligand:
    """
    Store and convert the data of one ligand in a :py:class:`pandas.Series`.

    """

    _observables = ["dg", "dh", "tds", "ki", "ic50", "pic50"]

    def __init__(self, d: dict, target: str = None):
        """
        Initialize :py:class:`PLBenchmarks.ligands.ligand` object from :py:class:`dict` and store data in a :py:class:`pandas.Series`.

        :param d: :py:class:`dict` with the ligand data
        :return None
        """
        self._target = target
        self._data = pd.Series(d)
        self._molecule = None
        self._name = self._data["name"][0]
        # Expand measurement dict into single fields
        if "measurement" in list(self._data.index):
            meas = pd.Series(self._data["measurement"])
            meas.index = ["measurement:" + c for c in meas.index]
            self._data.drop(["measurement"], inplace=True)
            self._data = pd.concat([self._data, meas])
            index = self._data.index.to_series().str.split(":", expand=True).fillna("")
            self._data.index = pd.MultiIndex.from_arrays(
                [index[0].tolist(), index[1].tolist()]
            )
            for obs in self._observables:
                if ("measurement", obs) in list(self._data.index):
                    self._data = self._data.append(
                        pd.Series(
                            [0],
                            index=pd.MultiIndex.from_tuples(
                                [("measurement", f"e_{obs}")]
                            ),
                        )
                    )
                    vals = self._data[("measurement", f"{obs}")]
                    u = utils.ureg("")
                    if vals[2] == "nM":
                        u = utils.ureg("nanomolar")
                    elif vals[2] == "uM":
                        u = utils.ureg("micromolar")
                    elif vals[2] == "kj/mol":
                        u = utils.ureg("kJ / mole")
                    elif vals[2] == "kcal/mol":
                        u = utils.ureg("kcal / mole")
                    else:
                        # let pint figure out what the unit means
                        u = utils.ureg(vals[2])
                    self._data[("measurement", f"e_{obs}")] = vals[1] * u
                    self._data[("measurement", obs)] = vals[0] * u

    def deriveObservables(
        self,
        derivedObs="dg",
        dest="DerivedMeasurement",
        outUnit=utils.ureg("kcal / mole"),
    ):
        """
        Derive observables from (stored) primary data, which is then stored in the :py:class:`pandas.DataFrame`

        :param derivedObs: type of derived observable, can be any of 'dg' 'ki', 'ic50' or 'pic50'
        :param dest: string with column name for 'pandas.DataFrame' where the derived observable should be stored.
        :param outUnit: unit of type :py:class:`pint` unit of derived coordinate
        :return: None
        """
        assert (
            derivedObs in self._observables
        ), "Observable to be derived not known. Should be any of dg, ki, ic50, or pic50"
        for obs in self._observables:
            if ("measurement", obs) in list(self._data.index):
                self._data = self._data.append(
                    pd.Series(
                        [
                            utils.convertValue(
                                self._data[("measurement", obs)],
                                obs,
                                derivedObs,
                                outUnit=outUnit,
                            ),
                            utils.convertError(
                                self._data[("measurement", f"e_{obs}")],
                                self._data[("measurement", obs)],
                                obs,
                                derivedObs,
                                outUnit=outUnit,
                            ),
                        ],
                        index=pd.MultiIndex.from_tuples(
                            [(dest, derivedObs), (dest, f"e_{derivedObs}")]
                        ),
                    )
                )
                break
        else:
            print(f"Conversion to observable {derivedObs} not possible.")

    def getName(self):
        """
        Access the name of the ligand.

        :return: name: string
        """
        return self._data["name"][0]

    def getDF(self, cols=None):
        """
        Access the ligand data as a :py:class:`pandas.Dataframe`

        :param cols: list of columns which should be returned in the :py:class:`pandas.Dataframe`
        :return: :py:class:`pandas.Dataframe`
        """
        if cols:
            return self._data[cols]
        else:
            return self._data

    def findLinks(self):
        """
        Processes primary data to have links in the html string of the ligand data

        :return: None
        """
        if ("measurement", "doi") in list(self._data.index):
            doi = self._data["measurement", "doi"]
            if str(doi) != "nan":
                res = []
                for ddoi in re.split(r"[; ]+", str(doi)):
                    res.append(utils.findDoiUrl(ddoi))
            self._data["measurement", "doi_html"] = (r"\n").join(res)
            self._data.drop([("measurement", "doi")], inplace=True)
            self._data.rename({"doi_html": "Reference"}, level=1, inplace=True)
        if ("pdb") in list(self._data.index):
            pdb = self._data["pdb"]
            self._data["pdb_html"] = utils.findPdbUrl(pdb)
            self._data.drop(["pdb"], inplace=True)
            self._data.rename({"pdb_html": "pdb"}, inplace=True)

    def getCoordFilePath(self):
        """
        Get file path relative to the PLBenchmarks repository of the SDF coordinate file of the docked ligand

        :return: file path as string
        """
        fname = self._data["docked"][0]
        return fname

    def getMol(self):
        """
        Get molecule object with coordinates of the docked ligand

        :return: file path as string
        """
        if self._molecule is not None:
            fname = self.getCoordFilePath()

            for td in target_list:
                if td["name"] == self._target:
                    sdfpath = ["PLBenchmarks", "data", td["dir"]] + fname.split("/")
                    break
            else:
                raise ValueError(f"Path for ligand {self._name} not found.")

            file = open_text(".".join(sdfpath[:-1]), sdfpath[-1])
            self._molecule = Molecule.from_file(file, "sdf")
        return self._molecule

    def addMolToFrame(self):
        """
        Adds a image file of the ligand to the :py:class:`pandas.Dataframe`

        :return: None
        """
        PandasTools.AddMoleculeColumnToFrame(
            self._data, smilesCol="smiles", molCol="ROMol", includeFingerprints=False
        )
        self._data["ROMol"].apply(lambda x: x[0])

    def getHTML(self, columns=None):
        """
        Access the ligand as a HTML string

        :param columns: list of columns which should be returned in the :py:class:`pandas.Dataframe`
        :return: HTML string
        """
        self.findLinks()
        if columns:
            html = pd.Dataframe(self._data[columns]).to_html()
        else:
            html = pd.Dataframe(self._data).to_html()
        html = html.replace("REP1", '<a target="_blank" href="')
        html = html.replace("REP2", '">')
        html = html.replace("REP3", "</a>")
        html = html.replace("\\n", "<br>")
        return html

    def getImg(self):
        """
        Creates a molecule image.

        :return: :py:class:`PIL.Image` object
        """
        dr = Draw.MolDraw2DCairo(200, 200)
        opts = dr.drawOptions()
        opts.clearBackground = True
        mol = Chem.MolFromSmiles(self._data["smiles"][0])
        Chem.rdDepictor.Compute2DCoords(mol)
        dr.DrawMolecule(mol, legend=self._data["name"][0])
        img = Image.open(io.BytesIO(dr.GetDrawingText())).convert("RGBA")
        datas = img.getdata()

        newData = []
        for item in datas:
            if item[0] == 255 and item[1] == 255 and item[2] == 255:
                newData.append((255, 255, 255, 0))
            else:
                newData.append(item)
        img.putdata(newData)
        return img


class ligandSet(dict):
    """
    Class inherited from dict to store all available ligands of one target.

    """

    def __init__(self, target, *arg, **kw):
        """
        Initializes :py:class:`~PLBenchmarks.ligands.ligandSet` class

        :param target: string name of target
        :param arg: arguments for :py:class:`dict` (base class)
        :param kw: keywords for :py:class:`dict` (base class)
        """
        super(ligandSet, self).__init__(*arg, **kw)
        tp = targets.getTargetDataPath(target)
        file = open_text(".".join(tp), "ligands.yml")
        data = yaml.full_load_all(file)
        for d in data:
            l = ligand(d, target)
            l.deriveObservables(derivedObs="dg")
            # l.findLinks()
            l.addMolToFrame()
            self[l.getName()] = l
        file.close()

    def getLigand(self, name):
        """
        Accesses one ligand of the :py:class:`~PLBenchmarks:ligands.ligandSet`

        :param name: string name of the ligand
        :return: :py:class:`PLBenchmarks.ligands.ligand` class
        """
        for key in self.keys():
            if key == name:
                return self[key]
                break
        else:
            raise ValueError(f"Ligand {name} is not part of set.")

    def getDF(self, columns=None):
        """
        Access the :py:class:`~PLBenchmarks:ligands.ligandSet` as a :py:class:`pandas.Dataframe`

        :param columns: :py:class:`list` of columns which should be returned in the :py:class:`pandas.Dataframe`
        :return: :py:class:`pandas.Dataframe`
        """
        dfs = []
        for key, item in self.items():
            dfs.append(item.getDF(columns))
        df = pd.concat(dfs, axis=1).T
        return df

    def getHTML(self, columns=None):
        """
        Access the :py:class:`PLBenchmarks:ligands.ligandSet` as a HTML string

        :param cols: :py:class:`list` of columns which should be returned in the :py:class:`pandas.Dataframe`
        :return: HTML string
        """
        for key, item in self.items():
            item.findLinks()
        df = self.getDF(columns)
        html = df.to_html()
        html = html.replace("REP1", '<a target="_blank" href="')
        html = html.replace("REP2", '">')
        html = html.replace("REP3", "</a>")
        html = html.replace("\\n", "<br>")
        return html
