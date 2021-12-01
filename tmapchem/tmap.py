import pickle
import numpy as np
import tmap as tm
import pandas as pd
import scipy.stats as ss
from rdkit.Chem import AllChem
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from collections import Counter
from matplotlib.colors import ListedColormap
from matplotlib import pyplot as plt
import matplotlib
import shutil
import os
import uuid
from matplotlib.cm import get_cmap


class Tmap(object):

    def __init__(self, input_file, output_folder, name=None):
        self.output_folder = os.path.abspath(output_folder)
        os.makedirs(self.output_folder, exist_ok=True)
        self.input_file = os.path.abspath(input_file)
        self.df = pd.read_csv(input_file)
        if name is None:
            self.name = os.path.basename(input_file).split(".")[0]
        self.cwd = os.getcwd()

    def process_data(self):
        self._remove_nans()
        self._to_ranks()

    def layout(self):
        os.chdir(self.output_folder)

        df = self.df

        enc = MHFPEncoder(1024)
        lf = tm.LSHForest(1024, 64)

        lf_file = os.path.join(self.output_folder, "lf.dat")
        props_file = os.path.join(self.output_folder, "props.pickle")

        if os.path.exists(lf_file):
            lf.restore(lf_file)
            hac, c_frac, ring_atom_frac, largest_ring_size = pickle.load(
             open(props_file, "rb")
            )

        else:
            fps = []
            hac = []
            c_frac = []
            ring_atom_frac = []
            largest_ring_size = []

            for i, row in df.iterrows():
                if i != 0 and i % 1000 == 0:
                    print(100 * i / len(df))
                mol = AllChem.MolFromSmiles(row["Smiles"])
                atoms = mol.GetAtoms()
                size = mol.GetNumHeavyAtoms()
                n_c = 0
                n_ring_atoms = 0
                for atom in atoms:
                    if atom.IsInRing():
                        n_ring_atoms += 1
                    if atom.GetSymbol().lower() == "c":
                        n_c += 1

                c_frac.append(n_c / size)
                ring_atom_frac.append(n_ring_atoms / size)
                sssr = AllChem.GetSymmSSSR(mol)
                if len(sssr) > 0:
                    largest_ring_size.append(max([len(s) for s in sssr]))
                else:
                    largest_ring_size.append(0)
                hac.append(size)
                fps.append(tm.VectorUint(enc.encode_mol(mol)))

            lf.batch_add(fps)
            lf.index()

            lf.store(lf_file)
            with open(props_file, "wb+") as f:
                pickle.dump(
                    (hac, c_frac, ring_atom_frac, largest_ring_size),
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )

        c_frak_ranked = ss.rankdata(np.array(c_frac) / max(c_frac)) / len(c_frac)

        cfg = tm.LayoutConfiguration()
        cfg.node_size = 1/26
        cfg.mmm_repeats = 2
        # cfg.sl_repeats = 2 # drugbank
        cfg.sl_extra_scaling_steps = 5
        cfg.k = 20
        cfg.sl_scaling_type = tm.RelativeToAvgLength
        x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)

        os.chdir(self.cwd)

        return x, y, s, t

    def _is_numeric(self, col):
        for r in self.df[self.df[col].notnull()][col].tolist():
            try:
                float(r)
            except:
                return False
        return True

    def _is_categorical(self, col):
        return not self._is_numeric(col)

    def _to_ranks(self):
        for col in list(self.df.columns):
            if self._is_numeric(col):
                vals = self.df[self.df[col].notnull()][col]
                rnks = ss.rankdata(vals)
                self.df.loc[self.df[col].notnull(), col] = rnks

    def _remove_nans(self):
        self.df = self.df.dropna(how='all', axis=1)

    def _get_listed_colormap(self, values):
        n = len(set(values))
        cmap = get_cmap("tab20", n)
        colors = []
        for i in range(cmap.N):
            rgba = cmap(i)
            colors += [matplotlib.colors.rgb2hex(rgba)]
        lcm = ListedColormap(
            colors,
            name = str(uuid.uuid4())
        )
        return lcm

    def _fill_continuous_nans(self, values):
        values_ = []
        for v in values:
            if np.isnan(v):
                values_ += [-1]
            else:
                values_ += [v]
        return values_

    def _get_continuous_colormap(self, values):
        values_ = []
        has_nan = False
        for v in values:
            if np.isnan(v):
                values_ += [-1]
                has_nan = True
            else:
                values_ += [v]
        values = values_
        cm = plt.get_cmap("viridis").copy()
        if has_nan:
            cm._init()
            cm._lut[0,: ] = (0,0,0,1)
        return cm

    def display(self):

        os.chdir(self.output_folder)

        name = self.name
        df = self.df

        labels = []
        for r in df[["CompoundId", "Smiles", "Source"]].values:
            labels += [r[1]+"__"+r[2]+"__"+r[0]]

        selected_labels = ["Source", "CompoundId"]

        columns = [c for c in list(self.df.columns) if c not in ["CompoundId", "Smiles"]]

        data = []
        categorical = []
        colormap = []
        legend_labels = []
        max_legend_label = []
        min_legend_label = []
        for c in columns:
            if self._is_categorical(c):
                vals = self.df[c].tolist()
                labels_groups, groups = Faerun.create_categories(vals)
                legend_labels += [labels_groups]
                data += [(labels_groups, groups)]
                categorical += [True]
                colormap += [self._get_listed_colormap(vals)]
                max_legend_label += [None]
                min_legend_label += [None]
            else:
                vals = self.df[c].tolist()
                categorical += [False]
                colormap += [self._get_continuous_colormap(vals)]
                data += [self._fill_continuous_nans(vals)]
                max_legend_label += [1]
                min_legend_label += [0]

        c_data = []
        for d in data:
            if type(d) is tuple:
                c_data += [d[1]]
            else:
                c_data += [d]

        x, y, s, t = self.layout()

        f = Faerun(view="front", coords=False)

        f.add_scatter(
            name,
            {
                "x": x,
                "y": y,
                "c": c_data,
                "labels": labels,
            },
            shader="smoothCircle",
            colormap=colormap,
            point_scale=2.5,
            categorical=categorical,
            has_legend=True,
            legend_labels=legend_labels,
            selected_labels=selected_labels,
            series_title=columns,
            max_legend_label=max_legend_label,
            min_legend_label=min_legend_label,
            title_index=2,
            legend_title="",
        )

        f.add_tree(name+"_tree", {"from": s, "to": t}, point_helper=name)
        f.plot(name, template="smiles")

        shutil.move(name+".html", os.path.join(self.output_folder, name+".html"))
        shutil.move(name+".js", os.path.join(self.output_folder, name+".js"))

        os.chdir(self.cwd)

    def run(self):
        self.process_data()
        self.display()
