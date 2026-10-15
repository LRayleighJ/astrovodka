# Examples

Runnable examples, in two parallel trees.

| | English | 中文 |
|---|---|---|
| **Tutorial notebook** — the full walkthrough on a real TNG-100 subhalo | [`en/tutorial_TNG100-HI-mock-observations.ipynb`](en/tutorial_TNG100-HI-mock-observations.ipynb) | [`zh/tutorial_TNG100-HI-mock-observations.ipynb`](zh/tutorial_TNG100-HI-mock-observations.ipynb) |
| **Demo script** — an analytic galaxy, or a real subhalo with `--basePath` | [`en/astrovodka_demo.py`](en/astrovodka_demo.py) | [`zh/astrovodka_demo.py`](zh/astrovodka_demo.py) |
| **Batch script** — one cube per object, over a subhalo catalog | [`en/batch_subhalos.py`](en/batch_subhalos.py) | [`zh/batch_subhalos.py`](zh/batch_subhalos.py) |

The two tutorials share byte-identical code, so they produce identical numbers;
only the prose and the figure labels differ.

---

## Running them

### The demo script needs no simulation data

```bash
python examples/en/astrovodka_demo.py
```

This builds an analytic rotating disc at z = 0.01 and writes a FITS cube to
`examples/outputs/`. It prints a mass-conservation check that should read 100%.

### The TNG examples need a snapshot

Point `TNG_BASE_PATH` at the directory *containing* `groups_<snap>/` and
`snapdir_<snap>/`:

```bash
export TNG_BASE_PATH=/path/to/TNG/TNG-100/

python examples/en/astrovodka_demo.py --snapNum 99 --subID 113
python examples/en/batch_subhalos.py --limit 10
```

The environment variable is the only place a data path appears — nothing in the
repository hard-codes one.

### The notebooks

Run them in Jupyter or the VS Code notebook UI. They ship with their outputs
cleared, so you will see them fill in as they go.

---

## Outputs

Everything written by these examples goes to `outputs/`, which is git-ignored.
Nothing here writes into the repository proper, and nothing needs cleaning up
between runs.
