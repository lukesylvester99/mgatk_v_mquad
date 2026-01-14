import os
import numpy as np
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


##############
# Load data #
##############

VIREO_OUT = "vireo"
os.makedirs(VIREO_OUT, exist_ok=True)

print("")
print("loading AD mat")
AD = mmread("cellSNP.tag.AD.mtx").tocsc()

print("loading DP mat")
DP = mmread("cellSNP.tag.DP.mtx").tocsc()

print("")
print("loading variants.tsv")
# cellSNP.variants.tsv is typically columns like: CHROM POS REF ALT (and possibly extra columns)
# We will build IDs like: 14783T>C  (matches your mquad passed_variant_names.txt style)
var_tab = np.genfromtxt("cellSNP.variants.tsv", dtype=str)

# genfromtxt returns 1D if single row; normalize to 2D
if var_tab.ndim == 1:
    var_tab = var_tab.reshape(1, -1)

if var_tab.shape[1] < 4:
    raise ValueError(f"cellSNP.variants.tsv has {var_tab.shape[1]} columns; expected >= 4")

pos = var_tab[:, 1]
ref = var_tab[:, 2]
alt = var_tab[:, 3]
mtSNP_ids = np.array([f"{p}{r}>{a}" for p, r, a in zip(pos, ref, alt)], dtype=str)

print("")
print("SECTION: Loading data...")
print("AD shape:", AD.shape, "DP shape:", DP.shape)
print("n_variants:", AD.shape[0], "n_cells:", AD.shape[1])
print("variant_names:", len(mtSNP_ids))

# sanity checks
assert AD.shape == DP.shape, "AD and DP shapes differ"
assert len(mtSNP_ids) == AD.shape[0], "variants.tsv-derived IDs length != number of matrix rows"


#######################
# Fixed clone labels  #
#######################

CLONE_ASSIGN_TSV = "clone_assignments.tsv"
if not os.path.exists(CLONE_ASSIGN_TSV):
    raise FileNotFoundError(
        f"Missing {CLONE_ASSIGN_TSV}. Put the existing vireo clone assignments here."
    )

SAMPLES_TSV = "cellSNP.samples.tsv"
cell_barcodes = np.genfromtxt(SAMPLES_TSV, dtype=str)

assert len(cell_barcodes) == AD.shape[1], \
    "Number of barcodes does not match number of cells in AD/DP"

# Load clone assignments and map to current barcode order
assign = {}
with open(CLONE_ASSIGN_TSV, "r") as f:
    header = f.readline().rstrip("\n").split("\t")
    if len(header) < 2 or header[0] != "cell_barcode" or header[1] != "clone_id":
        raise ValueError(f"Unexpected header in {CLONE_ASSIGN_TSV}: {header}")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        bc, cid = parts[0], parts[1]
        assign[bc] = int(cid)

clone_ids = np.array([assign.get(bc, -1) for bc in cell_barcodes], dtype=int)

n_missing = int(np.sum(clone_ids < 0))
if n_missing > 0:
    missing_examples = cell_barcodes[clone_ids < 0][:10]
    raise ValueError(
        f"{n_missing} cells in cellSNP.samples.tsv are missing from {CLONE_ASSIGN_TSV}. "
        f"Examples: {missing_examples}"
    )

n_clones = int(clone_ids.max() + 1)
print("")
print("SECTION: Loaded fixed clone assignments")
print("n_clones:", n_clones)
print("Clone sizes:")
clone_counts = np.bincount(clone_ids, minlength=n_clones)
for k in range(n_clones):
    print(f"  clone {k}: {int(clone_counts[k])}")


######################################
# Mean AF (beta_mu) from fixed clones#
######################################

print("")
print("SECTION: Computing beta_mu (mean AF) per clone from fixed assignments...")

beta_mu = np.zeros((AD.shape[0], n_clones), dtype=float)

for k in range(n_clones):
    idx = np.where(clone_ids == k)[0]
    if len(idx) == 0:
        continue
    ad_sum = np.asarray(AD[:, idx].sum(axis=1)).ravel()
    dp_sum = np.asarray(DP[:, idx].sum(axis=1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        beta_mu[:, k] = np.divide(
            ad_sum, dp_sum,
            out=np.zeros_like(ad_sum, dtype=float),
            where=(dp_sum > 0)
        )

print("beta_mu shape:", beta_mu.shape)


######################################
# Mean AF Plot
######################################

print("")
print("SECTION: Generating mean AF plot...")

raw_col = cm.get_cmap('pink_r', 200)
new_col = np.vstack((raw_col(np.linspace(0, 0.7, 10)),
                     raw_col(np.linspace(0.7, 1, 90))))
segpink = ListedColormap(new_col, name='segpink')

fig = plt.figure(figsize=(5, 4), dpi=300)

im = plt.imshow(
    beta_mu,
    aspect="auto",
    cmap=segpink,
    interpolation="nearest",
    vmin=np.nanmin(beta_mu),
    vmax=np.nanmax(beta_mu)
)
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.title("Mean allelic ratio")
plt.xlabel("Clone")
plt.ylabel("%d SNPs" % (AD.shape[0]))
plt.xticks(range(n_clones))

plt.tight_layout()
plt.savefig(f"{VIREO_OUT}/mean_AF.png", dpi=300)
plt.close(fig)

print("")
print("Completed mean AF plot (mean_AF.png).")


#######################
# Heatmap
#######################

print("")
print("SECTION: Starting clone visualization and tsv file generation...")

AD_dense = AD.toarray()
DP_dense = DP.toarray()
AF = np.divide(
    AD_dense,
    DP_dense,
    out=np.zeros_like(AD_dense, dtype=float),
    where=(DP_dense > 0)
)

cell_label = np.array([f"clone{c}" for c in clone_ids])
id_uniq = [f"clone{i}" for i in range(n_clones)]

order = np.concatenate([np.where(cell_label == cid)[0] for cid in id_uniq])
AF_sorted = AF[:, order]
clone_sorted = cell_label[order]

clone_to_int = {cid: i for i, cid in enumerate(id_uniq)}
clone_int = np.array([clone_to_int[c] for c in clone_sorted])
clone_cmap = plt.get_cmap("tab10", len(id_uniq))

fig = plt.figure(figsize=(12, 10), dpi=300)
gs = fig.add_gridspec(
    nrows=2, ncols=2,
    height_ratios=[0.35, 8],
    width_ratios=[30, 1.2],
    hspace=0.05,
    wspace=0.10
)

ax_bar = fig.add_subplot(gs[0, 0])
ax_bar.imshow(clone_int[np.newaxis, :], aspect="auto", cmap=clone_cmap, interpolation="nearest")
ax_bar.set_axis_off()

ax = fig.add_subplot(gs[1, 0])
cmap_af = segpink.copy()
cmap_af.set_bad(color="white")

im = ax.imshow(AF_sorted, aspect="auto", cmap=cmap_af, interpolation="nearest", vmin=0.0, vmax=1.0)
ax.set_xlabel("Cells")
ax.set_ylabel("mtDNA variants")
ax.set_yticks(np.arange(len(mtSNP_ids)))
ax.set_yticklabels(mtSNP_ids, fontsize=6)

cax = fig.add_subplot(gs[1, 1])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Allele Frequency")

handles = [Patch(facecolor=clone_cmap(i), label=id_uniq[i]) for i in range(len(id_uniq))]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.18),
          ncol=min(len(id_uniq), 6), frameon=False)

plt.savefig(f"{VIREO_OUT}/AF_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("")
print("Completed heatmap.")


print("")
print("Saving output TSV files...")

# One-hot clone probabilities for compatibility
ID_prob = np.zeros((AD.shape[1], n_clones), dtype=float)
ID_prob[np.arange(AD.shape[1]), clone_ids] = 1.0
np.savetxt(f"{VIREO_OUT}/clone_id_prob.tsv", ID_prob, delimiter="\t")

# Re-emit clone assignments in the current barcode order
out_path = f"{VIREO_OUT}/clone_assignments.tsv"
with open(out_path, "w") as f:
    f.write("cell_barcode\tclone_id\n")
    for bc, cid in zip(cell_barcodes, clone_ids):
        f.write(f"{bc}\t{cid}\n")

# Mean allelic ratios per clone
np.savetxt(f"{VIREO_OUT}/clone_beta_mu_estimated_AF.tsv", beta_mu, delimiter="\t")

print("")
print("Pipeline completed successfully.")
