#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Generate six PDF figures from Proton.csv
#   • three for pair_type = Del
#   • three for pair_type = OS + SS (plotted together)
# Each figure has two panels (delta / gamma) and one curve per diff_bin.

import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------------------
# CONFIGURATION
# --------------------------------------------------
CSV_FILE   = "Proton.csv"        # input data
OUT_DIR    = "."                 # where PDFs are written
DIFF_BINS  = [0.5, 1.5, 2.5]     # bins to draw (order on the legend)
PAIRS_DEL  = ["Del"]             # Del alone
PAIRS_OSSS = ["OS", "SS"]        # OS + SS together
DIFF_TYPES = ["Intg", "DEta", "SPt"]

# Line/marker style per diff_bin
MARKER_STYLE = {0.5: "o-", 1.5: "s--", 2.5: "d-."}

# Colours for Del (one colour per diff_bin)
COLOUR_DEL = {0.5: "blue", 1.5: "orange", 2.5: "green"}

# Colours for OS / SS (one palette per pair, still separated by diff_bin)
COLOUR_OS  = {0.5: "teal",     1.5: "mediumturquoise", 2.5: "darkcyan"}
COLOUR_SS  = {0.5: "magenta",  1.5: "hotpink",         2.5: "purple"}

# --------------------------------------------------
# HELPER
# --------------------------------------------------
def plot_group(df, title, colour_map, outfile):
    """
    df          : DataFrame already filtered to a single (pair_type OR pair_group) & diff_type
    title       : title prefix for panels
    colour_map  : function (pair, bin) -> colour     (pair == "Del" or "OS"/"SS")
    outfile     : PDF file name
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

    for pair in sorted(df["pair_type"].unique()):
        for db in DIFF_BINS:
            sel = df[(df["pair_type"] == pair) & (df["diff_bin"] == db)]
            if sel.empty:
                continue
            x  = sel["centrality"]
            dy = sel["delta"]
            gy = sel["gamma"]
            dstat = sel["delta_err"]
            gstat = sel["gamma_err"]
            dsyst = sel["delta_syst_err"]
            gsyst = sel["gamma_syst_err"]

            colour = colour_map(pair, db)
            style  = MARKER_STYLE[db]

            # -------------------------------- delta
            axes[0].errorbar(x, dy, yerr=dstat,
                             fmt=style,  capsize=4,
                             color=colour,
                             label=f"{pair} bin={db}")
            axes[0].fill_between(x, dy-dsyst, dy+dsyst,
                                 color=colour, alpha=0.20)
            # -------------------------------- gamma
            axes[1].errorbar(x, gy, yerr=gstat,
                             fmt=style,  capsize=4,
                             color=colour,
                             label=f"{pair} bin={db}")
            axes[1].fill_between(x, gy-gsyst, gy+gsyst,
                                 color=colour, alpha=0.20)

    axes[0].set_title(f"{title}  –  delta")
    axes[0].set_ylabel("delta")
    axes[0].grid(True)
    axes[0].legend(fontsize="small", ncol=3)

    axes[1].set_title(f"{title}  –  gamma")
    axes[1].set_ylabel("gamma")
    axes[1].grid(True)
    axes[1].legend(fontsize="small", ncol=3)

    axes[1].set_xlabel("centrality")      # shared x-axis label

    plt.tight_layout()
    fig.savefig(f"{OUT_DIR}/{outfile}")
    plt.close(fig)

# --------------------------------------------------
# MAIN
# --------------------------------------------------
if __name__ == "__main__":
    # load & pre-filter data (drop centrality 65)
    data = pd.read_csv(CSV_FILE)
    data = data[data["centrality"] != 65.0]

    # ---------- Del plots (three files)
    for dtyp in DIFF_TYPES:
        d_del = data[(data["pair_type"] == "Del") & (data["diff_type"] == dtyp)]
        plot_group(
            d_del,
            title=f"Del  {dtyp}",
            colour_map=lambda _pair, db: COLOUR_DEL[db],
            outfile=f"Del_{dtyp}.pdf"
        )

    # ---------- OS + SS plots (three files)
    for dtyp in DIFF_TYPES:
        d_osss = data[(data["pair_type"].isin(PAIRS_OSSS)) &
                      (data["diff_type"] == dtyp)]
        plot_group(
            d_osss,
            title=f"OS + SS  {dtyp}",
            colour_map=lambda pair, db: (COLOUR_OS if pair == "OS" else COLOUR_SS)[db],
            outfile=f"OS_SS_{dtyp}.pdf"
        )

    print("✓ Six PDF files written.")