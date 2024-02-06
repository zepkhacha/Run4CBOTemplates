import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import gm2fr.src.io as io
import gm2fr.src.style as style
style.set_style()

import reference as ref

import argparse
import uproot
import itertools
import scipy.optimize as opt

# =================================================================================================

_std_vars = [
  "chi2",
  "chi2ndf"
]

for osc in ("CBO", "2CBO", "y", "VW", "VWmCBO"):
  _std_vars.append(f"f_{osc}")
  for band in ("", "_p", "_m"):
    for symbol in ("A", "phi", "alpha", "beta"):
      _std_vars.append(f"{symbol}{band}_{osc}")

# =================================================================================================
        
  def _format_fft():

    def label_freq(label, freq):
      plt.annotate(label, (freq, 0), (0, -3), textcoords = "offset pixels", rotation = "vertical", ha = "center", va = "top", fontsize = 6)

    for label, freq in ref.freq.items():
      label = ref.format_oscillation(label)
      label_freq(label, freq)
      if label != "a":
        label_freq(f"{label} + a", freq + ref.freq["a"])
        label_freq(f"{label} - a", freq - ref.freq["a"])

    style.draw_horizontal()
    ylow, yhigh = plt.ylim()
    plt.ylim(ylow - 0.2 * (yhigh - ylow), yhigh)
    plt.gca().set_xmargin(0)

# =================================================================================================

def plot(
  input_filename,
  output_filename,
  x_labels,
  y_labels,
  c_labels = None,
  tree = "tree",
  err_suffix = "_err"
):

  with uproot.open(input_filename) as file:
    file[tree].show()
    data = file[tree].arrays(library = "np")
    # data = file[tree].arrays(("fitStart", "alpha_CBO", "caloNum", "alpha_CBOerr"), library = "np")
  print(data)

  cut = (data['fitStart'] < 300)
  # cut = np.ones(len())
  x_labels = io.force_list(x_labels)
  y_labels = io.force_list(y_labels)
  c_labels = io.force_list(c_labels)

  with style.make_pdf(output_filename) as pdf:

    try:

      for c_label in c_labels:

        colors = None
        cbar_labels = None, None
        discrete_colors = False

        if c_label is not None:

          # !!!!!!!!!!!!!!!!!!!!!!!!
          c = data[c_label][cut]
          c_unique = np.unique(c)
          num_unique = len(c_unique)

          if len(c_unique) < 30 and c_label not in ("mid", "start", "end"):
            discrete_colors = True

            cmap = plt.get_cmap("coolwarm", num_unique)
            discrete_colors = [cmap(i) for i in range(num_unique)]
            cbar_labels = list(c_unique)
            if c_label == "calo":
              discrete_colors[0] = mpl.colors.to_rgba("k")
              cbar_labels[0] = "All"
            cmap = mpl.colors.LinearSegmentedColormap.from_list("", discrete_colors, N = len(discrete_colors))
            cmap_norm = mpl.colors.BoundaryNorm(np.arange(-0.5, num_unique), num_unique)

          else:

            cmap = plt.get_cmap("coolwarm")
            cmap_norm = mpl.colors.Normalize(np.min(c), np.max(c))
            # colors = cmap(cmap_norm(c))

          cmap_sm = plt.cm.ScalarMappable(cmap_norm, cmap)

        for x_label, y_label in itertools.product(x_labels, y_labels):

          print(f"Working on ({x_label}, {y_label}) plot.")

          x = data[x_label][cut]
          x_err = data[f"{x_label}{err_suffix}"][cut] if f"{x_label}{err_suffix}" in data else None

          y = data[y_label][cut]
          y_err = data[f"{y_label}{err_suffix}"][cut] if f"{y_label}{err_suffix}" in data else None

          # y_min, y_max = min(y), max(y)
          # x_min, x_max = min(x), max(x)

          # y_range = y_max - y_min
          # x_range = x_max - x_min

          # margin = 0.2
          # plt.xlim(x_min - margin * x_range, x_max + margin * x_range)
          # plt.ylim(y_min - margin * y_range, y_max + margin * y_range)

          def fit_double_exp(x, y, err):
            double_exp = lambda t, A, tau_A, B, tau_B, C: A*np.exp(-t/tau_A) + B*np.exp(-t/tau_B) + C
            try:
              p_opt, p_cov = opt.curve_fit(
                double_exp,
                x,
                y,
                p0 = [0, 50, 0, 10, 0], # one long lifetime, one short
                sigma = err,
                absolute_sigma = True,
                maxfev = 10000,
                bounds = (
                  [-np.inf, 0, -np.inf, 0, -np.inf],
                  [np.inf, np.inf, np.inf, np.inf, np.inf]
                )
              )
              plt.plot(x, double_exp(x, *p_opt), c = "k", alpha = 0.2, label = "Two-Exponential Fit")
            except:
              print("Failed double-exponential fit.")

          # TODO: add --line option to connect points, and --band option to draw bands instead of errorbars. independent options.
          # TODO: this logic is becoming messy, clean it up!
          # TODO: check if data type is a list of arrays, in which case plot each one as separate curve
          #           should automatically handle all such cases like fft, and also the original wiggle + fit, etc.
          if discrete_colors:
            for i, c_value in enumerate(c_unique):
              select = (c == c_value)
              x_select, y_select = np.squeeze(x[select]), np.squeeze(y[select])
              if y_label.endswith("fft") and c_label == "calo" and c_value == 0:
                y_select *= 1/np.sqrt(24)
              err_select = np.squeeze(y_err[select]) if y_err is not None else None
              zorder = 10 if (c_label == "calo" and i == 0) else 1
              # plt.plot(x_select, y_select, c = cmap(i), zorder = zorder, lw = 0.75)
              if err_select is not None:
                plt.fill_between(x_select, y_select - err_select, y_select + err_select, fc = cmap(i), alpha = 0.5, zorder = zorder - 1)
              if x_label == "mid" and y_label.startswith(("alpha", "beta", "A", "phi")):
                fit_double_exp(x_select, y_select, err_select)
          else:
            colors = cmap(cmap_norm(c)) if c_label is not None else None
            if not y_label.endswith("fft"):
              plt.scatter(x, y, c = colors)
              style.errorbar(x, y, y_err, x_err, ls = "", marker = "none", ecolor = colors, zorder = 0)
              if x_label == "mid" and y_label.startswith(("alpha", "beta", "A", "phi")):
                fit_double_exp(x, y, y_err)
            else:
              for i, c_value in enumerate(c):
                plt.plot(x[i], y[i], c = colors[i], lw = 0.5)

          if y_label.endswith("fft"):
            _format_fft()

          if c_label is not None:
            cbar = style.colorbar(mappable = cmap_sm, label = ref.format_label(c_label), ax = plt.gca())
            cbar.ax.tick_params(labelsize = 10)
            if discrete_colors:
              cbar.set_ticks(np.arange(num_unique))
              cbar.minorticks_off()
              if cbar_labels is not None:
                cbar.set_ticklabels(cbar_labels)

          # if calorimeter on the x or y axis, set the tick label for calo "0" as "All"
          for label, tick_format in [(x_label, plt.xticks), (y_label, plt.yticks)]:
            if label == "calo":
              ticks, labels = np.arange(0, 26, 5), list(np.arange(0, 26, 5))
              labels[0] = "All"
              tick_format(ticks, labels)

          for label, draw_line in [(x_label, style.draw_vertical), (y_label, style.draw_horizontal)]:
            if label.startswith(("alpha", "beta", "A", "phi")):
              draw_line(0)
            elif label == "chi2ndf":
              draw_line(1)

          # TODO: can also move textbox left/right and alter size, but hard to find best configuration
          # offset x axis label if there is a multiplier textbox
          # offset_limits = plt.rcParams["axes.formatter.limits"]
          # is_x_offset = (np.log10(abs(x_max if x_max != 0 else 1)) > offset_limits[1]) or (np.log10(abs(x_min if x_min != 0 else 1)) < offset_limits[0])

          # plt.xlim(right = 300)

          style.ylabel(ref.format_label(y_label))
          style.xlabel(ref.format_label(x_label))#, offset = 0.12 if is_x_offset else 0)
          style.make_unique_legend()
          pdf.savefig()
          plt.clf()

    except KeyboardInterrupt:
      print("Saving plots made so far and exiting...")

# =================================================================================================

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument("--input", required = True)
  parser.add_argument("--output", required = True)
  parser.add_argument("--tree", default = "fitresults")
  parser.add_argument("--err_suffix", default = "_err")
  parser.add_argument("-x", required = True, nargs = "+")
  parser.add_argument("-y", nargs = "+", default = _std_vars)
  parser.add_argument("-c", nargs = "*", default = None)
  args = parser.parse_args()

  plot(args.input, args.output, args.x, args.y, args.c, tree = args.tree, err_suffix = args.err_suffix)
