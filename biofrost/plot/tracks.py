import tempfile
import matplotlib.pyplot as plt
import pygenometracks.tracksClass
from pathlib import Path
from biofrost.utils import generate_random_key


class Section(object):
    def __init__(self, section_name, section_dict):
        self.name = section_name
        self.values = section_dict

    def __str__(self):
        lines = [f"[{self.name}]", ]
        for k, v in self.values.items():
            lines.append(f"{k} = {v}")
        return "\n".join(lines) + "\n"


class PyGenomeTracksIni(object):

    def __init__(self, filename=None):
        if filename is None:
            self.filename = Path(tempfile.gettempdir()) / f"pygenometracks_{generate_random_key(6)}.ini"
        else:
            self.filename = filename
        self.fh = open(self.filename, 'w')
        self.sections = []

    def add_section(self, section_name, section_dict):
        s = Section(section_name, section_dict)
        self.fh.write(str(s) + "\n")

    def plot(self, filename, contig, start, end, flank=500, title="", fig_width=15, fig_height=10):
        self.fh.close()

        p_start, p_end = start-flank, end+flank
        trp = pygenometracks.tracksClass.PlotTracks(
            tracks_file=self.filename,
            fig_width=fig_width,
            fig_height=fig_height,
            fontsize=7,
            dpi=300,
            track_label_width=0.1,
            plot_regions=[(contig, p_start, p_end)],
            plot_width=10,
        )
        fig = trp.plot(filename, contig, p_start, p_end, title)
        plt.show()


bed_template = {
    "file": "",
    "title": "",
    "height": 3,
    "color": "darkblue",
    "height": 3,
    "labels": "false",
    "fontsize": 7,
    "file_type": "bed",
}

bw_template = {
    "file": "",
    "title": "",
    "height": 1.5,
    "color": "#ef6548",
    "alpha": 1,
    "min_value": 0,
    "number_of_bins": 700,
    "nans_to_zeros": "true",
    "summary_method": "mean",
    "y_axis_values": "original",
    "file_type": "bigwig",
}

gtf_template = {
    "file": "",
    "titlte": "",
    "height": 3,
    "color": "darkblue",
    "height": 3,
    "labels": "false",
    "fontsize": 7,
    "file_type": "gtf",
}