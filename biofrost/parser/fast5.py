import numpy as np
from ont_fast5_api.fast5_interface import get_fast5_file

__all__ = [
    "Fast5Data",
]


class Fast5Data(object):

    def __init__(self, workspace, seq_summary, batch_size=64, chunk_len=4_000, overlap=200, stride=10):
        self.workspace = workspace
        self.summary = seq_summary
        self.fast5_col = 'filename_fast5' if 'filename_fast5' in seq_summary.columns else "filename"
        # Iter to find all fast5s
        self.f5_files = {f5_file.name: f5_file for f5_file in workspace.rglob("*.fast5")}
        if len(self.f5_files) == 0:
            raise IndexError("No fast5 found")

        self.fasta = {}
        self.batch_size = batch_size
        self.chunk_len = chunk_len
        self.overlap = overlap
        self.stride = stride


    def __len__(self):
        return self.summary.shape[0]


    def get_read(self, read_id, f5=None, filter_outlier=False, trim_adapter=True):
        _f5 = None
        if f5 is None:
            row = self.summary.loc[read_id]
            f5_file = self.f5_files[row[self.fast5_col]]
            _f5 = get_fast5_file(f5_file, 'r')
            read = _f5.get_read(read_id)
        else:
            read = f5.get_read(read_id)

        # Determine analysis ID
        analysis_id = read.get_latest_analysis("Basecall_1D").split("_")[-1]
        move = read.get_analysis_dataset(f"Basecall_1D_{analysis_id}/BaseCalled_template", "Move")
        if move is None:
            raise KeyError(f"Move not found for read: {read_id}, please make sure the fast5s are generated using guppy --fast5_out")

        # Get Rawdata
        raw_data = read.get_raw_data()
        fastq = read.get_analysis_dataset(f'Basecall_1D_{analysis_id}/BaseCalled_template', 'Fastq').split('\n')[
                    1].replace("U", "T")[::-1]
        trace = read.get_analysis_dataset(f'Basecall_1D_{analysis_id}/BaseCalled_template', 'Trace')
        segments = read.get_analysis_attributes(f'Segmentation_{analysis_id}/Summary/segmentation')
        channel_info = read.get_channel_info()

        # Get scaling parameters
        raw_unit = channel_info['range'] / channel_info['digitisation']
        offset = channel_info['offset']
        current_level = np.array((raw_data + offset) * raw_unit, dtype=np.float32)
        if trim_adapter:
            current_level = current_level[segments['first_sample_template']:]

        # Filter outliers
        if filter_outlier:
            current_level = current_level[(current_level > 0) & (current_level < 200)]
            current_mean, current_stdv = np.mean(current_level), np.std(current_level)
            current_level = current_level[np.abs((current_level - current_mean) / current_stdv) <= 3]

        if _f5 is not None:
            _f5.close()

        # Return normalized level
        return current_level, fastq, move


    def yield_read(self):
        for fast5_name, df in self.summary.groupby(self.fast5_col):
            fast5_path = self.f5_files[fast5_name]
            read_names = list(df.index)

            with get_fast5_file(fast5_path, 'r') as f5:
                for read_id in read_names:
                    current_level, fastq, move = self.get_read(read_id, f5)
                    yield read_id, current_level, fastq

    def get_chunk(self):
        chunks, infos = [], []
        for read_id, current, fasta in self.yield_read():
            self.fasta[read_id] = fasta
            x_len = current.shape[0]

            breaks = list(np.arange(0, x_len, self.chunk_len-2*self.overlap))
            if breaks[-1] != x_len:
                breaks += [x_len, ]

            for start, end in zip(breaks[:-1], breaks[1:]):
                x_en = min(x_len, end+self.overlap)
                x_st = max(0, start-self.overlap)
                clip_st = int(start-x_st) // self.stride
                clip_en = int(x_en-end) // self.stride
                chunks.append(torch.from_numpy(scale_events(current[x_st:x_en])))
                infos.append((read_id, x_st, x_en, clip_st, int(x_en-x_st) // self.stride - clip_en))
                if len(chunks) >= self.batch_size:
                    chunks = pad_sequence(chunks)
                    yield infos, chunks.reshape(chunks.shape[0], chunks.shape[1], 1).permute(1, 0, 2)
                    chunks, infos = [], []
        if len(chunks) > 0:
            chunks = pad_sequence(chunks)
            yield infos, chunks.reshape(chunks.shape[0], chunks.shape[1], 1).permute(1, 0, 2)
