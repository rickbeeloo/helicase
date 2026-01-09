use crate::{
    Config, FastxParser,
    input::{FromSlice, InputData, MmapInput},
};
use paraseq::ProcessError;
use rayon::ThreadPoolBuilder;
use rayon::iter::IntoParallelIterator;

use rayon::prelude::*;

pub struct ParallelHelicaseReader<'a, const CONFIG: Config> {
    mmap: MmapInput<'a>,
    batch_size: usize,
}

impl<'a, const CONFIG: Config> ParallelHelicaseReader<'a, CONFIG> {
    /// Batch size is the number of records between calls to `ParallelProcessor::on_batch_complete`.
    pub fn new(path: &std::path::Path, batch_size: usize) -> Self {
        Self {
            mmap: MmapInput::new(path).unwrap(),
            batch_size,
        }
    }
}

#[inline(always)]
fn find_fasta_boundary(data: &[u8], mut pos: usize) -> usize {
    while let Some(offset) = memchr::memchr(b'>', &data[pos..]) {
        let abs_pos = pos + offset;
        // Must be preceded by \n
        if abs_pos > 0 && data[abs_pos - 1] == b'\n' {
            return abs_pos;
        }
        pos = abs_pos + 1;
    }
    data.len()
}

#[inline(always)]
fn find_fastq_boundary(data: &[u8], mut pos: usize) -> usize {
    while let Some(offset) = memchr::memchr(b'\n', &data[pos..]) {
        let abs_pos = pos + offset;
        let next_pos = abs_pos + 1;

        // If we are at the end, this was quality line (continue)
        if next_pos >= data.len() || data[next_pos] != b'@' {
            pos = next_pos;
            continue;
        }

        // If we have a '@' at the next line, this was quality line (continue)
        let remaining = &data[next_pos..];
        if let Some(newline_offset) = memchr::memchr(b'\n', &remaining[1..]) {
            let next_line_start = newline_offset + 2;
            if next_line_start < remaining.len() && remaining[next_line_start] == b'@' {
                pos = next_pos;
                continue;
            }
        }

        // We made it here, so valid fastq returns positions
        return next_pos;
    }
    data.len()
}

impl<'b, const CONFIG: Config> paraseq::parallel::ParallelReader
    for ParallelHelicaseReader<'b, CONFIG>
{
    type Rf<'a> = &'a FastxParser<'a, CONFIG>;

    fn process_parallel<T>(self, processor: &mut T, num_threads: usize) -> paraseq::Result<()>
    where
        T: for<'a> paraseq::prelude::ParallelProcessor<Self::Rf<'a>> + Clone + Send,
    {
        let data = self.mmap.data();
        assert!(data[0] == b'>' || data[0] == b'@');
        let is_fastq = data[0] == b'@';

        let chunk_size = data.len().div_ceil(num_threads);

        let mut splits = (0..=num_threads)
            .map(|i| (i * chunk_size).min(data.len()))
            .collect::<Vec<usize>>();

        let len = splits.len();
        for split in &mut splits[1..len - 1] {
            *split = if is_fastq {
                find_fastq_boundary(data, *split)
            } else {
                find_fasta_boundary(data, *split)
            };
        }

        splits.dedup();

        // Could end up with less threads than splits
        let actual_threads = splits.len().saturating_sub(1);
        assert!(actual_threads > 0);

        let ranges: Vec<(usize, usize)> = splits.windows(2).map(|w| (w[0], w[1])).collect();

        let mut processors = Vec::with_capacity(actual_threads);
        for _ in 0..actual_threads {
            processors.push(processor.clone());
        }

        let pool = ThreadPoolBuilder::new()
            .num_threads(actual_threads)
            .build()
            .map_err(|_| ProcessError::Process("failed to build rayon pool".into()))?;

        pool.install(|| {
            // parallel iterate ranges zipped with processors
            ranges
                .into_par_iter()
                .zip(processors.into_par_iter())
                .map(|((start, end), mut worker_processor)| {
                    let slice = &data[start..end];
                    let mut parser = FastxParser::<CONFIG>::from_slice(slice);
                    let mut count = 0usize;

                    while let Some(_event) = parser.next() {
                        worker_processor.process_record(&parser)?;
                        count += 1;
                        if count == self.batch_size {
                            worker_processor.on_batch_complete()?;
                            count = 0;
                        }
                    }
                    if count > 0 {
                        worker_processor.on_batch_complete()?;
                    }
                    worker_processor.on_thread_complete()?;
                    Ok::<(), paraseq::ProcessError>(())
                })
                .collect::<Result<Vec<_>, _>>()
                .map(|_v| ())
        })
    }

    fn process_parallel_paired<T>(
        self,
        _r2: Self,
        _processor: &mut T,
        _num_threads: usize,
    ) -> paraseq::Result<()>
    where
        T: for<'a> paraseq::prelude::PairedParallelProcessor<Self::Rf<'a>>,
    {
        todo!()
    }

    fn process_parallel_interleaved<T>(
        self,
        _processor: &mut T,
        _num_threads: usize,
    ) -> paraseq::Result<()>
    where
        T: for<'a> paraseq::prelude::PairedParallelProcessor<Self::Rf<'a>>,
    {
        todo!()
    }

    fn process_parallel_multi<T>(
        self,
        _rest: Vec<Self>,
        _processor: &mut T,
        _num_threads: usize,
    ) -> paraseq::Result<()>
    where
        T: for<'a> paraseq::prelude::MultiParallelProcessor<Self::Rf<'a>>,
        Self: Sized,
    {
        todo!()
    }

    fn process_parallel_multi_interleaved<T>(
        self,
        _arity: usize,
        _processor: &mut T,
        _num_threads: usize,
    ) -> paraseq::Result<()>
    where
        T: for<'a> paraseq::prelude::MultiParallelProcessor<Self::Rf<'a>>,
    {
        todo!()
    }
}
