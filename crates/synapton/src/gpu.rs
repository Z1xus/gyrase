use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, PoisonError, mpsc};
use std::time::Duration;

use wgpu::util::DeviceExt;

use crate::config::{AlignConfig, SeedConfig};
use crate::front::FrontBatch;
use crate::index::Index;
use crate::nucleotide::encode;

const MAX_READ: usize = 160;
const POOL_WORDS_PER_READ: usize = 128;
const ROW_BYTES: usize = 32;
const GROUP: usize = 64;
const UPLOAD_WORDS: usize = 1 << 20;
const TIMEOUT: Duration = Duration::from_secs(30);

#[derive(Debug, Default)]
struct Staging {
	codes: Vec<u8>,
	starts: Vec<u8>,
}

#[derive(Debug)]
pub(crate) struct Gpu {
	device: wgpu::Device,
	queue: wgpu::Queue,
	pipeline: wgpu::ComputePipeline,
	group: wgpu::BindGroup,
	params: wgpu::Buffer,
	bases: wgpu::Buffer,
	starts: wgpu::Buffer,
	table: wgpu::Buffer,
	pool: wgpu::Buffer,
	counter: wgpu::Buffer,
	table_out: wgpu::Buffer,
	counter_out: wgpu::Buffer,
	pool_out: wgpu::Buffer,
	bucket_bits: u32,
	cutoff: u32,
	capacity: usize,
	failed: Arc<AtomicBool>,
	staging: Mutex<Staging>,
}

fn words_buffer<T: Copy, const N: usize>(
	device: &wgpu::Device,
	values: &[T],
	bytes_of: impl Fn(T) -> [u8; N],
) -> wgpu::Buffer {
	let buffer = device.create_buffer(&wgpu::BufferDescriptor {
		label: None,
		size: (values.len() * N) as u64,
		usage: wgpu::BufferUsages::STORAGE,
		mapped_at_creation: true,
	});
	if let Ok(mut view) = buffer.slice(..).get_mapped_range_mut() {
		for (number, part) in values.chunks(UPLOAD_WORDS).enumerate() {
			let bytes: Vec<u8> = part.iter().flat_map(|&value| bytes_of(value)).collect();
			let start = number * UPLOAD_WORDS * N;
			view.slice(start..start + bytes.len())
				.copy_from_slice(&bytes);
		}
	}
	buffer.unmap();
	buffer
}

fn pipeline(device: &wgpu::Device) -> wgpu::ComputePipeline {
	let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
		label: None,
		source: wgpu::ShaderSource::Wgsl(include_str!("gpu.wgsl").into()),
	});
	device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
		label: None,
		layout: None,
		module: &module,
		entry_point: Some("main"),
		compilation_options: wgpu::PipelineCompilationOptions::default(),
		cache: None,
	})
}

fn buffer(device: &wgpu::Device, size: usize, usage: wgpu::BufferUsages) -> wgpu::Buffer {
	device.create_buffer(&wgpu::BufferDescriptor {
		label: None,
		size: size as u64,
		usage,
		mapped_at_creation: false,
	})
}

impl Gpu {
	pub(crate) fn new(index: &Index, config: &AlignConfig, capacity: usize) -> Option<Self> {
		if index.config() != SeedConfig::short_read() {
			return None;
		}
		let (bucket_bits, buckets, entries) = index.tables();
		let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
		let adapter = pollster::block_on(instance.enumerate_adapters(wgpu::Backends::PRIMARY))
			.into_iter()
			.filter(|adapter| {
				matches!(
					adapter.get_info().device_type,
					wgpu::DeviceType::DiscreteGpu | wgpu::DeviceType::IntegratedGpu
				) && adapter.features().contains(wgpu::Features::SHADER_INT64)
			})
			.min_by_key(|adapter| {
				adapter.get_info().device_type != wgpu::DeviceType::DiscreteGpu
			})?;
		let limits = adapter.limits();
		if entries.len() * 8 > usize::try_from(limits.max_storage_buffer_binding_size).ok()? {
			return None;
		}
		let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
			required_features: wgpu::Features::SHADER_INT64,
			required_limits: limits,
			..Default::default()
		}))
		.ok()?;
		let failed = Arc::new(AtomicBool::new(false));
		let flag = Arc::clone(&failed);
		device.on_uncaptured_error(Arc::new(move |_| flag.store(true, Ordering::Relaxed)));
		let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
		let pipeline = pipeline(&device);
		let storage = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST;
		let readable = wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC;
		let staging = wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST;
		let params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
			label: None,
			contents: &[0; 16],
			usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
		});
		let bases = buffer(&device, (capacity * MAX_READ).next_multiple_of(4), storage);
		let starts = buffer(&device, (capacity + 1) * 4, storage);
		let table = buffer(&device, capacity * ROW_BYTES, readable);
		let pool = buffer(&device, capacity * POOL_WORDS_PER_READ * 8, readable);
		let counter = buffer(&device, 4, readable | wgpu::BufferUsages::COPY_DST);
		let table_out = buffer(&device, capacity * ROW_BYTES, staging);
		let counter_out = buffer(&device, 4, staging);
		let pool_out = buffer(&device, capacity * POOL_WORDS_PER_READ * 8, staging);
		let bucket_data = words_buffer(&device, buckets, u32::to_le_bytes);
		let entry_data = words_buffer(&device, entries, u64::to_le_bytes);
		let group = device.create_bind_group(&wgpu::BindGroupDescriptor {
			label: None,
			layout: &pipeline.get_bind_group_layout(0),
			entries: &[
				params.as_entire_binding(),
				bases.as_entire_binding(),
				starts.as_entire_binding(),
				bucket_data.as_entire_binding(),
				entry_data.as_entire_binding(),
				table.as_entire_binding(),
				pool.as_entire_binding(),
				counter.as_entire_binding(),
			]
			.into_iter()
			.zip(0..)
			.map(|(resource, binding)| wgpu::BindGroupEntry { binding, resource })
			.collect::<Vec<_>>(),
		});
		if pollster::block_on(scope.pop()).is_some() || failed.load(Ordering::Relaxed) {
			return None;
		}
		Some(Self {
			device,
			queue,
			pipeline,
			group,
			params,
			bases,
			starts,
			table,
			pool,
			counter,
			table_out,
			counter_out,
			pool_out,
			bucket_bits,
			cutoff: u32::try_from(config.max_seed_hits.max(index.filter_cutoff())).ok()?,
			capacity,
			failed,
			staging: Mutex::new(Staging::default()),
		})
	}

	pub(crate) fn front(&self, sequences: &[&[u8]]) -> Option<FrontBatch> {
		let result = self.compute(sequences);
		if result.is_none() {
			self.failed.store(true, Ordering::Relaxed);
		}
		result
	}

	fn compute(&self, sequences: &[&[u8]]) -> Option<FrontBatch> {
		if self.failed.load(Ordering::Relaxed) || sequences.len() > self.capacity {
			return None;
		}
		let reads = sequences.len();
		let mut staging = self.staging.lock().unwrap_or_else(PoisonError::into_inner);
		let Staging { codes, starts } = &mut *staging;
		codes.clear();
		starts.clear();
		starts.extend_from_slice(&0u32.to_le_bytes());
		for sequence in sequences {
			if sequence.len() <= MAX_READ {
				encode(sequence, codes);
			}
			starts.extend_from_slice(&u32::try_from(codes.len()).ok()?.to_le_bytes());
		}
		codes.resize(codes.len().next_multiple_of(4), 0);
		let pool_words = u32::try_from(reads * POOL_WORDS_PER_READ).ok()?;
		let params = [
			u32::try_from(reads).ok()?,
			self.bucket_bits,
			self.cutoff,
			pool_words,
		];
		self.queue.write_buffer(
			&self.params,
			0,
			&params
				.iter()
				.flat_map(|value| value.to_le_bytes())
				.collect::<Vec<_>>(),
		);
		self.queue.write_buffer(&self.bases, 0, codes);
		self.queue.write_buffer(&self.starts, 0, starts);
		self.queue.write_buffer(&self.counter, 0, &[0; 4]);
		let table_bytes = (reads * ROW_BYTES) as u64;
		let mut encoder = self
			.device
			.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
		{
			let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
			pass.set_pipeline(&self.pipeline);
			pass.set_bind_group(0, &self.group, &[]);
			pass.dispatch_workgroups(u32::try_from(reads.div_ceil(GROUP)).ok()?, 1, 1);
		}
		encoder.copy_buffer_to_buffer(&self.table, 0, &self.table_out, 0, table_bytes);
		encoder.copy_buffer_to_buffer(&self.counter, 0, &self.counter_out, 0, 4);
		self.queue.submit([encoder.finish()]);
		let table = self.read(&self.table_out, table_bytes)?;
		let counter = self.read(&self.counter_out, 4)?;
		let used = u32::from_le_bytes(counter.as_chunks::<4>().0.first().copied()?).min(pool_words);
		let pool_bytes = u64::from(used) * 8;
		let pool = if pool_bytes == 0 {
			Vec::new()
		} else {
			let mut encoder = self
				.device
				.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
			encoder.copy_buffer_to_buffer(&self.pool, 0, &self.pool_out, 0, pool_bytes);
			self.queue.submit([encoder.finish()]);
			self.read(&self.pool_out, pool_bytes)?
		};
		if self.failed.load(Ordering::Relaxed) {
			return None;
		}
		let table = table
			.as_chunks::<4>()
			.0
			.iter()
			.map(|&word| u32::from_le_bytes(word))
			.collect();
		Some(FrontBatch::new(table, pool))
	}

	fn read(&self, source: &wgpu::Buffer, bytes: u64) -> Option<Vec<u8>> {
		let (sender, receiver) = mpsc::channel();
		source
			.slice(..bytes)
			.map_async(wgpu::MapMode::Read, move |result| {
				let _ = sender.send(result);
			});
		self.device
			.poll(wgpu::PollType::Wait {
				submission_index: None,
				timeout: Some(TIMEOUT),
			})
			.ok()?;
		receiver.recv().ok()?.ok()?;
		let data = source.slice(..bytes).get_mapped_range().ok()?.to_vec();
		source.unmap();
		Some(data)
	}
}
